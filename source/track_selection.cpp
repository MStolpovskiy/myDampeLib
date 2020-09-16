#include "track_selection.hpp"
#include "etacorr.hpp"

#include <algorithm>
#include <math.h>


myDampeLib::DmpTrackSelector::DmpTrackSelector() :
    DmpTrackSelector("")
{;}

myDampeLib::DmpTrackSelector::DmpTrackSelector(const char * file)
{
    setBadChannelsFile(file);
    mSelectTypes.resize(0);

    mXYnonoverlaps = 2;
    mMaxDistStk2Bgo = 10.;
    mMaxAngStk2Bgo = 0.15;
    mSTKtrackChi2Max = 15.;
    mMinNPoints = 6;
    mSTKProtonMipE = 60.;
    mSTKXProtonMipRange = 0.1;
    mSTKYProtonMipRange = 0.1;
    mDistFromPrimCut = 0.;
}

myDampeLib::DmpTrackSelector::~DmpTrackSelector()
{;}

void myDampeLib::DmpTrackSelector::setBadChannelsFile(string file)
{
    mBadChannelsFile = file;

    mBadChannelList.resize(NLADDERS);
    for (vector<vector<bool> >::iterator it = mBadChannelList.begin();
         it != mBadChannelList.end();
         ++it) {
        it->resize(NCHANNELS);
    }
}

void myDampeLib::DmpTrackSelector::setSelectTypes(vector<Select> types)
{
    for (vector<Select>::iterator it = types.begin(); it != types.end(); it++) {
        addSelect(*it);
    }
}

void myDampeLib::DmpTrackSelector::addSelect(Select type){
    vector<Select>::iterator it;
    it = find(mSelectTypes.begin(), mSelectTypes.end(), type);
    if (it == mSelectTypes.end()) {   
        switch (type) {
            case stk_bad_channel :
                readBadChannelsFile();
                break;
            default :
                break;
        }
        mSelectTypes.push_back(type);
    }
}

void myDampeLib::DmpTrackSelector::removeSelect(Select type) {
    vector<Select>::iterator it;
    it = find(mSelectTypes.begin(), mSelectTypes.end(), type);
    if (it != mSelectTypes.end()) {
        mSelectTypes.erase(it);
    }
}

bool myDampeLib::DmpTrackSelector::selected(DmpStkTrack* track, DmpEvent * event) const {
    for (vector<Select>::const_iterator it = mSelectTypes.begin();
         it != mSelectTypes.end(); 
         ++it) {
        if(!pass(track, event, *it)) return false;
    }
    return true;
}

bool myDampeLib::DmpTrackSelector::pass(DmpStkTrack * track, DmpEvent * event, Select type) const {
    switch (type) {
        case stk_bad_channel :
            return (hasBadChannel(track, event) != -1);
            break;
        case psd_match :
            return psdMatch(track, event);
        case stk_xy_overlap :
	    return stkXYoverlap(track);
        case stk_missing_impact_point :
	    return stkNoMissingImpactPoint(track);
        case stk_impact_first_layer :
	    return stkImpactFirstLayer(track);
        case stk_bgo_dist_high :
            return stk2bgoCut(track, event);
        case stk_chi2_cut :
	    return stkChi2Cut(track);
        case stk_npoints :
	    return stkNPoints(track);
        case stk_nXYhits :
	    return stkNXYhits(track);
        case stk_no_1strip_clusters :
	    return stkNo1stripClusters(track, event);
        case proton_mip :
	    return isProtonMip(track, event);
        case match_primary :
            return dist_from_prim_cut(track, event);
        default :
            return false;
    }
}

bool myDampeLib::DmpTrackSelector::readBadChannelsFile() {
    ifstream infile(mBadChannelsFile.c_str());  
    if (!infile)
    {
        std::cout << "Can't open Bad Channels file: " << mBadChannelsFile;
        std::cout << " ==> throwing exception!" << std::endl;
        return false;
    }

    std::string line;
    while (getline(infile, line)) {
        if (line.c_str() == std::string("")) break;
        std::istringstream lineStream(line.c_str());
        std::string number;
        int i = 0;
        int ildr   = -1;
        int bdchnl = -1;
        while (getline(lineStream, number, ',')) {
            switch (i++) {
            case 0:
                ildr = atoi(number.c_str());
                break;
            case 1:
                bdchnl = atoi(number.c_str());
                break;          
            }
        }
        if(ildr>=0 && bdchnl>=0 && ildr<NLADDERS && bdchnl<NCHANNELS){
            mBadChannelList[ildr][bdchnl] = true;
        }
    }
    return true;
}

int myDampeLib::DmpTrackSelector::hasBadChannel(DmpStkTrack * track, DmpEvent * pev) const{
    vector<int> layers = layerBadChannel(track, pev);
    if (layers.size() == 0) return -1;
    else return layers[0];
}

vector<int> myDampeLib::DmpTrackSelector::layerBadChannel(DmpStkTrack * track, DmpEvent * pev) const {
    TClonesArray * stkclusters = pev->GetStkSiClusterCollection();
    vector <int> ret(0);

    // Loop over clusters
    for (int ipoint=0; ipoint<track->GetNPoints(); ipoint++) {
        for (int ixy=0; ixy<2; ixy++) {
            DmpStkSiCluster* cluster;
	    if(ixy == 0)
	        cluster = track -> GetClusterX(ipoint, stkclusters);
            else
	        cluster = track -> GetClusterY(ipoint, stkclusters);
	    if(!cluster) continue;

	    //ladder num
	    int ladder = cluster->getLadderHardware();

	    //address of the first strip in the cluster (0-363)
	    int minc = cluster->getIndex1();

	    //last strip of the cluster
	    int maxc = cluster->getIndex1() + cluster->getNstrip() -1;

//          // Check if a bad cluster with 0 signal is on the edge of the cluster
//	    if (minc != 0 && mBadChannelList[ladder][minc - 1])
//	        ret.push_back(ipoint);
//	    if (maxc < NCHANNELS && mBadChannelList[ladder][maxc + 1])
//		ret.push_back(ipoint);

            // Check if a bad cluster is inside the cluster
	    for(int i = minc; i <= maxc; i++){
	        if(mBadChannelList[ladder][i]) {
		    ret.push_back(ipoint);
		    break;
		}
	    }
        }
    }

    return ret;
}

bool myDampeLib::DmpTrackSelector::psdMatch(DmpStkTrack * track, DmpEvent * event) const {
    DmpEvtPsdRec *psdRec = event->pEvtPsdRec();

    TVector3 impact = track->getImpactPoint();
    TVector3 direction = track->getDirection();

    bool mc = (event->pEvtSimuHeader() != nullptr);

    // Loop over psd layers
    for (int ilayer = 0; ilayer < NPSDLAYERS; ilayer++) {
        // Loop over psd bars
        for (int ibar = 0; ibar < NPSDBARS; ibar++) {
            double e = psdRec->GetEdep(ilayer, ibar);
            if (e == 0) continue;
            double * len = new double[2];
	    bool check_fd = !mc && gPsdECor->GetPathLengthPosition(ilayer, ibar, direction, impact, len);
	    bool check_mc = mc && gPsdECor->GetPathLPMC(ilayer, ibar, direction, impact, len);
	    if (check_fd || check_mc) return true;
        }
    }
    return false;
}

bool myDampeLib::DmpTrackSelector::stkXYoverlap(DmpStkTrack * track) const {
    int nonoverlaps = track->getNhitX() + track->getNhitY() - 2 * track->getNhitXY();
    return nonoverlaps <= mXYnonoverlaps;
}

bool myDampeLib::DmpTrackSelector::stkNoMissingImpactPoint(DmpStkTrack * track) const {
    bool ret = true;
    bool hasX = track->getImpactPointHasX();
    bool hasY = track->getImpactPointHasY();
    if (mNimpactPoints == 1) ret = hasX || hasY;
    if (mNimpactPoints == 2) ret = hasX && hasY;
    return ret;
}

bool myDampeLib::DmpTrackSelector::stkImpactFirstLayer(DmpStkTrack * track) const {
    if (track->getImpactPointPlane() > 0) return false;
}

float myDampeLib::DmpTrackSelector::stk2bgoAngl(DmpStkTrack * track, DmpEvent * event) const {
    TVector3 trackdirection = track->getDirection();
    DmpEvtBgoRec* bgoRec = event->pEvtBgoRec();
    if( (bgoRec->GetSlopeXZ()==0 && bgoRec->GetInterceptXZ()==0) ||
        (bgoRec->GetSlopeYZ()==0 && bgoRec->GetInterceptYZ()==0)) return false;

    TVector3 bgodirection = TVector3(bgoRec->GetSlopeXZ(),
                                     bgoRec->GetSlopeYZ(),
                                     1);

    float angle = trackdirection.Angle(bgodirection);
    return angle;
}

float myDampeLib::DmpTrackSelector::stk2bgoDist(DmpStkTrack * track, DmpEvent * event) const {
    DmpEvtBgoRec* bgoRec = event->pEvtBgoRec();

    float x = track->getImpactPoint().x();
    float y = track->getImpactPoint().y();
    float z = track->getImpactPoint().z();

    float px = track->getDirection().x();
    float py = track->getDirection().y();
    float pz = track->getDirection().z();

    double bgoRec_slope[2];
    double bgoRec_intercept[2];
    double BGO_TopZ = 46;

    bgoRec_slope[1] = bgoRec->GetSlopeXZ();
    bgoRec_slope[0] = bgoRec->GetSlopeYZ();
    bgoRec_intercept[1] = bgoRec->GetInterceptXZ();
    bgoRec_intercept[0] = bgoRec->GetInterceptYZ();

    double topX = bgoRec_slope[1]*BGO_TopZ + bgoRec_intercept[1];
    double topY = bgoRec_slope[0]*BGO_TopZ + bgoRec_intercept[0];
    TVector3 bgoImpact(topX, topY, BGO_TopZ);

    float x_bgo = x + (bgoImpact.z() - z) * px/pz;
    float y_bgo = y + (bgoImpact.z() - z) * py/pz;

    float dist = pow(pow((bgoImpact.x() - x_bgo), 2) +
                     pow((bgoImpact.y() - y_bgo), 2), 0.5);
    return dist;
}

bool myDampeLib::DmpTrackSelector::stk2bgoCut(DmpStkTrack * track, DmpEvent * event) const {
    // Check angle match
    float angle = stk2bgoAngl(track, event);
    bool angle_ok = (angle <= mMaxAngStk2Bgo);

    // Check distance match
    float dist = stk2bgoDist(track, event);
    bool dist_ok = (dist <= mMaxDistStk2Bgo);

    return angle_ok && dist_ok;
}

bool myDampeLib::DmpTrackSelector::stkChi2Cut(DmpStkTrack * track) const {
    double chisqndof = track -> getChi2NDOF();
    return chisqndof <= mSTKtrackChi2Max;
}

bool myDampeLib::DmpTrackSelector::stkNPoints(DmpStkTrack * track) const {
    return track->GetNPoints() >= mMinNPoints;
}

bool myDampeLib::DmpTrackSelector::stkNXYhits(DmpStkTrack * track) const {
    return track->getNhitXY() >= mMinNXYhits;
}

bool myDampeLib::DmpTrackSelector::stkNo1stripClusters(DmpStkTrack * track, DmpEvent * event) const {
    TClonesArray * stkclusters = event->GetStkSiClusterCollection();

    for (int ixy=0; ixy<2; ixy++) {
        for (int ipoint=0; ipoint<track->GetNPoints(); ipoint++) {
            DmpStkSiCluster* cluster;
            if(ixy == 0)
                cluster = track -> GetClusterX(ipoint, stkclusters);
            else
                cluster = track -> GetClusterY(ipoint, stkclusters);
            if(!cluster) continue;
	    if(cluster->getNstrip() == 1) return false;
        }
    }

    return true;
}

void myDampeLib::DmpTrackSelector::getSMean(DmpStkTrack * track, DmpEvent * event, float * s_mean_xy) const{
    s_mean_xy[0] = 0.;
    s_mean_xy[1] = 0.;
    int n_points_xy[] = {0, 0};
    TClonesArray * stkclusters = event->GetStkSiClusterCollection();

    double costheta = track->getDirection().CosTheta();

    for (int ixy=0; ixy<2; ixy++) {
	for (int ipoint=0; ipoint<track->GetNPoints(); ipoint++) {
	    DmpStkSiCluster* cluster;
	    if(ixy == 0)
		cluster = track -> GetClusterX(ipoint, stkclusters);
	    else
		cluster = track -> GetClusterY(ipoint, stkclusters);
	    if(!cluster) continue;
	    if (cluster->getNstrip() == 1) continue;

	    double e = cluster -> getEnergy() * costheta;

	    s_mean_xy[ixy] += e;
	    n_points_xy[ixy]++;
	}

	s_mean_xy[ixy] *= 1. / float(n_points_xy[ixy]) / mSTKProtonMipE;
	s_mean_xy[ixy] = sqrt(s_mean_xy[ixy]);
    }
}

bool myDampeLib::DmpTrackSelector::isProtonMip(DmpStkTrack * track, DmpEvent * event) const{
    float s_mean_xy[2];
    getSMean(track, event, s_mean_xy);

    bool ret = true;

    float min[] = {1. - mSTKXProtonMipRange, 1. - mSTKYProtonMipRange};
    float max[] = {1. + mSTKXProtonMipRange, 1. + mSTKYProtonMipRange};

    for (int ixy=0; ixy<2; ixy++) {
	if (s_mean_xy[ixy] < min[ixy] || s_mean_xy[ixy] > max[ixy]) ret = false;
    }
    return ret;
}

float myDampeLib::DmpTrackSelector::charge_consistency(DmpStkTrack * track, DmpEvent * pev, int charge) const {
    // Returns the sum of (e - charge**2)**2 
    // where e are the charge measurements from PSD and first two layers of STK
    // and charge is for example 1 for proton and 2 for helium

    float e[6] = {0., 0., 0., 0., 0., 0.};

    // PSD charge
    int * igbar = new int(-1);
    bool mc = (pev->pEvtSimuHeader() != nullptr);
    double psde = psdEnergy(track, pev->pEvtPsdHits(), igbar, mc, 1);
    e[0] = psde / 2.07;
    psde = psdEnergy(track, pev->pEvtPsdHits(), igbar, mc, 2);
    e[1] = psde / 2.07;
    
    // STK charge
    EtaCorr etaCorr;
    etaCorr.setTargetParameters(64.27, 6., 257.75, 6.*257.75/64.27);
    TClonesArray * stkclusters = pev->GetStkSiClusterCollection();
    double costheta = track->getDirection().CosTheta();
    for (int ipoint=0; ipoint<track->GetNPoints(); ipoint++) {
	for (int ixy=0; ixy<2; ixy++) {
	    // Check if cluster is associated to a hit                                                                       
	    DmpStkSiCluster* cluster;
	    if(ixy == 0)
		cluster = track -> GetClusterX(ipoint, stkclusters);
	    else
		cluster = track -> GetClusterY(ipoint, stkclusters);
	    if(!cluster) continue;
	    double estk = cluster -> getEnergy() * costheta;
	    double eta = etaCorr.calcEta(cluster);
	    int index = ixy + (cluster -> getPlane()) * 2;
	    if (index < 4) {
		e[index + 2] = etaCorr.corrEnergy(estk, eta, costheta) / 60.; // 60 is the proton mip position in STK
	    }
	}
    }

    float cc = 0.;
    for(int i = 0; i < 6; i++) {
	if (e[i] != 0.) {
	    cc += pow(e[i] - charge*charge, 2);
	}
    }

    return cc;
}

float myDampeLib::DmpTrackSelector::dist_from_prim(DmpStkTrack * track, DmpEvent * pev) const {
    // Reconstruction of the primary track based on the linear model of :
    // - charge consistency of the track (only proton and helium!)
    // - inpact plane
    // - chi2

    // Charge consistency for the proton
    float ccp = charge_consistency(track, pev, 1);

    // Charge consistency for the helium
    float cch = charge_consistency(track, pev, 2);

    // Impact plane (all harcoded values are from the training I did in python)
    float ip = track->getImpactPointPlane() == 0 ? -1.31115115 : 1.21147072;

    // chi2
    double chi2 = track -> getChi2NDOF();

    // Coefs of the linear (Ridge) model:
    // -2.72187012,  2.8741397 ,  0.83496084,  0.12752761
    float dist = -2.72187012 * log10(ccp) + 2.8741397 * log10(cch) + 0.83496084 * ip + 0.12752761 * chi2;

    return dist;
}

bool myDampeLib::DmpTrackSelector::dist_from_prim_cut(DmpStkTrack * track, DmpEvent * pev) const {
    float dist = dist_from_prim(track, pev);
    return dist < mDistFromPrimCut;
}

bool myDampeLib::DmpTrackSelector::first_is_better(DmpStkTrack * lhs, DmpStkTrack * rhs, DmpEvent * pev) const {
    // Ternary comparison used.
    // -1 means lhs track is better
    // 0  -- tracks are equal
    // +1 means rhs track is better
    // Final result is the the sum of all the comparisons. If it is negative, the lhs is better.

    // Bad channel comparison
    short bad_chan_comp = 0;
    vector<Select>::const_iterator it;
    it = find(mSelectTypes.begin(), mSelectTypes.end(), stk_bad_channel);
    if (it != mSelectTypes.end()) {
	vector<int> lhs_bad_chan_layers = layerBadChannel(lhs, pev);
	vector<int> rhs_bad_chan_layers = layerBadChannel(rhs, pev);
        bad_chan_comp = lhs_bad_chan_layers.size() - rhs_bad_chan_layers.size();
        if(bad_chan_comp != 0) bad_chan_comp /= bad_chan_comp;
    }

    // STK XY non overlaps
    int lhs_nonoverlaps = lhs->getNhitX() + lhs->getNhitY() - 2 * lhs->getNhitXY();
    int rhs_nonoverlaps = rhs->getNhitX() + rhs->getNhitY() - 2 * rhs->getNhitXY();
    short nonoverlaps_comp = lhs_nonoverlaps - rhs_nonoverlaps;
    if (nonoverlaps_comp != 0) nonoverlaps_comp /= nonoverlaps_comp;

    // Missing impact point
    short missing_impact_comp = 0;
    it = find(mSelectTypes.begin(), mSelectTypes.end(), stk_missing_impact_point);
    if (it == mSelectTypes.end()) {
	missing_impact_comp = -(short)stkNoMissingImpactPoint(lhs) + (short)stkNoMissingImpactPoint(rhs);
	if(missing_impact_comp != 0) missing_impact_comp /= missing_impact_comp;
    }

    // BGO - STK angle
    float lhs_angle = stk2bgoAngl(lhs, pev);
    float rhs_angle = stk2bgoAngl(rhs, pev);
    short bgo_stk_angle_comp = (lhs_angle > rhs_angle) ? 1 : -1;
    if (lhs_angle == rhs_angle) bgo_stk_angle_comp = 0;

    // BGO - STK distance
    float lhs_dist = stk2bgoDist(lhs, pev);
    float rhs_dist = stk2bgoDist(rhs, pev);
    short bgo_stk_dist_comp = (lhs_dist > rhs_dist) ? 1 : -1;
    if (lhs_dist == rhs_dist) bgo_stk_dist_comp = 0;

    // Chi2 comparison
    double lhs_chisqndof = lhs -> getChi2NDOF();
    double rhs_chisqndof = rhs -> getChi2NDOF();
    short chi2_comp = (lhs_chisqndof > rhs_chisqndof) ? 1 : -1;
    if (lhs_chisqndof == rhs_chisqndof) chi2_comp = 0;

    // N points comparison
    short npoints_comp = - lhs->GetNPoints() + rhs->GetNPoints();
    if(npoints_comp != 0) npoints_comp /= npoints_comp;

    // Proton MIP comparison
    short proton_mip_comp = 0;
    it = find(mSelectTypes.begin(), mSelectTypes.end(), proton_mip);
    if (it != mSelectTypes.end()) {
	float lhs_s_mean_xy[2], rhs_s_mean_xy[2];
	getSMean(lhs, pev, lhs_s_mean_xy);
	getSMean(rhs, pev, rhs_s_mean_xy);

        float lhs_s_mean = sqrt(pow(lhs_s_mean_xy[0] - 1, 2) + pow(lhs_s_mean_xy[1] - 1, 2));
	float rhs_s_mean = sqrt(pow(rhs_s_mean_xy[0] - 1, 2) + pow(rhs_s_mean_xy[1] - 1, 2));

	proton_mip_comp = (lhs_s_mean > rhs_s_mean) ? 1 : -1;
	if (lhs_s_mean == rhs_s_mean) proton_mip_comp = 0;
    }

    short ret = bad_chan_comp + nonoverlaps_comp + 
	        missing_impact_comp + bgo_stk_angle_comp + bgo_stk_dist_comp +
                chi2_comp + npoints_comp + proton_mip_comp;

    return ret < 0;
}

bool myDampeLib::DmpTrackSelector::first_is_closer(DmpStkTrack * lhs, DmpStkTrack * rhs, DmpEvent * pev) const {
    // Track corresponds to the primary particle                                                                                
    float lscore = dist_from_prim(lhs, pev);
    float rscore = dist_from_prim(rhs, pev);

    return lscore < rscore;
}
