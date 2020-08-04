#include "event_selection.hpp"
#include "definitions.hpp"
#include <algorithm>
#include <DmpEvtPsdHits.h>

myDampeLib::DmpEventSelector::DmpEventSelector(bool check_all/*=false*/)
{
    mSelectTypes.resize(0);
    mHselect = new TH1I("Hselect", "Selection criteria", SELECT_N_ITEMS + 1, 0, SELECT_N_ITEMS + 1);
    mHselect->GetXaxis()->SetBinLabel(1, "All");
    for (int i=1; i <= SELECT_N_ITEMS; i++) {
        Select type = (Select)(i-1);
	string type_str;
	switch (type) {
            case bgo_skim_cuts : type_str = "Skim"; break;
            case het : type_str = "HET"; break;
            case let : type_str = "LET"; break;
            case has_STK_track : type_str = "has STK track"; break;
            case has_only_one_STK_track : type_str = "has only one STK track"; break;
            case bgo_psd_match : type_str= "BGO to PSD match"; break;
            default: type_str = "";
	}
        mHselect->GetXaxis()->SetBinLabel(i+1, type_str.c_str());
    }
    if (check_all) checkAll();
}

myDampeLib::DmpEventSelector::~DmpEventSelector()
{
    delete mTrackSelector;
}

void myDampeLib::DmpEventSelector::setSelectTypes(vector<Select> types)
{
    for (vector<Select>::iterator it = types.begin(); it != types.end(); it++) {
        addSelect(*it);
    }
}

void myDampeLib::DmpEventSelector::addSelect(Select type)
{
    vector<Select>::iterator it;
    it = find(mSelectTypes.begin(), mSelectTypes.end(), type);
    if (it == mSelectTypes.end())
        mSelectTypes.push_back(type);
}

void myDampeLib::DmpEventSelector::removeSelect(Select type)
{
    vector<Select>::iterator it;
    it = find(mSelectTypes.begin(), mSelectTypes.end(), type);
    if (it != mSelectTypes.end())
        mSelectTypes.erase(it);
}

bool myDampeLib::DmpEventSelector::selected(DmpEvent * event)
{
    bool ret;
    mHselect->Fill(0);
    for (vector<Select>::const_iterator it = mSelectTypes.begin();
         it != mSelectTypes.end(); 
         ++it) {
        ret = pass(event, *it);
	if (ret) mHselect->Fill(*it + 1);
	else break;
    }
    return ret;
}

bool myDampeLib::DmpEventSelector::pass(DmpEvent * event, Select type)
{
    bool ret = true;
    switch (type) {
        case het :
            ret = event->pEvtHeader()->GeneratedTrigger(3); break;
        case let :
            ret = event->pEvtHeader()->GeneratedTrigger(2); break;
        case has_STK_track :
            ret = hasSTKtrack(event); break;
        case has_only_one_STK_track:
	    ret = hasOnlyOneSTKtrack(event); break;
        case bgo_skim_cuts :
            ret = notSideIn(event) && bgoContainment(event) && bgoMaxBar(event); break;
        case bgo_psd_match :
	    ret = bgoPsdMatch(event); break;
        default :
            ret = false;
    }
    return ret;
}

bool myDampeLib::DmpEventSelector::hasSTKtrack(DmpEvent * pev)
{
    TClonesArray * tracks = pev->GetStkKalmanTrackCollection();
    if (pev->NStkKalmanTrack () < 1) return false;

    mNtracks = 0;
    int best_track_index = 0;

    for(int i=0; i < pev->NStkKalmanTrack(); i++)
    {
	DmpStkTrack * track = (DmpStkTrack*)tracks->ConstructedAt(i);
	if (mTrackSelector->selected(track, pev)) {
	    if (mNtracks == 0) {
		best_track_index = i;
	    }
	    else 
	    {
		DmpStkTrack * best_track = (DmpStkTrack*)tracks->ConstructedAt(best_track_index);
		if(mTrackSelector->first_is_better(track, best_track, pev)) {
		    best_track_index = i;
		}
	    }
	    mNtracks++;
	}
    }
    if(mNtracks > 0) mSTKtrack = (DmpStkTrack*)tracks->ConstructedAt(best_track_index);
    return mNtracks > 0;
}

bool myDampeLib::DmpEventSelector::hasOnlyOneSTKtrack(DmpEvent * pev)
{
    vector<Select>::iterator it;
    it = find(mSelectTypes.begin(), mSelectTypes.end(), has_STK_track);
    if (it == mSelectTypes.end())
	hasSTKtrack(pev); // calculate mNtracks
    return mNtracks == 1;
}

bool myDampeLib::DmpEventSelector::notSideIn(DmpEvent * pev) const
{
    DmpEvtBgoRec* bgoRec = pev->pEvtBgoRec();
    double bgo_layer_sum[NLAYERSTOT];
    for (int ilayer=0; ilayer<NLAYERSXY; ilayer++){
        bgo_layer_sum[ilayer * 2] = 0.;
        bgo_layer_sum[ilayer * 2 + 1] = 0.;
        for (int ibar=0; ibar<NBARSL; ibar++) {
            bgo_layer_sum[ilayer * 2] += bgoRec->GetEdep(ilayer * 2, ibar);
            bgo_layer_sum[ilayer * 2 + 1] += bgoRec->GetEdep(ilayer * 2 + 1, ibar);
        }
    }

    // Energy_L2 + Energy_L3 > Energy_L0 + Energy_L1 To exclude Down-to-Up events
    if (bgo_layer_sum[0] + bgo_layer_sum[1] >
        bgo_layer_sum[2] + bgo_layer_sum[3]) return false;

    // MaxLayerEnergyRatio < 0.35 To exclude Side-In events
    double max_layer_ratio = 0.;
    double bgo_energy = bgoRec->GetTotalEnergy();
    for (int ilayer=0; ilayer<NLAYERSTOT; ilayer++) {
        double r = bgo_layer_sum[ilayer] / bgo_energy;
        if (max_layer_ratio < r) max_layer_ratio = r;
    }
    if (max_layer_ratio >= mMaxLayerRatio) return false;

    return true;
}

bool myDampeLib::DmpEventSelector::bgoSlopeIntercept(DmpEvent * pev) const
{
    DmpEvtBgoRec* bgorec = pev->pEvtBgoRec();
    if( (bgorec->GetSlopeXZ()==0 && bgorec->GetInterceptXZ()==0) ||
	(bgorec->GetSlopeYZ()==0 && bgorec->GetInterceptYZ()==0)) return false;
    return true;
}

bool myDampeLib::DmpEventSelector::bgoContainment(DmpEvent * event) const
{
    DmpEvtBgoRec* bgorec = event->pEvtBgoRec();
    double BGO_TopZ = 46;
    double BGO_BottomZ = 448;
    double bgoRec_slope[2];
    double bgoRec_intercept[2];
    if(!bgoSlopeIntercept(event)) return false;
    bgoRec_slope[1] = bgorec->GetSlopeXZ();
    bgoRec_slope[0] = bgorec->GetSlopeYZ();
    bgoRec_intercept[1] = bgorec->GetInterceptXZ();
    bgoRec_intercept[0] = bgorec->GetInterceptYZ();

    double topX = bgoRec_slope[1]*BGO_TopZ + bgoRec_intercept[1];
    double topY = bgoRec_slope[0]*BGO_TopZ + bgoRec_intercept[0];
    double bottomX = bgoRec_slope[1]*BGO_BottomZ + bgoRec_intercept[1];
    double bottomY = bgoRec_slope[0]*BGO_BottomZ + bgoRec_intercept[0];
    if(fabs(topX)<280 && 
       fabs(topY)<280 && 
       fabs(bottomX)<280 && 
       fabs(bottomY)<280) return true;
    return false;
}

bool myDampeLib::DmpEventSelector::bgoMaxBar(DmpEvent * event) const
{
    DmpEvtBgoHits * bgohits = event -> pEvtBgoHits ();
    int nBgoHits = bgohits->GetHittedBarNumber();
    short  barNumberMaxEBarLay1_2_3[3] = {-1}; //bar number of maxE bar in layer 1, 2, 3                                                
    double MaxEBarLay1_2_3[3] = {0};           //E of maxE bar in layer 1, 2, 3                                                         
    for (int ihit = 0; ihit <nBgoHits; ihit++){
	double hitE = (bgohits->fEnergy)[ihit];
	short   lay = bgohits->GetLayerID(ihit);
	if(lay==1 || lay==2 || lay==3) {
	    if(hitE > MaxEBarLay1_2_3[lay-1]) {
		int iBar = ((bgohits -> fGlobalBarID)[ihit] >> 6) & 0x1f;
		MaxEBarLay1_2_3[lay-1] = hitE;
		barNumberMaxEBarLay1_2_3[lay-1] = iBar;
	    }
	}
    }
    for (int j = 0; j <3; j++){
	if(barNumberMaxEBarLay1_2_3[j] <=0 || barNumberMaxEBarLay1_2_3[j] == 21) return false;
    }
    return true;
}

bool myDampeLib::DmpEventSelector::bgoPsdMatch(DmpEvent * event) const {
    DmpEvtBgoRec* bgoRec = event->pEvtBgoRec();

    float psdXZ = -291.5;
    float psdYZ = -317.7;

    double bgoXproj = bgoRec->GetSlopeXZ()*psdXZ + bgoRec->GetInterceptXZ();
    double bgoYproj = bgoRec->GetSlopeYZ()*psdYZ + bgoRec->GetInterceptYZ();

    // DmpEvtPsdRec *psdRec = event->pEvtPsdRec();
    DmpEvtPsdHits * psdHits = event -> pSimuPsdHits ();
    bool matchx = false;
    bool matchy = false;

    // Loop over psd hits
    for (int ihit = 0; ihit < event -> NEvtPsdHits(); ihit++ ) {
	if (!matchx && psdHits->IsHitMeasuringX (ihit)) {
	    double x = psdHits->GetHitX(ihit);
	    matchx = abs(bgoXproj - x) < mBgoPsdDist;
	}
	if (!matchy && psdHits->IsHitMeasuringY (ihit)) {
	    double y = psdHits->GetHitY(ihit);
	    matchy = abs(bgoYproj - y) < mBgoPsdDist;
	}
    }
    return matchx && matchy;
}
