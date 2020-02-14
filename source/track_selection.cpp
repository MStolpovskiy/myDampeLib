#include "track_selection.hpp"

#include <algorithm>


myDampeLib::DmpTrackSelector::DmpTrackSelector() :
    DmpTrackSelector("bad_chan.txt")
{;}

myDampeLib::DmpTrackSelector::DmpTrackSelector(const char * file)
{
    mBadChannelsFile = file;
    mSelectTypes.resize(0);

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
            case psd_match :
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
        break;
        default :
            return true;
    }
}

void myDampeLib::DmpTrackSelector::readBadChannelsFile() {
    ifstream infile(mBadChannelsFile.c_str());  
    if (!infile)
    {
        std::cout << "Can't open Bad Channels file: " << mBadChannelsFile;
        std::cout << " ==> throwing exception!" << std::endl;
        throw;
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
}

int myDampeLib::DmpTrackSelector::hasBadChannel(DmpStkTrack * track, DmpEvent * pev) const{
    TClonesArray * stkclusters = pev->GetStkSiClusterCollection();

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

            for(int i = minc; i <= maxc; i++){
                if(mBadChannelList[ladder][i]) return ipoint;
            }
        }
    }
    return -1;
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
            if(!mc && !gPsdECor->GetPathLengthPosition(ilayer, ibar, direction, impact, len))
                continue;
            if(mc && !gPsdECor->GetPathLPMC(ilayer, ibar, direction, impact, len))
                continue;
            return true;
        }
    }
    return false;
}
