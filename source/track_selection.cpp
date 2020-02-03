#include "track_selection.hpp"


DmpTrackSelector::DmpTrackSelector() 
{
    mSelectTypes.resize(0);

    mBadChannelList.resize(NLADDERS);
    for (vector<vector<bool>>::iterator it = mBadChannelList.begin();
         it != mBadChannelList.end();
         ++it) {
        it->resize(NCHANNELS);
    }
}

DmpTrackSelector::DmpTrackSelector(const char * file) :
    DmpTrackSelector()
{
    mBadChannelsFile = file;
} 

void DmpTrackSelector::addSelect(Select type){
    switch (type) {
        case stk_bad_channel :
            readBadChannelsFile(mBadChannelsFile);
            break;
        case psd_match :
        default :
            break;
    }
    mSelectTypes.push_bash type;
}

bool DmpTrackSelector::selected(DmpStkTrack* track) const {
    for (vector<Select>::iterator it = mSelectTypes.begin();
         it != mSelectTypes.end(); 
         ++it) {
        if(!pass(track, *it)) return false;
    }
    return true
}

bool DmpTrackSelector::pass(DmpStkTrack * track, DmpEvent * event, Select type) const {
    switch (type) {
        case stk_bad_channel :
            return hasBadChannel(track);
            break;
        case psdMatch :
            return psdMatch(track);
        default :
            return true;
    }
}

void DmpTrackSelector::readBadChannelsFile() {
    ifstream infile(mBadChannelsFile);  
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

bool DmpTrackSelector::hasBadChannel(DmpStkTrack * track) const{
    // Loop over clusters
    for (int ipoint=0; ipoint<stktrack->GetNPoints(); ipoint++) {
        for (int ixy=0; ixy<2; ixy++) {
            DmpStkSiCluster* cluster;
            if(ixy == 0)
                cluster = stktrack -> GetClusterX(ipoint, stkclusters);
            else
                cluster = stktrack -> GetClusterY(ipoint, stkclusters);
            if(!cluster) continue;

            //ladder num
            int ladder = cluster->getLadderHardware();

            //address of the first strip in the cluster (0-363)
            int minc = cluster->getIndex1();

            //last strip of the cluster
            int maxc = cluster->getIndex1() + cluster->getNstrip() -1;

            for(int i = minc; i <= maxc; i++){
                if(mBadChannelList[ladder][i]) return true;
            }
        }
    }
    return false;
}

bool DmpTrackSelector::psdMatch(DmpStkTrack * tracl, DmpEvent * event) const {
    DmpEvtPsdRec *psdRec = pev->pEvtPsdRec();

    TVector3 impact = stktrack->getImpactPoint();
    TVector3 direction = stktrack->getDirection();

    bool mc = (pev->pEvtSimuHeader() != nullptr);

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
