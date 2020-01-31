#include "track_selection.hpp"


bool** read_bad_channels_file(const char* filename){
	bool** badchannels = new bool*[NLADDERS];
	for(int i=0; i<NLADDERS; i++){
		badchannels[i] = new bool[NCHANNELS];
		for(int j=0; j<NCHANNELS; j++){
			badchannels[i][j] = false;
		}
	}

	ifstream infile(filename);	
	if (!infile)
	{
		std::cout<<"Can't open Bad Channels file: "<<filename<<" ==> throwing exception!"<<std::endl;
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
			badchannels[ildr][bdchnl] = true;
		}
	}
	std::cout<<"Done reading bad channels"<<std::endl;
	return badchannels;
} // end read_bad_channels_file


bool is_cluster_bad_channel(DmpStkSiCluster * cluster, bool** badchannels) {
    int ladder = cluster->getLadderHardware(); //ladder num
    int minc = cluster->getIndex1();//address of the first strip in the cluster (0-363)
    int maxc = cluster->getIndex1() + cluster->getNstrip() -1;//last strip of the cluster

	for(int i = minc; i <= maxc; i++){
		if(badchannels[ladder][i]) return true;
	}
	return false;
} // end is_cluster_bad_channel


int reaction_in_psd(DmpStkTrackHelper * stk_helper, DmpEvtPsdRec *psdRec, bool mc) {
	// returns 0 if it seems like the first inelastic interaction
	// took place in the first PSD layer or above it
	// returns 1 if the first II took place in the second layer of the PSD
	// If the first interaction was in STK or below, return 2
	// If there is no track, return 3
    int ntracks = stk_helper->GetSize();

    if (ntracks == 0) return 3;
    if (ntracks == 1) return 2;

    int count_fired_bars[2] = {0, 0};

    for(int tck = 0; tck < stk_helper->GetSize(); tck++) {
        DmpStkTrack* stktrack = stk_helper->GetTrack(tck);

	    TVector3 impact = stktrack->getImpactPoint();
	    TVector3 direction = stktrack->getDirection();

	    for (int ilayer = 0; ilayer < 2; ilayer++) {
	        // Loop over psd bars
	        for (int ibar = 0; ibar < 41; ibar++) {
	        	double e = psdRec->GetEdep(ilayer, ibar);
	        	if (e == 0) continue;
	            double * len = new double[2];
	            if(!mc && !gPsdECor->GetPathLengthPosition(ilayer, ibar, direction, impact, len))
	                continue;
	            if(mc && !gPsdECor->GetPathLPMC(ilayer, ibar, direction, impact, len))
	                continue;
                count_fired_bars[ilayer] ++;
			}
		}
		if (count_fired_bars[0] > 3) return 0;
		if (count_fired_bars[1] > 3) return 1;
    }
    return 2;
}