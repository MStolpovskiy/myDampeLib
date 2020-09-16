#ifndef _PSDCHARGE
#define _PSDCHARGE

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>
#include <sstream>
#include <math.h>
#include <cmath>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TChain.h"
#include "DmpEvtHeader.h"
#include "DmpEvtPsdRec.h"
#include "DmpRootEvent.h"
#include "DmpChain.h"
#include "DmpEvtGlobTrack.h"
#include "DmpSvcPsdEposCor.h"
#include "TClonesArray.h"
#include "TSystem.h"
#include "TROOT.h"

#include "DmpStkTrack.h"
#include "DmpStkTrackHelper.h"
#include "DmpStkSiCluster.h"
#include "DmpEvtPsdRec.h"
#include "DmpEvtPsdHits.h"

namespace myDampeLib {
        double psdEnergy(DmpStkTrack * track, DmpEvtPsdHits * psdHits,
	                 int * igbar, bool mc=false, int layer=0);

	double psdCharge(double e, double proton_peak=2.07);
	double psdCharge(DmpStkTrack * track, DmpEvtPsdHits *psdHits,
	                 int * igbar, bool mc=false);
}

#endif
