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

// #include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"

#include "DmpStkTrack.h"
#include "DmpStkTrackHelper.h"
#include "DmpStkSiCluster.h"
#include "DmpEvtPsdRec.h"
#include "TStopwatch.h"
#include "DmpStkClusterCalibration.h"

#include "track_selection.hpp"
#include "etacorr.hpp"

using namespace std;
using namespace TMVA;

namespace myDampeLib {
    class DmpAnalysis {
    public:
        DmpAnalysis();
        DmpAnalysis(const char * filename);
        DmpAnalysis(string filename);
        ~DmpAnalysis();

        void setOutputFilename(const char * filename) {mOutputFilename = filename;}
        void setOutputFilename(string filename) {mOutputFilename = filename;}
        string outputFilename() const {return mOutputFilename;}

        /*
        Set TChain of the input files from a list.txt file
        */
        void setTChain(const char * filename, bool verbose=true);

        /*
        Add an element to the TChain of the input files
        */
        void add2TChain(string filename, bool verbose=true)

        void openOutputFile(string option="RECREATE");
        void closeOutputFile();

    private:
        string mOutputFilename;
        TFile * mOutputFile;

        TChain * mChain;
        int mNFilesChained;
        int mNEvents;
    };
}