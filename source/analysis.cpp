#include "analysis.hpp"

#define GET_VARIABLE_NAME(Variable) (void(Variable),#Variable)

myDampeLib::DmpAnalysis::DmpAnalysis():
    DmpAnalysis("default.root")
{
}

myDampeLib::DmpAnalysis::DmpAnalysis(const char * filename, const char * option) :
    DmpAnalysis(string(filename), string(option))
{
}

myDampeLib::DmpAnalysis::DmpAnalysis(string filename, string option/*=RECREATE*/) :
    mOutputFilename(filename),
    mNFilesChained(0)
{
    openOutputFile(option);
}

myDampeLib::DmpAnalysis::~DmpAnalysis()
{
    closeOutputFile();
}

void myDampeLib::DmpAnalysis::openOutputFile(string option/*="RECREATE"*/)
{
    mOutputFile = new TFile(mOutputFilename.c_str(), option.c_str());
    if (mOutputFile -> IsZombie()) {
        std::cout << "Error opening output file" << std::endl;
        exit(-1);
    }
    std::cout << "Output file open " << mOutputFilename << std::endl;
}

void myDampeLib::DmpAnalysis::closeOutputFile()
{
    mOutputFile->Close();
    // delete mOutputFile;
}

void myDampeLib::DmpAnalysis::setTChain(const char * filename, bool verbose/*=true*/)
{
    mChain = new DmpChain("CollectionTree");
    ifstream runlist(filename);
    string line;
    mNFilesChained = 0;

    while(getline(runlist, line)) {
        add2TChain(line, verbose);
    }

    runlist.close();
}

void myDampeLib::DmpAnalysis::add2TChain(string filename, bool verbose/*=true*/)
{
    TFile * f = TFile::Open(filename.c_str(), "READ");
    if (f) {
        if (verbose) std::cout << "File found " << filename << endl;
        mChain -> Add(filename.c_str());
        mNFilesChained++;
    }

    mNEvents = mChain->GetEntries();
}

void myDampeLib::DmpAnalysis::addBranch(int var)
{
    string t = string(GET_VARIABLE_NAME(var)) + "/I";
    mTree->Branch(GET_VARIABLE_NAME(var), &var, t.c_str());
}

void myDampeLib::DmpAnalysis::addBranch(float var)
{
    string t = string(GET_VARIABLE_NAME(var)) + "/F";
    mTree->Branch(GET_VARIABLE_NAME(var), &var, t.c_str());
}

void myDampeLib::DmpAnalysis::addBranch(double var)
{
    string t = string(GET_VARIABLE_NAME(var)) + "/D";
    mTree->Branch(GET_VARIABLE_NAME(var), &var, t.c_str());
}

void myDampeLib::DmpAnalysis::addBranch(int var[])
{
    int len = sizeof(var) / sizeof(*var);
    char * lenStr = itoa(len)
    string t = string(GET_VARIABLE_NAME(var)) + "[" + lenStr + "]/I";
    mTree->Branch(GET_VARIABLE_NAME(var), &var, t.c_str());
}

void myDampeLib::DmpAnalysis::addBranch(float var[])
{
    int len = sizeof(var) / sizeof(*var);
    char * lenStr = itoa(len)
    string t = string(GET_VARIABLE_NAME(var)) + "[" + lenStr + "]/F";
    mTree->Branch(GET_VARIABLE_NAME(var), &var, t.c_str());
}

void myDampeLib::DmpAnalysis::addBranch(double var[])
{
    int len = sizeof(var) / sizeof(*var);
    char * lenStr = itoa(len)
    string t = string(GET_VARIABLE_NAME(var)) + "[" + lenStr + "]/D";
    mTree->Branch(GET_VARIABLE_NAME(var), &var, t.c_str());
}

void myDampeLib::DmpAnalysis::run(int n/*=-1*/)
{ 
    int nnn = (n < 0)? mNEvents : n;
    for(mCurrentEvent = 0; mCurrentEvent < nnn; mCurrentEvent++) {
        if ((mCurrentEvent % (nnn / 10 + 1)) == int(nnn / 10)) {
            int percentage = 10. * int((mCurrentEvent + 1) / int(nnn / 10));
            cout << "Processing percentage: " << percentage << "% \n";
        }

        analyseOneEvent();
        mTree -> Fill();
    }
}
