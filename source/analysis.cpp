#include "analysis.hpp"

myDampeLib::DmpAnalysis::DmpAnalysis():
    DmpAnalysis("default.root")
{
}

myDampeLib::DmpAnalysis::DmpAnalysis(const char * filename) :
    DmpAnalysis(string(filename))
{
}

myDampeLib::DmpAnalysis::DmpAnalysis(string filename) :
    mOutputFilename(filename),
    mNFilesChained(0)
{
    openOutputFile(mOutputFilename);
}

myDampeLib::DmpAnalysis::~DmpAnalysis()
{
    closeOutputFile();
}

void myDampeLib::DmpAnalysis::openOutputFile(string option/*="RECREATE"*/)
{
    mOutputFile = new TFile(mOutputFilename.c_str(), option.c_str());
}

void myDampeLib::DmpAnalysis::closeOutputFile()
{
    mOutputFile->Close();
}

void myDampeLib::sDmpAnalysis::etTChain(const char * filename, bool verbose/*=true*/)
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
}