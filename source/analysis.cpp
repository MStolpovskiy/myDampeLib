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
    mOutputFilename(filename)
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