#include "analysis.hpp"

myDampeLib::DmpAnalysis::DmpAnalysis():
    mOutputFilename("")
{
}

myDampeLib::DmpAnalysis::DmpAnalysis(const char * filename) :
    DmpAnalysis()
{
    DmpAnalysis(string(filename));
}

myDampeLib::DmpAnalysis::DmpAnalysis(string filename) :
    DmpAnalysis(),
    mOutputFilename(filename)
{
}

myDampeLib::DmpAnalysis::~DmpAnalysis()
{
    if (mOutputFile) {
        closeOutputFile();
    }
}

void myDampeLib::DmpAnalysis::openOutputFile(string option="RECREATE")
{
    if (mOutputFilename != "") {
        mOutputFile = new TFile(mOutputFilename.c_str(), option.c_str());
    }
    else {
        std::cout << "mOuputFilename is not defined" << std::endl;
    }
}

void myDampeLib::DmpAnalysis::closeOutputFile()
{
    mOutputFile.Close();
}