#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include "TH2D.h"
#include "TF1.h"
#include "DmpStkTrack.h"
#include "DmpStkSiCluster.h"

class EtaCorr
{
public:
    EtaCorr();
    EtaCorr(std::string t_file);
    ~EtaCorr();
    
    double calcEta(DmpStkSiCluster * t_cluster);
    void fit_pr(TH2D * t_h);
    void fit_pr_He(TH2D * t_h);
    bool writeCorrFile(std::string t_filename);
    bool readCorrFile(std::string t_filename);

    // Double_t corrEnergy(DmpStkTrack* t_track, DmpStkSiCluster* t_cluster);
    Double_t corrEnergy(double e, double eta, double theta);
    // Double_t corrEnergy_(DmpStkTrack* t_track, DmpStkSiCluster* t_cluster);

    void setTargetParameters(Double_t t_mp_p,  Double_t t_width_p,
							 Double_t t_mp_He, Double_t t_width_He);
    void setTargetParameters(Double_t t_mp_p,  Double_t t_mp_He);
    TF1 * getFitFunction(int t_ibin);

private:
    std::fstream m_filestr;
    int m_nbinsx;
    std::vector<Double_t> m_mp_p;
    std::vector<Double_t> m_width_p;
    std::vector<Double_t> m_mp_He;
    std::vector<Double_t> m_width_He;
    std::vector<Double_t> m_bincenters;

    bool m_p_fitted;
    bool m_He_fitted;
    bool m_p_target_set;
    bool m_He_target_set;

    std::vector<Double_t> getEstimatedLangauPar(Double_t t_eta);

    std::vector<TF1 *> m_ffit;

    // Desired parameters of the proton Landau peak after correction
    Double_t m_mp_proton;
    Double_t m_width_proton;
    Double_t m_mp_Helium;
    Double_t m_width_Helium;
};
