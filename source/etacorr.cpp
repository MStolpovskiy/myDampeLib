#include "etacorr.hpp"
#include "landaufit.hpp"
#include "DmpStkClusterCalibration.h"

using namespace std;

EtaCorr::EtaCorr()
    : m_nbinsx(0),
    m_p_fitted(false),
    m_He_fitted(false),
    m_p_target_set(false),
    m_He_target_set(false)
{
    ;
}

EtaCorr::EtaCorr(string t_file)
    : EtaCorr::EtaCorr()
{
    this->readCorrFile(t_file);
}

EtaCorr::~EtaCorr(){
    m_filestr.close();
}

double EtaCorr::calcEta(DmpStkSiCluster* cluster) {
    double eta;
    double s1 = 0;
    double s2 = 0;
    int argmax = -1;
    for(int istrip = 0; istrip < cluster->getNstrip(); istrip++){
        double s = cluster->GetSignal(istrip);
        if(s > s1) {
            s1 = s;
            argmax = istrip;
        }
    }
    double sm1 = 0.;
    if (argmax != 0) sm1 = cluster->GetSignal(argmax - 1);
    double sp1 = 0.;
    if (argmax != cluster->getNstrip() - 1) sp1 = cluster->GetSignal(argmax + 1);
    s2 = max(sm1, sp1);
    eta = s1 / (s1 + s2);
    return eta;
}

void EtaCorr::setTargetParameters(Double_t t_mp_p,  Double_t t_width_p,
                                  Double_t t_mp_He, Double_t t_width_He){
    this->setTargetParameters(t_mp_p, t_mp_He);
    m_width_proton = t_width_p;
    m_width_Helium = t_width_He;
    m_He_target_set = true;
}

void EtaCorr::setTargetParameters(Double_t t_mp_p, Double_t t_mp_He){
    m_mp_proton = t_mp_p;
    m_mp_Helium = t_mp_He;
    m_p_target_set = true;
}



void EtaCorr::fit_pr(TH2D * t_h){
    m_nbinsx = t_h->GetNbinsX();
    m_mp_p.resize(m_nbinsx);
    m_width_p.resize(m_nbinsx);
    m_bincenters.resize(m_nbinsx);
    m_ffit.resize(m_nbinsx);
    for (int ibin=1; ibin<=m_nbinsx; ibin++) {
        TH1D * slice = t_h->ProjectionY("", ibin, ibin);
        Double_t p0[4];
        p0[0] = slice->GetRMS();
        p0[1] = slice->GetMean();
        p0[2] = slice->GetEntries();
        p0[3] = 1.;
        m_ffit[ibin-1]  = langauFit1(slice, p0);
        m_mp_p[ibin-1]    = m_ffit[ibin-1]->GetParameter("MP");
        m_width_p[ibin-1] = m_ffit[ibin-1]->GetParameter("Width");
        m_bincenters[ibin-1] = t_h->GetXaxis()->GetBinCenter(ibin);
    }
    m_p_fitted = true;
}

void EtaCorr::fit_pr_He(TH2D * t_h) {
    // this->fit_pr(t_h);

    m_nbinsx = t_h->GetNbinsX();
    m_mp_p.resize(m_nbinsx);
    m_width_p.resize(m_nbinsx);
    m_mp_He.resize(m_nbinsx);
    m_width_He.resize(m_nbinsx);
    m_bincenters.resize(m_nbinsx);
    m_ffit.resize(m_nbinsx);

    int npar = 8;

    for (int ibin=1; ibin<=m_nbinsx; ibin++) {
        TH1D * slice = t_h->ProjectionY("", ibin, ibin);

        // Now fit all together
        Double_t p0[npar];
        bool extreme_bin = false;
        if (ibin == 1 || ibin == m_nbinsx) {
            // In the extreme bins there is no proton peak
            extreme_bin = true;
            Double_t p0_ext[4];
            p0_ext[0] = slice->GetRMS();
            p0_ext[1] = slice->GetMean();
            p0_ext[2] = slice->GetEntries();
            p0_ext[3] = 1.;
            m_ffit[ibin-1]  = langauFit1(slice, p0_ext);
            m_mp_He[ibin-1]  = m_ffit[ibin-1]->GetParameter("MP");
            m_width_He[ibin-1] = m_ffit[ibin-1]->GetParameter("Width");
            m_mp_p[ibin-1] = m_mp_He[ibin-1] / 4.;
            m_width_p[ibin-1] = m_width_He[ibin-1] / 4.;
        }
        else if (ibin == 1 || ibin == m_nbinsx) {
            // In the second-to-extreme bin the proton peak is barely visible
            p0[0] = 50.;
            p0[1] = 130.;
            p0[2] = slice->GetEntries() / 10.;
            p0[3] = 1.;
            p0[4] = 40.;
            p0[5] = 240.;
            p0[6] = slice->GetEntries();
            p0[7] = 1.;
            m_ffit[ibin-1] = langauFit2(slice, p0);
        }
        else {
            p0[0] = 10.;
            p0[1] = 50.;
            p0[2] = slice->GetEntries() * 2.;
            p0[3] = 5.;
            p0[4] = 20.;
            p0[5] = 200.;
            p0[6] = p0[2];
            p0[7] = 20.;
            m_ffit[ibin-1] = langauFit2(slice, p0);
        }

        if(!extreme_bin) {
            m_mp_p[ibin-1]  = m_ffit[ibin-1]->GetParameter("MP_p");
            m_width_p[ibin-1] = m_ffit[ibin-1]->GetParameter("Width_p");
            m_mp_He[ibin-1] = m_ffit[ibin-1]->GetParameter("MP_He");
            m_width_He[ibin-1] = m_ffit[ibin-1]->GetParameter("Width_He");
        }
        m_bincenters[ibin-1] = t_h->GetXaxis()->GetBinCenter(ibin);
    }
    m_p_fitted = true;
    m_He_fitted = true;
}

TF1 * EtaCorr::getFitFunction(int t_ibin){
    ///\brief Returns fit function of i-th bin. Note that bins are numbered from 0 to Nbins-1.
    if (t_ibin >= m_nbinsx || t_ibin < 0) throw "Wrong bin number!";
    return m_ffit[t_ibin];
}

vector<Double_t> EtaCorr::getEstimatedLangauPar(Double_t t_eta){
    // Return the proton and helium peak maximum and width
    // expected for the given value of eta

    // Simplest approach for the moment:
    // no interpolation, just take values in bins
    if (!m_p_fitted) {
        throw "Fit at least the proton peak, or load the parameter file";
    }
    vector<Double_t>::iterator const it = lower_bound(m_bincenters.begin(), m_bincenters.end(), t_eta);
    int ind = it - m_bincenters.begin();
    int npar = m_He_fitted ? 4 : 2;
    vector<Double_t> par(npar);
    par[0] = m_mp_p[ind];
    par[1] = m_width_p[ind];
    if (m_He_fitted){
        par[2] = m_mp_He[ind];
        par[3] = m_width_He[ind];
    }
    return par;
}

//Double_t EtaCorr::corrEnergy(DmpStkTrack* t_track, DmpStkSiCluster* t_cluster) {
//    if (!m_p_target_set) {
//        throw "Set the target value at least for the proton";
//    }
//    double costheta = t_track->getDirection().CosTheta();
//    double e = t_cluster -> getEnergy() * costheta;
//    double eta = CalcEta(t_cluster);
//    vector<Double_t> par = this->getEstimatedLangauPar(eta);
//    Double_t mp_p = par[0];
//    Double_t width_p = par[1];
//    Double_t mp_He;
//    Double_t width_He;
//    if (m_He_fitted) {
//        mp_He = par[2];
//        width_He = par[3];
//    }
//
//    // Correction for both proton and helium peak positions
//    Double_t proportion;
//    if (m_He_target_set) proportion = (m_mp_proton / mp_p + m_mp_Helium / mp_He) / 2.;
//    else proportion = m_mp_proton / mp_p;
//    Double_t ecorr = e * proportion;
//    return ecorr;
//}

Double_t EtaCorr::corrEnergy(double e, double eta, double theta) {
    // returns the cluster energy corrected
    // eta -- s_max / (s_max + s_neighbour)
    // theta -- abs of either theta_x or theta_y in degrees
    double pal[3] = { 5.58971156e+01, -1.19714996e+00,  2.51431862e-02 };
    double pbl[3] = {-4.27545738e+01,  2.55264953e+00, -3.98720934e-02};
    double par[2] = {-26.31069995,   1.70947931};
    double pbr[2] = {92.57422316, -2.16547319};
    double cnstl = pal[0] + pal[1] * theta + pal[2] * theta * theta;
    double slopl = pbl[0] + pbl[1] * theta + pbl[2] * theta * theta;
    double cnstr = par[0] + par[1] * theta;
    double slopr = pbr[0] + pbr[1] * theta;

    double fl = cnstl + slopl * eta;
    double fr = cnstr + slopr * eta;

    double f = max(fl, fr);

    double corr = this->m_mp_proton / f;
    return corr * e;
}

// Double_t EtaCorr::corrEnergy_(DmpStkTrack* t_track, DmpStkSiCluster* t_cluster) {
//     double costheta = t_track->getDirection().CosTheta();
//     double e = t_cluster -> getEnergy() * costheta;
//     double eta = CalcEta(t_cluster);
//     vector<Double_t> par = this->getEstimatedLangauPar(eta);
//     Double_t mp_p = par[0];
//     Double_t width_p = par[1];
//     Double_t mp_He = par[2];
//     Double_t width_He = par[3];
//     // Correction for proton peak position and width
//     Double_t ecorr = (e - mp_p) / width_p * m_width_proton + m_mp_proton;
//     return ecorr;
// }

bool EtaCorr::writeCorrFile(string t_filename){
    m_filestr.open(t_filename.c_str(), fstream::out);
    if(!m_filestr.is_open()) return false;
    double eta, pr, pr_sig, he, he_sig;
    for(int i=0; i<m_bincenters.size(); i++) {
        pr = m_mp_p[i];
        pr_sig = m_width_p[i];
        he = m_mp_He[i];
        he_sig = m_width_He[i];
        eta = m_bincenters[i];

        m_filestr << eta << " " << pr << " " << pr_sig << " " << he << " " << he_sig << endl;
    }
    return true;
}

bool EtaCorr::readCorrFile(string t_filename){
    m_filestr.open(t_filename.c_str(), fstream::in);
    if(!m_filestr.is_open()) return false;
    double eta, pr, pr_sig, he, he_sig;

    while(m_filestr >> eta >> pr >> pr_sig >> he >> he_sig) {
        m_mp_p.push_back(pr);
        m_width_p.push_back(pr_sig);
        m_mp_He.push_back(he);
        m_width_He.push_back(he_sig);
        m_bincenters.push_back(eta);
    }
    m_nbinsx = m_bincenters.size();
    m_p_fitted = true;
    m_He_fitted = true;

    return true;
}

