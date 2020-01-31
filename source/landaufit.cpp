#include "landaufit.hpp"
#include <iostream>

using namespace std;

Double_t langau(Double_t *x, Double_t *par) {
    //Fit parameters:
    //par[0]=Width (scale) parameter of Landau density
    //par[1]=Most Probable (MP, location) parameter of Landau density
    //par[2]=Total area (integral -inf to inf, normalization constant)
    //par[3]=Width (sigma) of convoluted Gaussian function
    //
    //In the Landau distribution (represented by the CERNLIB approximation),
    //the maximum is located at x=-0.22278298 with the location parameter=0.
    //This shift is corrected within this function, so that the actual
    //maximum is identical to the MP parameter.

    // Numeric constants
    Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
    Double_t mpshift  = -0.22278298;       // Landau maximum location

    // Control constants
    Double_t np = 100.0;      // number of convolution steps
    Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
    // Variables
    Double_t xx;
    Double_t mpc;
    Double_t fland;
    Double_t sum = 0.0;
    Double_t xlow,xupp;
    Double_t step;
    Double_t i;


    // MP shift correction
    mpc = par[1] - mpshift * par[0];

    // Range of convolution integral
    xlow = x[0] - sc * par[3];
    xupp = x[0] + sc * par[3];

    step = (xupp-xlow) / np;

    // Convolution integral of Landau and Gaussian by sum
    for(i=1.0; i<=np/2; i++) {
        xx = xlow + (i-.5) * step;
        fland = TMath::Landau(xx,mpc,par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0],xx,par[3]);

        xx = xupp - (i-.5) * step;
        fland = TMath::Landau(xx,mpc,par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0],xx,par[3]);
    }

    return (par[2] * step * sum * invsq2pi / par[3]);
}

Double_t langau2(Double_t *x, Double_t *par) {
    Double_t par1[4], par2[4];
    for(int i=0; i<4; i++) {
        par1[i] = par[i];
        par2[i] = par[i+4];
    }
    return langau(x, par1) + langau(x, par2);
}

// Try it later: fit with two langau and a polynomial background
// Double_t langau2_p2(Double_t *x, Double_t *par) {
//     Double_t pol[3];
//     for(int i=4; i<7; i++) {
//         pol[i-4] = par[i];
//     }
//     TF1 * p = TMath::Polynomial(2);
//     p->SetParameters(pol);
//     return langau2(x, par) + p->Eval(x);
// }

TF1 * langauFit1(TH1D *h, Double_t * p0) {
    TF1 *ffit = new TF1("f", langau, 0., 500., 4);
    ffit->SetParNames("Width","MP","Area","GSigma");
    ffit->SetParameters(p0);
    h->Fit(ffit, "Q0"); // silent fit, no plotting
    return ffit;
}

TF1 * langauFit2(TH1D *h, Double_t * p0) {
    TF1 *ffit = new TF1("f", langau2, 0., 500., 8);
    ffit->SetParNames("Width_p" ,"MP_p" ,"Area_p" ,"GSigma_p",
                       "Width_He","MP_He","Area_He","GSigma_He");
    ffit->SetParameters(p0);
    h->Fit(ffit, "Q0"); // silent fit, no plotting
    return ffit;
}