#pragma once

#include "TROOT.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"

Double_t langau(Double_t *x, Double_t *par);
Double_t langau2(Double_t *x, Double_t *par);
TF1 * langauFit1(TH1D *h, Double_t *p0);
TF1 * langauFit2(TH1D *h, Double_t *p0);
