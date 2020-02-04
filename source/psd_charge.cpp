#include "psd_charge.hpp"

double myDampeLib::psdEnergy(DmpStkTrack * track, DmpEvtPsdRec *psdRec,
                 int * ibar1, int * ibar2, bool mc/*=false*/, int layer/*=0*/) {
    // layer -- get energy from a specific psd layer

    TVector3 globImpact = track->getImpactPoint();
    TVector3 globDirection = track->getDirection();

    double e = 0.;
    int nhits = 0;

    // Loop over psd layers
    for (int ilayer = 0; ilayer < 2; ilayer++) {
        if (layer != 0 && ilayer != layer-1) continue;
        // Loop over psd bars
        for (int ibar = 0; ibar < 41; ibar++) {
            // psd energy, no correction
            double etemp = psdRec->GetEdep(ilayer, ibar);

            // Get hit point and path length
            double * len = new double[2];
            if(!mc && !gPsdECor->GetPathLengthPosition(ilayer, ibar, globDirection, globImpact, len))
                continue;
            if(mc && !gPsdECor->GetPathLPMC(ilayer, ibar, globDirection, globImpact, len))
                continue;
            // len[0] -- hit point
            // len[1] -- path length

            switch (ilayer) {
                case 0 : *ibar1 = ibar; break;
                case 1 : *ibar2 = ibar; break;
            }

            // correction for the path length
            etemp *= 10. / len[1]; // 10. is the psd bar thickness in mm

            // correction for the attenuation
            double corr = gPsdECor -> GetPsdECorSp3(ilayer, ibar, len[0]);
            // double corr = gPsdECor -> GetPsdECor(ilayer, ibar, len[0]);
            etemp *= corr;

            e += etemp;
            nhits += 1;
        } // end bar loop
    } // end layer loop
    double ret = e / nhits;
    if(!std::isnormal(ret)) ret = 0.;
    return ret;
}

double myDampeLib::psdCharge(double e, double proton_peak/*=2.07*/) {
    // proton_peak -- position of the proton peak in the dE/dx distribution
    return TMath::Sqrt(e / proton_peak);
}

double myDampeLib::psdCharge(DmpStkTrack * track, DmpEvtPsdRec *psdRec,
                 int * ibar1, int * ibar2, bool mc/*=false*/) {
    double e = psdEnergy(track, psdRec, ibar1, ibar2, mc);
    return psdCharge(e);
}
