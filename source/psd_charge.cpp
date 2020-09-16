#include "psd_charge.hpp"

double myDampeLib::psdEnergy(DmpStkTrack * track, DmpEvtPsdHits * psdHits,
                             int * igbar, bool mc/*=false*/, int layer/*=0*/) {
    // layer -- get energy from a specific psd layer (1 or 2)

    TVector3 globImpact = track->getImpactPoint();
    TVector3 globDirection = track->getDirection();

    double e = 0.;
    int nhits = 0;

    // Loop over psd hits
    for (int ihit = 0; ihit < psdHits->GetHittedBarNumber(); ihit++) {
	int ilayer = psdHits->GetLayerID(ihit);
        if (layer != 0 && ilayer != layer-1) continue;

	short gid = (psdHits->fGlobalBarID)[ihit];

	// psd energy, no correction
	double etemp = (psdHits->fEnergy)[ihit];

	// Get hit point and path length
        // len[0] -- hit point
        // len[1] -- path length
	double * len = new double[2];
	TVector3 hpos(-1000,-1000,-1000);

	bool isHit = false;
	if (!mc) {
	    isHit = gPsdECor->GetPathLengthPoint(gid, globDirection, globImpact, len, hpos);
	} else {
	    isHit = gPsdECor->GetPathLPMCPoint  (gid, globDirection, globImpact, len, hpos);
	}
	if(!isHit) continue; // no match between PSD and STK
            
	*igbar = gid;

        // correction for the path length
        etemp *= 10. / len[1]; // 10. is the psd bar thickness in mm

        // correction for the attenuation
	float pos = (ilayer==0) ? hpos.x() : hpos.y();
	double corr = gPsdECor->GetPsdECorSp3(gid, pos);
	//double corr = gPsdECor->GetPsdMipAttE(gid, pos);
	etemp *= corr;

	e += etemp;
	nhits += 1;
    } // end hit loop
    double ret = e / nhits;
    if(!std::isnormal(ret)) ret = 0.;
    return ret;
}

double myDampeLib::psdCharge(double e, double proton_peak/*=0.07*/) {
    // proton_peak -- position of the proton peak in the dE/dx distribution
    return TMath::Sqrt(e / proton_peak);
}

double myDampeLib::psdCharge(DmpStkTrack * track, DmpEvtPsdHits *psdHits,
                             int * igbar, bool mc/*=false*/) {
    double e = psdEnergy(track, psdHits, igbar, mc);
    return psdCharge(e);
}
