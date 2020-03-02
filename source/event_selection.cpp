#include "event_selection.hpp"
#include "psd_charge.hpp"
#include "definitions.hpp"

myDampeLib::DmpEventSelector::DmpEventSelector(const char * bad_chan_file,
                                               bool check_all/*=false*/)
{
    mSelectTypes.resize(0);
    mHselect = new TH1I("Hselect", "Selection criteria", SELECT_N_ITEMS + 1, 0, SELECT_N_ITEMS + 1);
    mHselect->GetXaxis()->SetBinLabel(1, "All");
    for (int i=1; i <= SELECT_N_ITEMS; i++) {
        Select type = (Select)(i-1);
	string type_str;
	switch (type) {
	    case het : type_str = "HET"; break;
	    case let : type_str = "LET"; break;
	    case has_STK_track : type_str = "has STK track"; break;
	    case has_PSD_track : type_str = "has PSD track"; break;
        case no_bad_clu_STK : type_str = "no bad STK channels"; break;
	    case not_side_in   : type_str = "not side-in"; break;
	    default: type_str = "";
	}
        mHselect->GetXaxis()->SetBinLabel(i+1, type_str.c_str());
    }
    if (check_all) checkAll();

    mTrackSelector = new DmpTrackSelector(bad_chan_file);
}

myDampeLib::DmpEventSelector::~DmpEventSelector()
{
    delete mTrackSelector;
}

void myDampeLib::DmpEventSelector::setSelectTypes(vector<Select> types)
{
    for (vector<Select>::iterator it = types.begin(); it != types.end(); it++) {
        addSelect(*it);
    }
}

void myDampeLib::DmpEventSelector::addSelect(Select type)
{
    vector<Select>::iterator it;
    it = find(mSelectTypes.begin(), mSelectTypes.end(), type);
    if (it == mSelectTypes.end())
        mSelectTypes.push_back(type);
}

void myDampeLib::DmpEventSelector::removeSelect(Select type)
{
    vector<Select>::iterator it;
    it = find(mSelectTypes.begin(), mSelectTypes.end(), type);
    if (it != mSelectTypes.end())
        mSelectTypes.erase(it);
}

bool myDampeLib::DmpEventSelector::selected(DmpEvent * event)
{
    bool ret;
    mHselect->Fill(0);
    for (vector<Select>::const_iterator it = mSelectTypes.begin();
         it != mSelectTypes.end(); 
         ++it) {
        ret = pass(event, *it);
	if (ret) mHselect->Fill(*it + 1);
	else break;
    }
    return ret;
}

bool myDampeLib::DmpEventSelector::pass(DmpEvent * event, Select type)
{
    bool ret = true;
    switch (type) {
        case het :
            ret = event->pEvtHeader()->GeneratedTrigger(3); break;
        case let :
            ret = event->pEvtHeader()->GeneratedTrigger(2); break;
        case has_STK_track :
            ret = hasSTKtrack(event); break;
        case has_PSD_track :
            ret = hasPSDtrack(event); break;
        case no_bad_clu_STK :
            ret = noBadCluSTK(event); break;
        case not_side_in :
            ret = notSideIn(event); break;
        default :
            ret = false;
    }
    return ret;
}

bool myDampeLib::DmpEventSelector::hasSTKtrack(DmpEvent * pev)
{
    DmpStkTrackHelper * stk_helper = new DmpStkTrackHelper(pev->GetStkKalmanTrackCollection (),
                                                           true, 
                                                           pev->pEvtBgoRec(),
                                                           pev->pEvtBgoHits());
    stk_helper -> SortTracks(3, true); // most strigent selection
    int ntracks = stk_helper->GetSize();
    if (ntracks > 0) mSTKtrack = stk_helper->GetTrack(0);
    delete stk_helper;
    return ntracks > 0;
}

bool myDampeLib::DmpEventSelector::hasPSDtrack(DmpEvent * pev) const
{
    return mTrackSelector->pass(mSTKtrack, pev, DmpTrackSelector::psd_match);
}

bool myDampeLib::DmpEventSelector::noBadCluSTK(DmpEvent * pev) const
{
    int bad_layer = mTrackSelector->hasBadChannel(mSTKtrack, pev);
    return (bad_layer == -1) || (bad_layer >= 4);
}

bool myDampeLib::DmpEventSelector::notSideIn(DmpEvent * pev) const
{
    DmpEvtBgoRec* bgoRec = pev->pEvtBgoRec();
    double bgo_layer_sum[NLAYERSTOT];
    for (int ilayer=0; ilayer<NLAYERSXY; ilayer++){
        bgo_layer_sum[ilayer * 2] = 0.;
        bgo_layer_sum[ilayer * 2 + 1] = 0.;
        for (int ibar=0; ibar<NBARSL; ibar++) {
            bgo_layer_sum[ilayer * 2] += bgoRec->GetEdep(ilayer * 2, ibar);
            bgo_layer_sum[ilayer * 2 + 1] += bgoRec->GetEdep(ilayer * 2 + 1, ibar);
        }
    }

    // Energy_L2 + Energy_L3 > Energy_L0 + Energy_L1 To exclude Down-to-Up events
    if (bgo_layer_sum[0] + bgo_layer_sum[1] >
        bgo_layer_sum[2] + bgo_layer_sum[3]) return false;

    // MaxLayerEnergyRatio < 0.35 To exclude Side-In events
    double max_layer_ratio = 0.;
    double bgo_energy = bgoRec->GetTotalEnergy();
    for (int ilayer=0; ilayer<NLAYERSTOT; ilayer++) {
        double r = bgo_layer_sum[ilayer] / bgo_energy;
        if (max_layer_ratio < r) max_layer_ratio = r;
    }
    if (max_layer_ratio >= mMaxLayerRatio) return false;

    return true;
}
