#pragma once

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

#include "TClonesArray.h"

#include "DmpEvtHeader.h"
#include "DmpRootEvent.h"
#include "DmpEvtPsdRec.h"
#include "DmpStkSiCluster.h"
#include "DmpStkTrackHelper.h"
#include "DmpSvcPsdEposCor.h"
#include "DmpStkTrack.h"
#include "DmpEvtPsdRec.h"

#define NLADDERS     192
#define NCHANNELS    384
#define NLAYERS      6
#define NPSDLAYERS   2
#define NPSDBARS     41
// #define VA           6
// #define ETA          2

// Define track cuts
#define MINTRACKHITS 6 // 5
#define MINCHISQ     6.0
#define CHI2KALMANCUT 5.


bool** read_bad_channels_file(const char* filename);
bool is_cluster_bad_channel(DmpStkSiCluster * cluster, bool** badchannels);
int reaction_in_psd(DmpStkTrackHelper * stk_helper, DmpEvtPsdRec *psdRec, bool mc=false);

class DmpTrackSelector {
public:
    DmpTrackSelector();
    DmpTrackSelector(const char * file);
    ~DmpTrackSelector();

    enum Select{stk_bad_channel,
                psd_match
                };

    void setSelectTypes(vector<Select> type) {mSelectTypes = type;}
    vector<Select> selectTypes() const {return mSelectTypes;}

    /**
     * Add selection criterium
     */
    void addSelect(Select type);

    /**
     * Returns true, if the track is selected.
     * event is the DmpEvent object, that contains the track
     */    
    bool selected(DmpStkTrack * track, DmpEvent * event) const;

    /**
     * Returns true, if the track passes the type selection
     * event is the DmpEvent object, that contains the track
     */    
    bool pass(DmpStkTrack * track, DmpEvent * event, Select type) const;

private:
    vector<Select> mSelectTypes;

    // Bad channel check
    vector<vector<bool> > mBadChannelList;
    string mBadChannelsFile;
    void readBadChannelsFile();
    bool hasBadChannel(DmpStkTrack * track, DmpEvent * event) const;

    // PSD match
    bool psdMatch(DmpStkTrack * track, DmpEvent * event) const;
};
