#pragma once

#include <fstream>
#include <iostream>
#include <sstream>
#include "TClonesArray.h"

#include "DmpEvtHeader.h"
#include "DmpRootEvent.h"
#include "DmpEvtPsdRec.h"
#include "DmpStkSiCluster.h"
#include "DmpStkTrackHelper.h"
#include "DmpSvcPsdEposCor.h"

#define NLADDERS     192
#define NCHANNELS    384
#define NLAYERS      6
#define VA           6
#define ETA          2

// Define track cuts
#define MINTRACKHITS 6 // 5
#define MINCHISQ     6.0
#define CHI2KALMANCUT 5.


bool** read_bad_channels_file(const char* filename);
bool is_cluster_bad_channel(DmpStkSiCluster * cluster, bool** badchannels);
int reaction_in_psd(DmpStkTrackHelper * stk_helper, DmpEvtPsdRec *psdRec, bool mc=false);
