#ifndef _EVENT_SELECTOR
#define _EVENT_SELECTOR

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

#include "TClonesArray.h"
#include "TH1I.h"

#include "DmpEvtHeader.h"
#include "DmpRootEvent.h"
#include "DmpEvtPsdRec.h"
#include "DmpStkSiCluster.h"
#include "DmpStkTrackHelper.h"
#include "DmpSvcPsdEposCor.h"
#include "DmpStkTrack.h"
#include "DmpEvtPsdRec.h"

#include "track_selection.hpp"


namespace myDampeLib {
    class DmpEventSelector {
    public:
        DmpEventSelector(const char * bad_chan_file, bool check_all=false);
        ~DmpEventSelector();

        enum Select{het, // High energy trigger
                    let, // Low energy trigger
                    has_STK_track,
                    has_PSD_track,
                    no_bad_clu_STK,
                    not_side_in,
                    SELECT_N_ITEMS
                   };

        void setSelectTypes(vector<Select> types);
        vector<Select> selectTypes() const {return mSelectTypes;}

        void checkAll() {
            vector<Select> types(SELECT_N_ITEMS);
            for (int i=0; i<SELECT_N_ITEMS; i++)
                types[i] = (Select)i;
            setSelectTypes(types);
        }

        /**
         * Add selection criterium
         */
        void addSelect(Select type);

        /**
         * Remove selection criterium
         */
        void removeSelect(Select type);

        /**
         * Returns true, if the track is selected.
         * event is the DmpEvent object, that contains the track
         */    
        bool selected(DmpEvent * event);

        /**
         * Returns true, if the track passes the type selection
         * event is the DmpEvent object, that contains the track
         */    
        bool pass(DmpEvent * event, Select type);

        TH1I * hSelect() const {return mHselect;}

    private:
        vector<Select> mSelectTypes;
        TH1I * mHselect;

        bool hasSTKtrack(DmpEvent * event) const;

        /*
         * Best STK track passes through both PSD layers
         */
        bool hasPSDtrack(DmpEvent * event) const;

        /*
         * Removes the side-in events by BGO 
         */
        bool notSideIn(DmpEvent * event) const;
        static constexpr float mMaxLayerRatio = 0.25;

        /*
         * Best STK track has no bad channels in first two XY layers
         * (or equivalently, first 4 X or Y layers)
         */
        bool noBadCluSTK(DmpEvent * pev) const;

        DmpTrackSelector * mTrackSelector;
    };
}

#endif
