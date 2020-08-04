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
#include "DmpStkSiCluster.h"
#include "DmpStkTrackHelper.h"
#include "DmpStkTrack.h"
#include "DmpEvtPsdRec.h"

#include "track_selection.hpp"

#define CHI2MAX 15.

namespace myDampeLib {
    class DmpEventSelector {
    public:
        DmpEventSelector(bool check_all=false);
        ~DmpEventSelector();

        enum Select{
                    bgo_skim_cuts, // Standard set of skim cuts
                    het, // High energy trigger
                    let, // Low energy trigger
                    has_STK_track,
		    has_only_one_STK_track,
		    bgo_psd_match,
                    SELECT_N_ITEMS
                   };

        void setSelectTypes(vector<Select> types);
        vector<Select> selectTypes() const {return mSelectTypes;}

	void setTrackSelector(DmpTrackSelector * trackSelector) {mTrackSelector = trackSelector;}
	DmpTrackSelector * trackSelector() const {return mTrackSelector;}

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

	/**
	 * Returns the best STK track
	 */
	DmpStkTrack * stkTrack() const {return mSTKtrack;}
        
        /**
	 * Return the number of tracks with chi2 lower than the specified number
	 */
        int Ntracks() const { return mNtracks; }

	/**
	 * Set the distance range in mm in which the
	 * BGO direction should match the PSD hit
	 */
	void setBgoPsdDist(float d) {mBgoPsdDist = d;}

	void setTrackSelectorTypes(vector<DmpTrackSelector::Select> v) { mTrackSelector -> setSelectTypes(v); }

    private:
        vector<Select> mSelectTypes;
        TH1I * mHselect;

        bool hasSTKtrack(DmpEvent * event);
        int mNtracks;

        /*
         * Removes the side-in events by BGO 
         */
        bool notSideIn(DmpEvent * event) const;
        static constexpr float mMaxLayerRatio = 0.25;

        /*
	 * Remove events with undefined slope and interception point
	 */
        bool bgoSlopeIntercept(DmpEvent * event) const;

        /*
	 * Remove events with shower not contained in BGO volume
	 */
        bool bgoContainment(DmpEvent * event) const;

        /*
	 * in first 3 layers the max bar is not on the edge
	 */
        bool bgoMaxBar(DmpEvent * event) const;

        /*
         * Returns true if there is only one track that satisfies
	 * the mTrackSelector criteria
	 */
	bool hasOnlyOneSTKtrack(DmpEvent * event);

	float mBgoPsdDist;
	bool bgoPsdMatch(DmpEvent * event) const;

        DmpTrackSelector * mTrackSelector;
	DmpStkTrack * mSTKtrack;
    };
}

#endif
