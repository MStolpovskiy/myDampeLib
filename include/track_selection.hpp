#ifndef _TRACKSELECTION
#define _TRACKSELECTION

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <math.h>

#include "TClonesArray.h"

#include "DmpEvtHeader.h"
#include "DmpRootEvent.h"
#include "DmpEvtPsdRec.h"
#include "DmpStkSiCluster.h"
#include "DmpStkTrackHelper.h"
#include "DmpSvcPsdEposCor.h"
#include "DmpStkTrack.h"
#include "DmpEvtPsdRec.h"

#include "psd_charge.hpp"

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

namespace myDampeLib {

    class DmpTrackSelector {
    public:
        DmpTrackSelector();
        DmpTrackSelector(const char * file);
        ~DmpTrackSelector();

	/**
	 * Available selection criteria:
	 * 
	 * stk_bad_channel -- remove tracks that pass through the bad channel
	 *                    in first two layers of STK
	 * stk_xy_overlap -- remove tracks with too many nonoverlaps.
	 *                   use setXYnonoverlaps to specify the limit
	 * stk_missing_impact_point -- remove tracks for which the impact point is measured on
	 *                             1 or 2 (set by setNimpactPoints()) planes.
	 * stk_impact_first_layer -- impact point must be in the first layer of STK
	 * stk_bgo_dist_high -- select tracks that match the BGO shower on the first BGO layer
	 *                      use setMaxDistStk2Bgo to specify the distance limit
	 *                      use setMaxAngStk2Bgo to specify the angle limit
	 * stk_chi2_cut -- select tracks with chi2 / ndf below limit
	 *                 use setSTKtrackChi2Max to specify the limit
	 * stk_npoints -- select tracks with at least N extrapolated points
	 *                use setMinNPoints to specify the limit
	 *                i.e. if setMinNPoints(6) the track by its direction passes
	 *                through all the layers of STK
         * stk_nXYhits -- select track with at least N real XY hits
	 *                use setMinNXYhits to specify the limit
	 * stk_no_1strip_clusters -- remove tracks that contain single-strip clusters
	 * proton_mip -- select proton MIP tracks
	 *               use setSTKProtonMipE to set the mpv ADC value for the proton Landau peak
	 *                   typical value is 60.
	 *               use setSTKProtonMipRange to set the relative width of the Landau peak
	 *                   typical value is 0.1
	 * psd_match -- STK track matches a PSD cluster
	 * match_primary -- Using a linear model trained on the Proton and Helium MC
	 *                  we reconstruct the distance from the track to the primary particle
	 *                  use setDistFromPrimCut to set the distance cut
	 *                      typical value is around 0.1-0.
	 */
        enum Select{stk_bad_channel,
		    stk_xy_overlap,
		    stk_missing_impact_point,
		    stk_impact_first_layer,
		    stk_bgo_dist_high,
		    stk_chi2_cut,
		    stk_npoints,
		    stk_nXYhits,
		    stk_no_1strip_clusters,
		    proton_mip,
                    psd_match,
		    match_primary
                    };

	enum Compare{standard,
		     prim_dist};

	/**
	 * Specify track selection criteria with a vector of type Select
	 */
        void setSelectTypes(vector<Select> types);
        vector<Select> selectTypes() const {return mSelectTypes;}

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
        bool selected(DmpStkTrack * track, DmpEvent * event) const;

        /**
         * Returns true, if the track passes the type selection
         * event is the DmpEvent object, that contains the track
         */    
        bool pass(DmpStkTrack * track, DmpEvent * event, Select type) const;

        /**
         * Returns the number of the first 
         * encountered STK layer (from 0 to 11) in which
         * the cluster contains a bad channel
         * If no bad channels, return -1
         */
        int hasBadChannel(DmpStkTrack * track, DmpEvent * event) const;

	/**
	 * Returns a vector of layers where cluster contains a bad channel
	 */
        vector<int> layerBadChannel(DmpStkTrack * track, DmpEvent * event) const;

        void setBadChannelsFile(string file);
        bool readBadChannelsFile();

	void setXYnonoverlaps(int n) {mXYnonoverlaps = n;}
	int XYnonoverlaps() const {return mXYnonoverlaps;}

	void setNimpactPoints(int n) {
            if (n == 1 || n == 2) mNimpactPoints = n; 
	    else std::cerr << "Error, N impact Points must be 1 or 2" << std::endl;
        }
	int nImpactPoints() const {return mNimpactPoints;}

	void setMaxDistStk2Bgo(float d) {mMaxDistStk2Bgo = d;}
	float maxDistStk2Bgo() const {return mMaxDistStk2Bgo;}

        void setMaxAngStk2Bgo(float d) {mMaxAngStk2Bgo = d;}
        float maxAngStk2Bgo() const {return mMaxAngStk2Bgo;}

	/*
	 * Return the angle between stk track and bgo shower
	 */
	float stk2bgoAngl(DmpStkTrack * track, DmpEvent * event) const;

        /*
         * Return distance on top of BGO between STK track and BGO shower
         */
	float stk2bgoDist(DmpStkTrack * track, DmpEvent * event) const;

	void setSTKtrackChi2Max(float chi2) { mSTKtrackChi2Max = chi2; }
        float STKtrackChi2Max() const { return mSTKtrackChi2Max; }

	void setMinNPoints(int n) { mMinNPoints = n; }
	int minNPoints() const { return mMinNPoints; }

	void setMinNXYhits(int n) { mMinNXYhits = n; }
	int minNXYhits() const { return mMinNXYhits; }
       
	void setSTKProtonMipE(float e) {mSTKProtonMipE = e;}
	float stkProtonMipE() const {return mSTKProtonMipE;}

	void setSTKXProtonMipRange(float r) {mSTKXProtonMipRange = r;}
	float stkXProtonMipRange() const {return mSTKXProtonMipRange;}

	void setSTKYProtonMipRange(float r) {mSTKYProtonMipRange = r;}
	float stkYProtonMipRange() const {return mSTKYProtonMipRange;}

	void setDistFromPrimCut(float d) {mDistFromPrimCut = d;}
	float distFromPrimCut() const {return mDistFromPrimCut;}

	void getSMean(DmpStkTrack * track, DmpEvent * event, float * s_mean_xy) const;

        /*
	 * Compare the track quality
	 * first_is_better -- comparison based on chi2, npoints etc.
	 * first_is_closer -- comparison based on the reconstructed distance from the primary track
	 */
	void setComparison(Compare c) {mCompare = c;}
	Compare comparison() const {return mCompare;}
	bool first_is_better(DmpStkTrack * lhs, DmpStkTrack * rhs, DmpEvent * event) const;
	bool first_is_closer(DmpStkTrack * lhs, DmpStkTrack * rhs, DmpEvent * event) const;

	/*
	 * Reconstruct the distance to the primary particle
	 * Reconstruction is based on the linear model trained separately in Python
	 */
	float dist_from_prim(DmpStkTrack * track, DmpEvent * pev) const;

    private:
        vector<Select> mSelectTypes;
	Compare mCompare;

        // Bad channel check
        vector<vector<bool> > mBadChannelList;
        string mBadChannelsFile;
        
        // PSD match
        bool psdMatch(DmpStkTrack * track, DmpEvent * event) const;

        // STK track quality criteria:
	int mXYnonoverlaps;
	bool stkXYoverlap(DmpStkTrack * track) const;
	int mNimpactPoints;
	bool stkNoMissingImpactPoint(DmpStkTrack * track) const;
	bool stkImpactFirstLayer(DmpStkTrack * track) const;
	float mMaxDistStk2Bgo;
	float mMaxAngStk2Bgo;
	bool stk2bgoCut(DmpStkTrack * track, DmpEvent * event) const;
	float mSTKtrackChi2Max;
	bool stkChi2Cut(DmpStkTrack * track) const;
	int mMinNPoints;
	bool stkNPoints(DmpStkTrack * track) const;
	int mMinNXYhits;
	bool stkNXYhits(DmpStkTrack * track) const;
	bool stkNo1stripClusters(DmpStkTrack * track, DmpEvent * event) const;

	// PSD-STK charge consistency
	float charge_consistency(DmpStkTrack * track, DmpEvent * pev, int charge) const;

        // Track distance from the primary particle. 
	// Reconstructed using a linear model, trained separately in Python
	float mDistFromPrimCut;
	bool dist_from_prim_cut(DmpStkTrack * track, DmpEvent * event) const;

        // Proton MIP selection
	float mSTKProtonMipE;
	float mSTKXProtonMipRange;
	float mSTKYProtonMipRange;
	bool isProtonMip(DmpStkTrack * track, DmpEvent * event) const;
    };
}

#endif
