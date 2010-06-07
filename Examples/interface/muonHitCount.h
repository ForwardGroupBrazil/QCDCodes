#ifndef UserCode_Examples_interface_muonHitCount_h
#define UserCode_Examples_interface_muonHitCount_h
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
namespace muon {
    /** Count the number of muon hits in a station from a given track
        subdet = 0 for all subdetectors, or 1/2/3 for DT/CSC/RPC (as defined in MuonSubdetId
	station = 0 for all stations, or 1,2,3,4 for stations
        validOnly = true (default) to include only valid hits, = false to include also bad ones 
        note: requires the TrackExtra and the TrackingRecHits to be available.
      */
  int muonHitCount(const reco::TrackRef track, int subdet=0, int station=0, bool validOnly=true) ;
}
#endif
