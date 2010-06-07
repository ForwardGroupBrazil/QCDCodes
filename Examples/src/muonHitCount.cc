#include "UserCode/Examples/interface/muonHitCount.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/DTLayerId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"

int muon::muonHitCount(const reco::TrackRef track, int subdet, int station, bool validOnly) {

  int hits = 0;
  //std::cout <<"Sono in MuonHitCountUtility " << subdet << " " << station << " " << validOnly << std::endl;
  if (track.isNull()) return 0;
  for (trackingRecHit_iterator it = track->recHitsBegin(), ed = track->recHitsEnd(); it != ed; ++it) {
    DetId id = (*it)->geographicalId();
    if (id.det() != DetId::Muon) continue;
    if (subdet != 0 && id.subdetId() != subdet) continue;
    if (validOnly && !(*it)->isValid()) continue;
    switch (id.subdetId()) {
    case MuonSubdetId::DT: 
      if(station == 0 || DTLayerId(id).station() == station) hits++; 
      break;
    case MuonSubdetId::CSC:
      if(station == 0 || CSCDetId(id).station()  == station) hits++; 
      break;
    case MuonSubdetId::RPC: 
      if(station == 0 || RPCDetId(id).station()  == station) hits++; 
      break;
    }
  }
  return hits;
}

