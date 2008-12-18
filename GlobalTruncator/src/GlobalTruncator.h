#ifndef UserCode_GlobalTruncator_GlobalTruncator_H
#define UserCode_GlobalTruncator_GlobalTruncator_H

/**  \class GlobalTruncator
 * 
 *   Global muon reconstructor:
 *   reconstructs muons using DT, CSC, RPC and tracker
 *   information,<BR>
 *   starting from a standalone reonstructed muon.
 *
 *
 *   $Date: 2008/05/13 03:31:44 $
 *   $Revision: 1.2 $
 *
 *   \author  R.Bellan - INFN TO
 *   \author A. Everett - Purdue University
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "UserCode/GlobalTruncator/interface/GlobalTruncRefitter.h"
#include "RecoMuon/TrackingTools/interface/MuonTrackLoader.h"

namespace edm {class ParameterSet; class Event; class EventSetup;}

class MuonTrackFinder;
class MuonServiceProxy;

class GlobalTruncator : public edm::EDProducer {

 public:

  /// constructor with config
  GlobalTruncator(const edm::ParameterSet&);
  
  /// destructor
  virtual ~GlobalTruncator(); 
  
  /// reconstruct muons
  virtual void produce(edm::Event&, const edm::EventSetup&);
  
 private:
    
  /// STA Label
  edm::InputTag theGLBCollectionLabel;

  /// the event setup proxy, it takes care the services update
  MuonServiceProxy* theService;
  
  GlobalTruncRefitter* theRefitter;

  MuonTrackLoader* theTrackLoader;
  
  std::string theAlias;
  std::vector<std::string> theRefits;
  std::vector<int> theRefitIndex;

  void setAlias( std::string alias ){
    alias.erase( alias.size() - 1, alias.size() );
    theAlias=alias;
  }
  
};

#endif
