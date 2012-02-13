// system include files
#include <memory>
#include <iostream>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"

#include "TMath.h"

class JetExtendedProducer : public edm::EDProducer {
   public:
      explicit JetExtendedProducer(const edm::ParameterSet&);
      ~JetExtendedProducer();

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data --------------------------
      edm::InputTag src_;
      std::string name_,payload_;
      bool isJecUncSet_;
      JetCorrectionUncertainty *jecUnc_;
};

JetExtendedProducer::JetExtendedProducer(const edm::ParameterSet& iConfig)
{
  src_     = iConfig.getParameter<edm::InputTag> ("jets");
  name_    = iConfig.getParameter<std::string>   ("result");
  payload_ = iConfig.getParameter<std::string>   ("payload");

  produces<pat::JetCollection>(name_);
  isJecUncSet_ = false;
}

JetExtendedProducer::~JetExtendedProducer()
{

}

// ------------ method called to produce the data  ------------
void JetExtendedProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace reco;
  
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByLabel(src_,jets);
  edm::View<pat::Jet> pat_jets = *jets;  

  std::auto_ptr<pat::JetCollection> result (new pat::JetCollection); //Extended jets
  //-------------- Vertex Info -----------------------------------
  Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel("goodOfflinePrimaryVertices",recVtxs);

  //-------------- Set the JEC uncertainty object ----------------
  if (!isJecUncSet_ && payload_ != "") {
    edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
    iSetup.get<JetCorrectionsRecord>().get(payload_,JetCorParColl); 
    JetCorrectorParameters const& par = (*JetCorParColl)["Uncertainty"];
    jecUnc_ = new JetCorrectionUncertainty(par);
    isJecUncSet_ = true;
  }
    
  unsigned int jetIdx = 0;
  for(edm::View<pat::Jet>::const_iterator ijet = pat_jets.begin(); ijet != pat_jets.end(); ++ijet) {
    //----- calculations based on the constituents -----
    std::vector<PFCandidatePtr> pfConst(ijet->getPFConstituents());
    int n_pf = pfConst.size();                                                                                          
    float phiJet = ijet->phi();                                                                                         
    float etaJet = ijet->eta();                                                                                         
    float deta,dphi,dR,weight,sum(0.0),sum2(0.0),sum_deta(0.0),sum_dphi(0.0),sum_deta2(0.0),sum_dphi2(0.0),sum_detadphi(0.0);                                                   
    // variables for the jet thrust computation               
    float weightT, Teta(0), Tphi(0);                                                                                                   
    float sumT=0; 
    float Ttheta = -9;                                                                                                 

    float jetPtMax = 0;                                                                                                           
    for(int j = 0; j < n_pf; j++) {                                                                                 
      deta = pfConst[j]->eta() - etaJet;                                                                                               
      dphi = 2*atan(tan((pfConst[j]->phi()-phiJet)/2));                                                                              
      dR = sqrt(deta*deta + dphi*dphi);                                                                                              
      weightT = pfConst[j]->pt();                                                                                                    
      sumT += weightT;                                                                                                               
      Teta += weightT * dR * deta;                                                                                                   
      Tphi += weightT * dR * dphi;                                                                                                   

      weight = weightT*weightT;                                                                                                      
      sum += weight;                                     
      sum2 += weight*weight;                                                                                                         
      sum_deta += deta*weight;                                                                                                       
      sum_dphi += dphi*weight;                                                                                                       
      sum_deta2 += deta*deta*weight;                                                                                                 
      sum_detadphi += deta*dphi*weight;                                                                                              
      sum_dphi2 += dphi*dphi*weight;                                                                                                 
      if (fabs(pfConst[j]->charge()) > 0) {                                                                                            
         //PFJetTracks[JetsPF]++;                                                                                                     
         jetPtMax = TMath::Max(jetPtMax,float(pfConst[j]->pt()));         
      }                                                                                                                              
    }                                                                                                                                
    Teta = Teta/sumT;                                                                                                                
    Tphi = Tphi/sumT;                                                                                                                
    if (Teta != 0 && Tphi !=0 ) {
      Ttheta = atan2(Tphi,Teta);                                                                           
    } 
    float ave_deta = sum_deta/sum;                                                                                                   
    float ave_dphi = sum_dphi/sum;                                                                                    
    float ave_deta2 = sum_deta2/sum;                                                                                                 
    float ave_dphi2 = sum_dphi2/sum;                                                                                                 

    float a =ave_deta2-ave_deta*ave_deta;                                                                                            
    float b =ave_dphi2-ave_dphi*ave_dphi;                                                                                            
    float c =-(sum_detadphi/sum-ave_deta*ave_dphi);                                                                                  
    float delta = sqrt(fabs((a-b)*(a-b)+4*c*c));                                                                                     

    float axis1 = sqrt((a+b+delta)/2);                                                                                               
    float axis2 = sqrt((a+b-delta)/2);                                                                                               

    float tana = (b-a+delta)/2/c;
    
    float ptD = sqrt(sum2)/sum;                                                                                                

    //---- vertex association -----------
    //---- get the vector of tracks -----
    const reco::PFJet& pfJet = dynamic_cast <const reco::PFJet&> (*(ijet->originalObject()));
    reco::TrackRefVector vTrks(pfJet.getTrackRefs());
    float sumTrkPt(0.0),sumTrkPtBeta(0.0),beta(0.0);
    //---- loop over the tracks of the jet ----
    for(reco::TrackRefVector::const_iterator i_trk = vTrks.begin(); i_trk != vTrks.end(); i_trk++) {
      sumTrkPt += (*i_trk)->pt();
      //---- loop over all vertices ----------------------------
      for(unsigned ivtx = 0;ivtx < recVtxs->size();ivtx++) {
        //---- loop over the tracks associated with the vertex ---
        for(reco::Vertex::trackRef_iterator i_vtxTrk = (*recVtxs)[ivtx].tracks_begin(); i_vtxTrk != (*recVtxs)[ivtx].tracks_end(); ++i_vtxTrk) {
          //---- match the jet track to the track from the vertex ----
          reco::TrackRef trkRef(i_vtxTrk->castTo<reco::TrackRef>());
          //---- check if the tracks match -------------------------
          if (trkRef == (*i_trk)) {
            if (ivtx > 0) {
              sumTrkPtBeta += (*i_trk)->pt();
            }   
            break;
          }
        }
      } 
    }
    if (sumTrkPt > 0) {
      beta = 1.-sumTrkPtBeta/sumTrkPt;  
    }
    pat::Jet * extendedJet = ijet->clone();
    jecUnc_->setJetEta(ijet->eta());
    jecUnc_->setJetPt(ijet->pt());
    extendedJet->addUserFloat("jecUnc",jecUnc_->getUncertainty(true));
    extendedJet->addUserFloat("beta",beta);
    extendedJet->addUserFloat("ptMax",jetPtMax);
    extendedJet->addUserFloat("ptD",ptD);
    extendedJet->addUserFloat("axis1",axis1);
    extendedJet->addUserFloat("axis2",axis2);
    extendedJet->addUserFloat("tana",tana);
    extendedJet->addUserFloat("ttheta",Ttheta);

    result->push_back(*extendedJet); 
    ++jetIdx;
  }

  iEvent.put(result,name_);
}

// ------------ method called once each job just before starting event loop  ------------
void JetExtendedProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void JetExtendedProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(JetExtendedProducer);
