#ifndef QCDEvent_h
#define QCDEvent_h
#include "KKousour/QCDAnalysis/interface/QCDTriggerObj.h"
#include "KKousour/QCDAnalysis/interface/QCDJet.h"
#include "KKousour/QCDAnalysis/interface/QCDMET.h"
#include "KKousour/QCDAnalysis/interface/QCDCaloJet.h"
#include "KKousour/QCDAnalysis/interface/QCDPFJet.h"
#include "KKousour/QCDAnalysis/interface/QCDEventHdr.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include <vector>

class QCDEvent 
{
    public:
      typedef reco::Particle::LorentzVector LorentzVector;
      //------------ Constructor ------------------------------
      QCDEvent();
      //------------ Destructor -------------------------------
      ~QCDEvent();
      void setCaloMET(const QCDMET& fCaloMET)                     {CaloMet_ = fCaloMET;}
      void setPFMET(const QCDMET& fPFMET)                         {PFMet_ = fPFMET;}
      void setEvtHdr(const QCDEventHdr& fEvtHdr)                  {EvtHdr_ = fEvtHdr;}
      void setCaloJets(const std::vector<QCDCaloJet>& fCaloJets);
      void setPFJets(const std::vector<QCDPFJet>& fPFJets);
      void setL1Obj(const std::vector<QCDTriggerObj>& fL1Obj);
      void setHLTObj(const std::vector<QCDTriggerObj>& fHLTObj);
      unsigned int nL1Obj()                            const {return L1Obj_.size();}
      unsigned int nHLTObj()                           const {return HLTObj_.size();}
      unsigned int nPFJets()                           const {return PFJets_.size();}
      unsigned int nCaloJets()                         const {return CaloJets_.size();}
      float pfmjj();
      float calomjj();
      float pfmjjcor(int k);
      float calomjjcor(int k);
      float pfmjjgen();
      float calomjjgen();
      const QCDMET&        calomet()      const {return CaloMet_;}
      const QCDMET&        pfmet()        const {return PFMet_;} 
      const QCDTriggerObj& hltobj(int i)  const {return HLTObj_[i];}  
      const QCDTriggerObj& l1obj(int i)   const {return L1Obj_[i];}   
      const QCDPFJet&      pfjet(int i)   const {return PFJets_[i];}
      const QCDCaloJet&    calojet(int i) const {return CaloJets_[i];}
      const QCDEventHdr&   evtHdr()       const {return EvtHdr_;}
 
    private:
      QCDEventHdr                EvtHdr_;
      QCDMET                     CaloMet_;
      QCDMET                     PFMet_; 
      std::vector<QCDTriggerObj> HLTObj_;
      std::vector<QCDTriggerObj> L1Obj_;
      std::vector<QCDCaloJet>    CaloJets_;
      std::vector<QCDPFJet>      PFJets_;
};
#endif
