#ifndef QCDEvent_h
#define QCDEvent_h
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
      void setGenJets(const std::vector<LorentzVector>& fGenJets);
      void setL1Obj(const std::vector<std::vector<LorentzVector> >& fL1Obj);
      void setHLTObj(const std::vector<std::vector<LorentzVector> >& fHLTObj);
      void setPrescales(const std::vector<int>& fPreL1, const std::vector<int>& fPreHLT) {L1Prescale_ = fPreL1; HLTPrescale_ = fPreHLT;}
      void setTrigDecision(const std::vector<int>& fTrigDecision) {TriggerDecision_ = fTrigDecision;}                           
      unsigned int nTriggers()                         const {return TriggerDecision_.size();}
      unsigned int nL1Obj(int i)                       const {return L1Obj_[i].size();}
      unsigned int nHLTObj(int i)                      const {return HLTObj_[i].size();}
      unsigned int nPFJets()                           const {return PFJets_.size();}
      unsigned int nCaloJets()                         const {return CaloJets_.size();}
      unsigned int nGenJets()                          const {return GenJets_.size();}
      int fired(int i)                                 const {return TriggerDecision_[i];}
      int preL1(int i)                                 const {return L1Prescale_[i];}
      int preHLT(int i)                                const {return HLTPrescale_[i];}
      float pfmjj();
      float calomjj();
      float genmjj(); 
      float pfmjjcor(int k);
      float calomjjcor(int k);
      float pfmjjgen();
      float calomjjgen();
      const QCDMET&        calomet()                   const {return CaloMet_;}
      const QCDMET&        pfmet()                     const {return PFMet_;} 
      const LorentzVector& hltobj(int itrig, int iobj) const {return (HLTObj_[itrig])[iobj];}  
      const LorentzVector& l1obj(int itrig, int iobj)  const {return (L1Obj_[itrig])[iobj];}   
      const LorentzVector& genjet(int i)               const {return GenJets_[i];}
      const QCDPFJet&      pfjet(int i)                const {return PFJets_[i];}
      const QCDCaloJet&    calojet(int i)              const {return CaloJets_[i];}
      const QCDEventHdr&   evtHdr()                    const {return EvtHdr_;}
 
    private:
     
      QCDEventHdr                              EvtHdr_;
      QCDMET                                   CaloMet_;
      QCDMET                                   PFMet_; 
      std::vector<int>                         TriggerDecision_;
      std::vector<int>                         L1Prescale_;
      std::vector<int>                         HLTPrescale_;
      std::vector<std::vector<LorentzVector> > HLTObj_;
      std::vector<std::vector<LorentzVector> > L1Obj_;
      std::vector<LorentzVector>               GenJets_;
      std::vector<QCDCaloJet>                  CaloJets_;
      std::vector<QCDPFJet>                    PFJets_;
};
#endif
