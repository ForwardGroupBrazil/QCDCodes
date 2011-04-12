#ifndef QCDEvent_h
#define QCDEvent_h
#include "KKousour/QCDAnalysis/interface/QCDTriggerObj.h"
#include "KKousour/QCDAnalysis/interface/QCDJet.h"
#include "KKousour/QCDAnalysis/interface/QCDMET.h"
#include "KKousour/QCDAnalysis/interface/QCDCaloJet.h"
#include "KKousour/QCDAnalysis/interface/QCDPFJet.h"
#include "KKousour/QCDAnalysis/interface/QCDEventHdr.h"
#include <vector>

class QCDEvent 
{
    public:
      //------------ Constructor ------------------------------
      QCDEvent() {}
      //------------ Destructor -------------------------------
      ~QCDEvent() {delete EvtHdr_; delete HLTObj_; delete L1Obj_; delete CaloJets_; delete PFJets_;}
      void setCaloMET(QCDMET *fCaloMET)                    {CaloMet_ = fCaloMET;}
      void setPFMET(QCDMET *fPFMET)                        {PFMet_ = fPFMET;}
      void setCaloJets(std::vector<QCDCaloJet> *fCaloJets) {CaloJets_ = fCaloJets;}
      void setPFJets(std::vector<QCDPFJet> *fPFJets)       {PFJets_ = fPFJets;}
      void setL1Obj(std::vector<QCDTriggerObj> *fL1Obj)    {L1Obj_ = fL1Obj;}
      void setHLTObj(std::vector<QCDTriggerObj> *fHLTObj)  {HLTObj_ = fHLTObj;}
      void setEvtHdr(QCDEventHdr *fEvtHdr)                 {EvtHdr_ = fEvtHdr;}
      int nL1Jets()                            {return L1Obj_->size();}
      int nHLTJets()                           {return HLTObj_->size();}
      int nPFJets()                            {return PFJets_->size();}
      int nCaloJets()                          {return CaloJets_->size();}
      const QCDMET& calomet()            const {return *CaloMet_;}
      const QCDMET& pfmet()              const {return *PFMet_;} 
      const QCDTriggerObj& hltobj(int i) const {return (*HLTObj_)[i];}  
      const QCDTriggerObj& l1obj(int i)  const {return (*L1Obj_)[i];}   
      const QCDPFJet& pfjet(int i)       const {return (*PFJets_)[i];}
      const QCDCaloJet& calojet(int i)   const {return (*CaloJets_)[i];}
      const QCDEventHdr& evtHdr()        const {return *EvtHdr_;}
 
    private:
      QCDEventHdr                *EvtHdr_;
      QCDMET                     *CaloMet_;
      QCDMET                     *PFMet_; 
      std::vector<QCDTriggerObj> *HLTObj_;
      std::vector<QCDTriggerObj> *L1Obj_;
      std::vector<QCDCaloJet>    *CaloJets_;
      std::vector<QCDPFJet>      *PFJets_;
};
#endif
