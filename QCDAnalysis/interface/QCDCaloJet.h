#ifndef QCDCaloJet_h
#define QCDCaloJet_h
#include "KKousour/QCDAnalysis/interface/QCDJet.h"

class QCDCaloJet : public QCDJet {
   public:
     //------------ Constructor ------------------------------
     QCDCaloJet() {emf_=0;fHPD_=0;fRBX_=0;n90hits_=0;nTrkCalo_=0;nTrkVtx_=0;}
     //------------ Destructor -------------------------------
     ~QCDCaloJet() {};
     void setVar(float femf,float ffHPD,float ffRBX,int fn90,int fnTrkCalo,int fnTrkVtx)  {emf_ =
     femf;fHPD_=ffHPD;fRBX_=ffRBX;n90hits_=fn90;nTrkCalo_=fnTrkCalo;nTrkVtx_=fnTrkVtx;}
     float emf()      {return emf_;} 
     float fHPD()     {return fHPD_;}
     float fRBX()     {return fRBX_;}
     int n90hits()    {return n90hits_;}
     int nTrkCalo()   {return nTrkCalo_;} 
     int nTrkVtx()    {return nTrkVtx_;}
     
   private:
     float emf_;
     float fHPD_;
     float fRBX_;
     int n90hits_;
     int nTrkCalo_;
     int nTrkVtx_;
    };
#endif    
