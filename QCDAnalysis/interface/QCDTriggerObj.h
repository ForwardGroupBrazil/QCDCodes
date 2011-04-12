#ifndef QCDTriggerObj_h
#define QCDTriggerObj_h

class QCDTriggerObj 
{
   public:
     //------------ Constructor ------------------------------
     QCDTriggerObj();
     QCDTriggerObj(float fPt, float fEta, float fPhi) {pt_  = fPt; eta_ = fEta; phi_ = fPhi;}
     //------------ Destructor -------------------------------
     ~QCDTriggerObj();
     float pt()    const {return pt_;}
     float eta()   const {return eta_;}
     float phi()   const {return phi_;}

   private:
     float pt_;
     float eta_;
     float phi_; 
};
#endif
