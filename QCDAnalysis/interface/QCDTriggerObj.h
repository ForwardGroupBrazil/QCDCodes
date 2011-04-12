#ifndef QCDTriggerObj_h
#define QCDTriggerObj_h

class QCDTriggerObj 
{
   public:
     //------------ Constructor ------------------------------
     QCDTriggerObj() {}
     //------------ Destructor -------------------------------
     ~QCDTriggerObj() {}
     void setPt(float fPt)   {pt_  = fPt;}
     void setEta(float fEta) {eta_ = fEta;}
     void setPhi(float fPhi) {phi_ = fPhi;}
     float pt()      {return pt_;}
     float eta()     {return eta_;}
     float phi()     {return phi_;}

   private:
     float pt_,eta_,phi_; 
};
#endif
