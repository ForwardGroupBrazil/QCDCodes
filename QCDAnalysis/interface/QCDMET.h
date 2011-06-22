#ifndef QCDMET_h
#define QCDMET_h

class QCDMET 
{
   public:
     //------------ Constructor ------------------------------
     QCDMET() {et_=0; sumEt_=0;}
     //------------ Destructor -------------------------------
     ~QCDMET() {}
     void setVar(float fEt, float fSumEt, float fPhi) {et_ = fEt; sumEt_ = fSumEt; phi_ = fPhi;} 
     float met()         const {return et_;}
     float phi()         const {return phi_;}
     float sumet()       const {return sumEt_;}
     float met_o_sumet() const {return et_/sumEt_;}
     
   private:
     float et_;
     float sumEt_;
     float phi_;
};
#endif
