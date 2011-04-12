#ifndef QCDMET_h
#define QCDMET_h

class QCDMET 
{
   public:
     //------------ Constructor ------------------------------
     QCDMET() {et_=0; sumEt_=0;}
     //------------ Destructor -------------------------------
     ~QCDMET() {}
     void setVar(float fEt, float fSumEt) {et_ = fEt;sumEt_=fSumEt;} 
     float met()         {return et_;}
     float sumet()       {return sumEt_;}
     float met_o_sumet() {return et_/sumEt_;}
     
   private:
     float et_;
     float sumEt_;
};
#endif
