#ifndef QCDPFJet_h
#define QCDPFJet_h
#include "KKousour/QCDAnalysis/interface/QCDJet.h"

class QCDPFJet : public QCDJet {
   public:
     //------------ Constructor ------------------------------
     QCDPFJet() {chf_=0;nhf_=0;phf_=0;elf_=0;muf_=0;chm_=0;nhm_=0;phm_=0;elm_=0;mum_=0;}
     //------------ Destructor -------------------------------
     ~QCDPFJet() {}
     //------------ Set methods ------------------------------
     void setFrac(float fchf, float fnhf, float fphf, float felf, float fmuf)  {chf_ = fchf; nhf_ = fnhf; phf_ = fphf; elf_ = felf; muf_ = fmuf;}
     void setMulti(int fncand, int fchm, int fnhm, int fphm, int felm, int fmum) {ncand_ = fncand; chm_ = fchm; nhm_ = fnhm; phm_ = fphm; elm_ = felm; mum_ = fmum;}
     //------------ Get methods ------------------------------ 
     float chf() const {return chf_;} 
     float nhf() const {return nhf_;}
     float phf() const {return phf_;} 
     float elf() const {return elf_;}
     float muf() const {return muf_;}
     int chm()   const {return chm_;}
     int nhm()   const {return nhm_;}
     int phm()   const {return phm_;}
     int elm()   const {return elm_;}
     int mum()   const {return mum_;}
     int ncand() const {return ncand_;}
   private:
     //---- charged hadron energy fraction ----
     float chf_;
     //---- neutral hadron energy fraction ----
     float nhf_;
     //---- photon energy fraction ------------
     float phf_;
     //---- electron energy fraction ----------
     float elf_;
     //---- muon energy fraction --------------
     float muf_;
     //---- charged hadron multiplicity -------
     int chm_;
     //---- neutral hadron multiplicity -------
     int nhm_;
     //---- photon multiplicity ---------------
     int phm_;
     //---- electron multiplicity -------------
     int elm_;
     //---- muon multiplicity -----------------
     int mum_;
     //---- number of PF candidates -----------
     int ncand_;
    };
#endif    
