#ifndef QCDPFJet_h
#define QCDPFJet_h
#include "KKousour/QCDAnalysis/interface/QCDJet.h"

class QCDPFJet : public QCDJet {
   public:
     //------------ Constructor ------------------------------
     QCDPFJet() {chf_=0;nhf_=0;phf_=0;elf_=0;chm_=0;nhm_=0;phm_=0;elm_=0;}
     //------------ Destructor -------------------------------
     ~QCDPFJet() {};
     void setFrac(float fchf, float fnhf, float fphf, float felf)  {chf_ = fchf; nhf_ = fnhf; phf_ = fphf; elf_ = felf;}
     void setMulti(int fchm, int fnhm, int fphm, int felm) {chm_ = fchm; nhm_ = fnhm; phm_ = fphm; elm_ = felm;}
     float chf() const {return chf_;} 
     float nhf() const {return nhf_;}
     float phf() const {return phf_;} 
     float elf() const {return elf_;}
     int chm()   const {return chm_;}
     int nhm()   const {return nhm_;}
     int phm()   const {return phm_;}
     int elm()   const {return elm_;}
   private:
     float chf_;
     float nhf_;
     float phf_;
     float elf_;
     int chm_;
     int nhm_;
     int phm_;
     int elm_;
    };
#endif    
