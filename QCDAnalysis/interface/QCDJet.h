#ifndef QCDJet_h
#define QCDJet_h
#include "DataFormats/JetReco/interface/Jet.h"

class QCDJet 
{
   public:
     typedef reco::Particle::LorentzVector LorentzVector;
     //------------ Constructor ------------------------------
     QCDJet() {}
     //------------ Destructor -------------------------------
     ~QCDJet() {}
     void setP4(LorentzVector fP4) {P4_ = fP4;}
     void setGen(LorentzVector fP4, float fgenR) {genP4_ = fP4;genR_ = fgenR;}
     void setCor(float fCor)       {cor_ = fCor;} 
     void setUnc(float fUnc)       {unc_ = fUnc;} 
     void setID(bool fID)          {id_ = fID;} 
     const LorentzVector& p4() const {return P4_;}
     const LorentzVector& genp4() const {return genP4_;}
     float pt()      const {return P4_.pt();}
     float genpt()   const {return genP4_.pt();}
     float geneta()  const {return genP4_.eta();} 
     float genR()    const {return genR_;} 
     float ptCor()   const {return cor_ * P4_.pt();}
     float e()       const {return P4_.energy();}
     float eCor()    const {return cor_ * P4_.energy();}
     float eta()     const {return P4_.eta();}
     float y()       const {return P4_.y();}
     float phi()     const {return P4_.phi();}
     float mass()    const {return P4_.mass();}
     float cor()     const {return cor_;}
     float unc()     const {return unc_;} 
     bool   id()     const {return id_;}

   private:
     LorentzVector P4_,genP4_;
     float genR_;
     float cor_;
     float unc_;
     bool  id_;
};
#endif
