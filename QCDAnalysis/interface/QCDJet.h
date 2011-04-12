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
     float pt()      {return P4_.pt();}
     float genpt()   {return genP4_.pt();}
     float geneta()  {return genP4_.eta();} 
     float genR()    {return genR_;} 
     float ptCor()   {return cor_ * P4_.pt();}
     float e()       {return P4_.energy();}
     float eCor()    {return cor_ * P4_.energy();}
     float eta()     {return P4_.eta();}
     float y()       {return P4_.y();}
     float phi()     {return P4_.phi();}
     float mass()    {return P4_.mass();}
     float cor()     {return cor_;}
     float unc()     {return unc_;} 
     bool   id()     {return id_;}

   private:
     LorentzVector P4_,genP4_;
     float genR_;
     float cor_;
     float unc_;
     bool  id_;
};
#endif
