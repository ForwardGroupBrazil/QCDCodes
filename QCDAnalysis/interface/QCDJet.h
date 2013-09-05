#ifndef QCDJet_h
#define QCDJet_h
#include "DataFormats/JetReco/interface/Jet.h"
//-------- Generic Jet class for QCD analyses ---------------
class QCDJet 
{
   public:
     typedef reco::Particle::LorentzVector LorentzVector;
     //------------ Constructor ------------------------------
     QCDJet() {}
     //------------ Destructor -------------------------------
     ~QCDJet() {}
     //------------ Sett methods -----------------------------
     void setP4(LorentzVector fP4) {P4_ = fP4;}
     void setGen(LorentzVector fP4, float fgenR) {genP4_ = fP4;genR_ = fgenR;}
     void setCor(float fCor)                     {cor_  = fCor;} 
     void setUnc(float fUnc)                     {unc_  = fUnc;} 
     void setUncSrc(std::vector<float> fUncSrc)  {uncSrc_ = fUncSrc;}
     void setArea(float fArea)                   {area_ = fArea;}
     void setLooseID(bool fLooseID)              {looseID_ = fLooseID;} 
     void setTightID(bool fTightID)              {tightID_ = fTightID;}
     void setBtag_tche(float fbtag_tche)         {btag_tche_ = fbtag_tche;}
     void setBtag_tchp(float fbtag_tchp)         {btag_tchp_ = fbtag_tchp;}
     void setBtag_csv(float fbtag_csv)           {btag_csv_ = fbtag_csv;}
     void setBtag_ssvhe(float fbtag_ssvhe)       {btag_ssvhe_ = fbtag_ssvhe;}
     void setBtag_ssvhp(float fbtag_ssvhp)       {btag_ssvhp_ = fbtag_ssvhp;}
     void setBtag_jp(float fbtag_jp)             {btag_jp_ = fbtag_jp;}
     void setFlavor(float fFlavor)               {flavor_ = fFlavor;}
     void setBstatus(float fStatus3, float fStatus2) {status3_ = fStatus3; status2_ = fStatus2;}
     void setPartonId(float fPartonId)           {PartonId_ = fPartonId;};
     //------------ Get methods ------------------------------
     const LorentzVector& p4()    const {return P4_;}
     const LorentzVector& genp4() const {return genP4_;}
     float pt()                   const {return P4_.pt();}
     float genpt()                const {return genP4_.pt();}
     float geneta()               const {return genP4_.eta();} 
     float genR()                 const {return genR_;} 
     float ptCor()                const {return cor_ * P4_.pt();}
     float e()                    const {return P4_.energy();}
     float eCor()                 const {return cor_ * P4_.energy();}
     float eta()                  const {return P4_.eta();}
     float y()                    const {return P4_.Rapidity();}
     float phi()                  const {return P4_.phi();}
     float mass()                 const {return P4_.mass();}
     float cor()                  const {return cor_;}
     float unc()                  const {return unc_;} 
     float uncSrc(int i)          const {return uncSrc_[i];}
     float area()                 const {return area_;} 
     bool  looseID()              const {return looseID_;}
     bool  tightID()              const {return tightID_;}
     float btag_tche()            const {return btag_tche_;}
     float btag_tchp()            const {return btag_tchp_;}
     float btag_csv()             const {return btag_csv_;}  
     float btag_ssvhe()           const {return btag_ssvhe_;}
     float btag_ssvhp()           const {return btag_ssvhp_;}
     float btag_jp()              const {return btag_jp_;}
     float flavor()               const {return flavor_;}
     float bstatus3()             const {return status3_;}
     float bstatus2()             const {return status2_;}
     float PartonId()             const {return PartonId_;}
   private:
     //------ jet 4-momentum vector------------------
     LorentzVector P4_;
     //------ matched genjet 4-momentum vector-------
     LorentzVector genP4_;
     //------ matching radius -----------------------
     float genR_;
     //------ jec factor ----------------------------
     float cor_;
     //------ jec uncertainty -----------------------
     float unc_;
     //------ jec uncertainty sources ---------------
     std::vector<float> uncSrc_;
     //------ jet area ------------------------------
     float area_;
     //------ loose ID flag -------------------------
     bool  looseID_;
     //------ tight ID flag -------------------------
     bool  tightID_;
     //------ Discriminator of TCHE -----------------
     float btag_tche_;
     //------ Discriminator of TCHP -----------------
     float btag_tchp_;
     //------ Discriminator of CSV  -----------------
     float btag_csv_;
     //------ Discriminator of SSVHE ----------------
     float btag_ssvhe_;
     //------ Discriminator of SSVHP ----------------
     float btag_ssvhp_;
     //------ Discriminator of JP   -----------------
     float btag_jp_;
     //----- Flavor of Jets -------------------------
     float flavor_;
     // ---- Status of b-quark inside jet with R<0.35 --
     float status3_;
     // ---- Status of b-quark inside jet with R<0.35 --
     float status2_;
     // ---- Hard scattering parton Id --------------
     float PartonId_;

};
#endif
