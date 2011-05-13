#include "KKousour/QCDAnalysis/interface/QCDEvent.h"
//---------------------------------------------------
QCDEvent::QCDEvent()
{  
  L1Obj_.clear();
  CaloJets_.clear();
  PFJets_.clear(); 
  GenJets_.clear();
}
//---------------------------------------------------
QCDEvent::~QCDEvent()
{
}
//---------------------------------------------------
void QCDEvent::setCaloJets(const std::vector<QCDCaloJet>& fCaloJets) 
{ 
  CaloJets_.clear();
  for(unsigned i=0;i<fCaloJets.size();i++) {
    CaloJets_.push_back(fCaloJets[i]);
  }
}
//---------------------------------------------------
void QCDEvent::setPFJets(const std::vector<QCDPFJet>& fPFJets) 
{ 
  PFJets_.clear();
  for(unsigned i=0;i<fPFJets.size();i++) {
    PFJets_.push_back(fPFJets[i]);
  }
}
//---------------------------------------------------
void QCDEvent::setGenJets(const std::vector<LorentzVector>& fGenJets)
{
  GenJets_.clear();
  for(unsigned i=0;i<fGenJets.size();i++) {
    GenJets_.push_back(fGenJets[i]);
  }
}
//---------------------------------------------------
void QCDEvent::setL1Obj(const std::vector<std::vector<LorentzVector> >& fL1Obj)       
{
  L1Obj_.clear();
  for(unsigned i=0;i<fL1Obj.size();i++) {
    std::vector<LorentzVector> vv;
    for(unsigned j=0;j<fL1Obj[i].size();j++) {
      vv.push_back(fL1Obj[i][j]);
    }
    L1Obj_.push_back(vv);
  }
}
//---------------------------------------------------
void QCDEvent::setHLTObj(const std::vector<std::vector<LorentzVector> >& fHLTObj)       
{
  HLTObj_.clear();
  for(unsigned i=0;i<fHLTObj.size();i++) {
    std::vector<LorentzVector> vv;
    for(unsigned j=0;j<fHLTObj[i].size();j++) {
      vv.push_back(fHLTObj[i][j]);
    }
    HLTObj_.push_back(vv);
  }
}
//---------------------------------------------------
float QCDEvent::genmjj()
{
  if (GenJets_.size() < 2)
    return 0.0;
  else {
    return (GenJets_[0]+GenJets_[1]).mass();
  }
}
//---------------------------------------------------
float QCDEvent::pfmjj()
{
  if (PFJets_.size() < 2)
    return 0.0;
  else {
    const LorentzVector& P0 = PFJets_[0].p4();
    const LorentzVector& P1 = PFJets_[1].p4();
    return (P0+P1).mass();
  }
}
//---------------------------------------------------
float QCDEvent::pfmjjcor(int k)
{
  int sign(0);
  if (PFJets_.size() < 2)
    return 0.0;
  else {
    if (k>0)
      sign = 1;
    if (k<0)
      sign = -1;  
    const LorentzVector& P0 = PFJets_[0].p4();
    const LorentzVector& P1 = PFJets_[1].p4();
    double cor0 = PFJets_[0].cor();
    double cor1 = PFJets_[1].cor();
    double unc0 = PFJets_[0].unc();
    double unc1 = PFJets_[1].unc();
    return (cor0*(1+sign*unc0)*P0+cor1*(1+sign*unc1)*P1).mass();
  }
}
//---------------------------------------------------
float QCDEvent::calomjj()
{
  if (CaloJets_.size() < 2)
    return 0.0;
  else {
    const LorentzVector& P0 = CaloJets_[0].p4();
    const LorentzVector& P1 = CaloJets_[1].p4();
    return (P0+P1).mass();
  }
}
//---------------------------------------------------
float QCDEvent::calomjjcor(int k)
{
  int sign(0);
  if (CaloJets_.size() < 2)
    return 0.0;
  else {
    if (k>0)
      sign = 1;
    if (k<0)
      sign = -1;
    const LorentzVector& P0 = CaloJets_[0].p4();
    const LorentzVector& P1 = CaloJets_[1].p4();
    double cor0 = CaloJets_[0].cor();
    double cor1 = CaloJets_[1].cor();
    double unc0 = CaloJets_[0].unc();
    double unc1 = CaloJets_[1].unc();
    return (cor0*(1+sign*unc0)*P0+cor1*(1+sign*unc1)*P1).mass();
  }
}
//---------------------------------------------------
float QCDEvent::pfmjjgen()
{
  if (PFJets_.size() < 2)
    return 0.0;
  else {
    const LorentzVector& P0 = PFJets_[0].genp4();
    const LorentzVector& P1 = PFJets_[1].genp4();
    return (P0+P1).mass();
  }
}
//---------------------------------------------------
float QCDEvent::calomjjgen()
{
  if (CaloJets_.size() < 2)
    return 0.0;
  else {
    const LorentzVector& P0 = CaloJets_[0].genp4();
    const LorentzVector& P1 = CaloJets_[1].genp4();
    return (P0+P1).mass();
  }
}


