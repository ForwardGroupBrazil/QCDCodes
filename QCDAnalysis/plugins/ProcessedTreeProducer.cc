#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <functional>
#include "TTree.h"
#include <vector>
#include <cassert>
#include <TLorentzVector.h>

#include "KKousour/QCDAnalysis/plugins/ProcessedTreeProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/JetExtendedAssociation.h"
#include "DataFormats/JetReco/interface/JetID.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

ProcessedTreeProducer::ProcessedTreeProducer(edm::ParameterSet const& cfg) 
{
  mPFJetsName        = cfg.getParameter<std::string>          ("pfjets");
  mCaloJetsName      = cfg.getParameter<std::string>          ("calojets");
  mPFJECservice      = cfg.getParameter<std::string>          ("pfjecService");
  mCaloJECservice    = cfg.getParameter<std::string>          ("calojecService");
  mPFPayloadName     = cfg.getParameter<std::string>          ("PFPayloadName");
  mCaloPayloadName   = cfg.getParameter<std::string>          ("CaloPayloadName");
  mCaloJetID         = cfg.getParameter<std::string>          ("calojetID");
  mCaloJetExtender   = cfg.getParameter<std::string>          ("calojetExtender");
  mGoodVtxNdof       = cfg.getParameter<double>               ("goodVtxNdof");
  mGoodVtxZ          = cfg.getParameter<double>               ("goodVtxZ");
  mMinCaloPt         = cfg.getParameter<double>               ("minCaloPt");
  mMinPFPt           = cfg.getParameter<double>               ("minPFPt");
  mMinNCaloJets      = cfg.getParameter<int>                  ("minNCaloJets");
  mMinNPFJets        = cfg.getParameter<int>                  ("minNPFJets");
  mIsMCarlo          = cfg.getUntrackedParameter<bool>        ("isMCarlo",false);
  mGenJetsName       = cfg.getUntrackedParameter<std::string> ("genjets","");
  processName_       = cfg.getParameter<std::string>          ("processName");
  triggerName_       = cfg.getParameter<std::string>          ("triggerName");
  triggerResultsTag_ = cfg.getParameter<edm::InputTag>        ("triggerResults");
  triggerEventTag_   = cfg.getParameter<edm::InputTag>        ("triggerEvent");
}
//////////////////////////////////////////////////////////////////////////////////////////
void ProcessedTreeProducer::beginJob() 
{
  mTree = fs->make<TTree>("ProcessedTree","ProcessedTree");
  mEvent = new QCDEvent();
  mTree->Branch("event","QCDEvent",&mEvent);
  //delete Event;
  //mEvtHdr     = new QCDEventHdr();
  //mCaloMet    = new QCDMET();
  //mPFMet      = new QCDMET();
  //mL1Objects  = new std::vector<QCDTriggerObj>(0);
  //mHLTObjects = new std::vector<QCDTriggerObj>(0);
  //mCaloJets   = new std::vector<QCDCaloJet>(0);
  //mPFJets     = new std::vector<QCDPFJet>(0);
}
//////////////////////////////////////////////////////////////////////////////////////////
void ProcessedTreeProducer::endJob() 
{/*
  delete mL1Jets;
  delete mHLTJets;
  delete mCaloJets;
  delete mPFJets; 
  delete mEvtHdr;
  delete mEvent;
  delete mCaloMet;
  delete mPFMet;
*/
}
//////////////////////////////////////////////////////////////////////////////////////////
void ProcessedTreeProducer::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{
  bool changed(true);
  if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
    if (changed) {
      // check if trigger name in (new) config
      if (triggerName_!="@") { // "@" means: analyze all triggers in config
        const unsigned int n(hltConfig_.size());
        const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName_));
        if (triggerIndex>=n) {
          cout << "ProcessedTreeProducer::analyze:"
               << " TriggerName " << triggerName_ 
               << " not available in (new) config!" << endl;
          cout << "Available TriggerNames are: " << endl;
          hltConfig_.dump("Triggers");
        }
      }
    }
  } else {
    cout << "ProcessedTreeProducer::analyze:"
         << " config extraction failure with process name "
         << processName_ << endl;
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void ProcessedTreeProducer::analyze(edm::Event const& event, edm::EventSetup const& iSetup) 
{ 
  /*
  mL1Objects->clear();
  mHLTObjects->clear();
  mCaloJets->clear();
  mPFJets->clear();
  */
  vector<LorentzVector> mL1Objects,mHLTObjects;
  vector<QCDCaloJet> mCaloJets;
  vector<QCDPFJet> mPFJets;
  QCDEventHdr mEvtHdr; 
  QCDMET mCaloMet,mPFMet;
  //-------------- Basic Event Info ------------------------------
  mEvtHdr.setRun(event.id().run());
  mEvtHdr.setEvt(event.id().event());
  mEvtHdr.setLumi(event.luminosityBlock());
  mEvtHdr.setBunch(event.bunchCrossing());
  //-------------- Trigger Info -----------------------------------
  event.getByLabel(triggerResultsTag_,triggerResultsHandle_);
  if (!triggerResultsHandle_.isValid()) {
    cout << "ProcessedTreeProducer::analyze: Error in getting TriggerResults product from Event!" << endl;
    return;
  }
  event.getByLabel(triggerEventTag_,triggerEventHandle_);
  if (!triggerEventHandle_.isValid()) {
    cout << "ProcessedTreeProducer::analyze: Error in getting TriggerEvent product from Event!" << endl;
    return;
  }
  // sanity check
  assert(triggerResultsHandle_->size()==hltConfig_.size());
  const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName_));
  assert(triggerIndex == event.triggerNames(*triggerResultsHandle_).triggerIndex(triggerName_));
  const std::pair<int,int> prescales(hltConfig_.prescaleValues(event,iSetup,triggerName_));
  mEvtHdr.setPrescales(prescales.first,prescales.second);
  // modules on this trigger path
  const vector<string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex));
  const unsigned int moduleIndex(triggerResultsHandle_->index(triggerIndex));
  bool foundL1(false);
  for (unsigned int j=0; j<=moduleIndex; ++j) {
    const string& moduleLabel(moduleLabels[j]);
    const string  moduleType(hltConfig_.moduleType(moduleLabel));
    // check whether the module is packed up in TriggerEvent product
    const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(moduleLabel,"",processName_)));
    if (filterIndex<triggerEventHandle_->sizeFilters()) {
      const Vids& VIDS (triggerEventHandle_->filterIds(filterIndex));
      const Keys& KEYS(triggerEventHandle_->filterKeys(filterIndex));
      const size_type nI(VIDS.size());
      const size_type nK(KEYS.size());
      assert(nI==nK);
      const size_type n(max(nI,nK));
      const TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());
      if (foundL1) {
        //cout<<"HLT: "<<endl; 
        for (size_type i=0; i!=n; ++i) {
          const TriggerObject& TO(TOC[KEYS[i]]);
          TLorentzVector P4;
          P4.SetPtEtaPhiM(TO.pt(),TO.eta(),TO.phi(),TO.mass());
          LorentzVector qcdhltobj(P4.Px(),P4.Py(),P4.Pz(),P4.E());
          mHLTObjects.push_back(qcdhltobj);
          //cout<<TO.pt()<<endl;
        }
      }
      if (!foundL1) { 
        //cout<<"L1: "<<endl;
        for (size_type i=0; i!=n; ++i) {
          const TriggerObject& TO(TOC[KEYS[i]]);
          TLorentzVector P4;
          P4.SetPtEtaPhiM(TO.pt(),TO.eta(),TO.phi(),TO.mass());
          LorentzVector qcdl1obj(P4.Px(),P4.Py(),P4.Pz(),P4.E());
          mL1Objects.push_back(qcdl1obj);
          //cout<<TO.pt()<<endl;  
        }
        foundL1 = true; 
      }
    }
  }
  //-------------- Vertex Info -----------------------------------
  Handle<reco::VertexCollection> recVtxs;
  event.getByLabel("offlinePrimaryVertices",recVtxs);
  int VtxGood(0);
  bool isPVgood(false);
  float PVx(0),PVy(0),PVz(0),PVndof(0);
  for(VertexCollection::const_iterator i_vtx = recVtxs->begin(); i_vtx != recVtxs->end(); i_vtx++) {
    int index = i_vtx-recVtxs->begin();
    if (index == 0) {
      PVx    = i_vtx->x();
      PVy    = i_vtx->y();
      PVz    = i_vtx->z();
      PVndof = i_vtx->ndof();
    }
    if (!(i_vtx->isFake()) && i_vtx->ndof() >= mGoodVtxNdof && fabs(i_vtx->z()) <= mGoodVtxZ) {
      if (index == 0) {
        isPVgood = true;
      }
      VtxGood++;
    }
  }
  mEvtHdr.setVertices(recVtxs->size(),VtxGood);
  mEvtHdr.setPV(isPVgood,PVndof,PVx,PVy,PVz);
  //-------------- Generator Info -------------------------------------
  Handle<GenEventInfoProduct> hEventInfo;
  if (mIsMCarlo) { 
    event.getByLabel("generator", hEventInfo);
    mEvtHdr.setPthat(hEventInfo->binningValues()[0]);
    mEvtHdr.setWeight(hEventInfo->weight());
  } 
  else {
    mEvtHdr.setPthat(0);
    mEvtHdr.setWeight(0); 
  }
  //---------------- Jets ---------------------------------------------
  mPFJEC   = JetCorrector::getJetCorrector(mPFJECservice,iSetup);
  mCALOJEC = JetCorrector::getJetCorrector(mCaloJECservice,iSetup);
  edm::ESHandle<JetCorrectorParametersCollection> PFJetCorParColl;
  if (mPFPayloadName != ""){
    iSetup.get<JetCorrectionsRecord>().get(mPFPayloadName,PFJetCorParColl); 
    JetCorrectorParameters const& PFJetCorPar = (*PFJetCorParColl)["Uncertainty"];
    mPFUnc = new JetCorrectionUncertainty(PFJetCorPar);
  }
  edm::ESHandle<JetCorrectorParametersCollection> CaloJetCorParColl;
  if (mPFPayloadName != ""){
    iSetup.get<JetCorrectionsRecord>().get(mCaloPayloadName,CaloJetCorParColl);    
    JetCorrectorParameters const& CaloJetCorPar = (*CaloJetCorParColl)["Uncertainty"];
    mCALOUnc = new JetCorrectionUncertainty(CaloJetCorPar);
  }
  Handle<GenJetCollection>  genjets;
  Handle<PFJetCollection>   pfjets;
  Handle<CaloJetCollection> calojets;
  Handle<JetExtendedAssociation::Container> calojetExtender;
  Handle<ValueMap<reco::JetID> > calojetID;
  event.getByLabel(mPFJetsName,pfjets);
  event.getByLabel(mCaloJetsName,calojets);
  event.getByLabel(mCaloJetExtender,calojetExtender);
  event.getByLabel(mCaloJetID,calojetID);
  if (mIsMCarlo)
    event.getByLabel(mGenJetsName,genjets);

  for(PFJetCollection::const_iterator i_pfjet = pfjets->begin(); i_pfjet != pfjets->end(); i_pfjet++) {
    int index = i_pfjet-pfjets->begin();
    edm::RefToBase<reco::Jet> pfjetRef(edm::Ref<PFJetCollection>(pfjets,index));
    double scale = mPFJEC->correction(*i_pfjet,pfjetRef,event,iSetup);
    double unc(0.0);
    if (mPFPayloadName != "") {
      mPFUnc->setJetEta(i_pfjet->eta());
      mPFUnc->setJetPt(scale * i_pfjet->pt());
      unc = mPFUnc->getUncertainty(true);
    }
    QCDPFJet qcdpfjet;
    qcdpfjet.setP4(i_pfjet->p4());
    qcdpfjet.setCor(scale);
    qcdpfjet.setUnc(unc);
    double chf   = i_pfjet->chargedHadronEnergyFraction();
    double nhf   = (i_pfjet->neutralHadronEnergy() + i_pfjet->HFHadronEnergy())/i_pfjet->energy();
    double phf   = i_pfjet->photonEnergyFraction();
    double elf   = i_pfjet->electronEnergyFraction();
    double chm   = i_pfjet->chargedHadronMultiplicity();
    int nhm   = i_pfjet->neutralHadronMultiplicity();
    int phm   = i_pfjet->photonMultiplicity();
    int elm   = i_pfjet->electronMultiplicity();
    int npr   = i_pfjet->chargedMultiplicity() + i_pfjet->neutralMultiplicity();
    bool looseID  = (npr>1 && phf<0.99 && nhf<0.99 && ((fabs(i_pfjet->eta())<=2.4 && elf<0.99 && chf>0 && chm>0) || fabs(i_pfjet->eta())>2.4)) ;
    qcdpfjet.setID(looseID);
    qcdpfjet.setFrac(chf,nhf,phf,elf);
    qcdpfjet.setMulti(chm,nhm,phm,elm);
    if (mIsMCarlo) {
      GenJetCollection::const_iterator i_matched;
      float rmin(999);
      for(GenJetCollection::const_iterator i_gen = genjets->begin(); i_gen != genjets->end(); i_gen++) {
        double deltaR = reco::deltaR(*i_pfjet,*i_gen);
        if (deltaR < rmin) {
          rmin = deltaR;
          i_matched = i_gen;
        }
      }
      qcdpfjet.setGen(i_matched->p4(),rmin);
    }
    else {
      LorentzVector tmpP4(0.0,0.0,0.0,0.0); 
      qcdpfjet.setGen(tmpP4,0);
    }  
    if (qcdpfjet.ptCor() >= mMinPFPt)
      mPFJets.push_back(qcdpfjet);
  }
  //----------- CaloJets -----------------------
  for(CaloJetCollection::const_iterator i_calojet = calojets->begin(); i_calojet != calojets->end(); i_calojet++) {
    int index = i_calojet-calojets->begin();
    edm::RefToBase<reco::Jet> calojetRef(edm::Ref<CaloJetCollection>(calojets,index));
    double scale = mCALOJEC->correction(*i_calojet,calojetRef,event,iSetup);
    double unc(0.0);
    if (mCaloPayloadName != "") {
      mCALOUnc->setJetEta(i_calojet->eta());
      mCALOUnc->setJetPt(scale * i_calojet->pt());
      unc = mCALOUnc->getUncertainty(true);
    } 
    QCDCaloJet qcdcalojet;
    qcdcalojet.setP4(i_calojet->p4());
    qcdcalojet.setCor(scale);
    qcdcalojet.setUnc(unc);
    double emf      = i_calojet->emEnergyFraction();
    int n90hits  = int((*calojetID)[calojetRef].n90Hits);
    double fHPD     = (*calojetID)[calojetRef].fHPD;
    double fRBX     = (*calojetID)[calojetRef].fRBX;
    int nTrkVtx  = JetExtendedAssociation::tracksAtVertexNumber(*calojetExtender,*i_calojet);
    int nTrkCalo = JetExtendedAssociation::tracksAtCaloNumber(*calojetExtender,*i_calojet);
    qcdcalojet.setVar(emf,fHPD,fRBX,n90hits,nTrkCalo,nTrkVtx);		   
    bool looseID  = ((emf > 0.01 || fabs(i_calojet->eta()) > 2.6) && (n90hits > 1) && (fHPD < 0.98));
    qcdcalojet.setID(looseID);
    if (mIsMCarlo) {
      GenJetCollection::const_iterator i_matched;
      float rmin(999);
      for(GenJetCollection::const_iterator i_gen = genjets->begin(); i_gen != genjets->end(); i_gen++) {
        double deltaR = reco::deltaR(*i_calojet,*i_gen);
        if (deltaR < rmin) {
          rmin = deltaR;
          i_matched = i_gen;
        }
      }
      qcdcalojet.setGen(i_matched->p4(),rmin);
    }
    else {
      LorentzVector tmpP4(0.0,0.0,0.0,0.0); 
      qcdcalojet.setGen(tmpP4,0);
    }
    if (qcdcalojet.ptCor() >= mMinCaloPt)
      mCaloJets.push_back(qcdcalojet);
  }
    
  //---------------- met ---------------------------------------------
  Handle<PFMETCollection> pfmet;
  Handle<CaloMETCollection> calomet;
  event.getByLabel("pfMet",pfmet);
  event.getByLabel("met",calomet);
  mPFMet.setVar((*pfmet)[0].et(),(*pfmet)[0].sumEt());
  mCaloMet.setVar((*calomet)[0].et(),(*calomet)[0].sumEt());
  //-------------- fill the tree -------------------------------------  
  sort(mCaloJets.begin(),mCaloJets.end(),sort_calojets);
  sort(mPFJets.begin(),mPFJets.end(),sort_pfjets);
  mEvent->setEvtHdr(mEvtHdr);
  mEvent->setCaloJets(mCaloJets);
  mEvent->setPFJets(mPFJets);
  mEvent->setCaloMET(mCaloMet);
  mEvent->setPFMET(mPFMet);
  mEvent->setL1Obj(mL1Objects);
  mEvent->setHLTObj(mHLTObjects);
  if ((mEvent->nPFJets() >= (unsigned)mMinNPFJets) && (mEvent->nCaloJets() >= (unsigned)mMinNCaloJets)) {
    mTree->Fill();
  }
  if (mPFPayloadName != "")
    delete mPFUnc;
  if (mCaloPayloadName != "")
    delete mCALOUnc;
}
//////////////////////////////////////////////////////////////////////////////////////////
ProcessedTreeProducer::~ProcessedTreeProducer() 
{
}

DEFINE_FWK_MODULE(ProcessedTreeProducer);
