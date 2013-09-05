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
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"


PFJetIDSelectionFunctor pfJetIDLoose( PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE );
PFJetIDSelectionFunctor pfJetIDTight( PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT );
pat::strbitset retpf = pfJetIDLoose.getBitTemplate();

JetIDSelectionFunctor jetIDLoose( JetIDSelectionFunctor::PURE09, JetIDSelectionFunctor::LOOSE );
JetIDSelectionFunctor jetIDTight( JetIDSelectionFunctor::PURE09, JetIDSelectionFunctor::TIGHT );
pat::strbitset ret = jetIDLoose.getBitTemplate();

ProcessedTreeProducer::ProcessedTreeProducer(edm::ParameterSet const& cfg)
{
  mPFPayloadName     = cfg.getParameter<std::string>               ("PFPayloadName");
  mCaloPayloadName   = cfg.getParameter<std::string>               ("CaloPayloadName");
  mGoodVtxNdof       = cfg.getParameter<double>                    ("goodVtxNdof");
  mGoodVtxZ          = cfg.getParameter<double>                    ("goodVtxZ");
  mMinCaloPt         = cfg.getParameter<double>                    ("minCaloPt");
  mMinPFPt           = cfg.getParameter<double>                    ("minPFPt");
  mMinPFFatPt        = cfg.getParameter<double>                    ("minPFFatPt");
  mMaxPFFatEta       = cfg.getParameter<double>                    ("maxPFFatEta");
  mMinJJMass         = cfg.getParameter<double>                    ("minJJMass");
  mMaxY              = cfg.getParameter<double>                    ("maxY");
  mMinNCaloJets      = cfg.getParameter<int>                       ("minNCaloJets");
  mMinNPFJets        = cfg.getParameter<int>                       ("minNPFJets");
  mCaloJetExtender   = cfg.getParameter<edm::InputTag>             ("calojetExtender");
  mOfflineVertices   = cfg.getParameter<edm::InputTag>             ("offlineVertices");
  mPFJetsName        = cfg.getParameter<edm::InputTag>             ("pfjets");
  mCaloJetsName      = cfg.getParameter<edm::InputTag>             ("calojets");
  mPFMETName         = cfg.getParameter<edm::InputTag>             ("pfmet");
  mCaloMETName       = cfg.getParameter<edm::InputTag>             ("calomet");
  mSrcCaloRho        = cfg.getParameter<edm::InputTag>             ("srcCaloRho");
  mSrcPFRho          = cfg.getParameter<edm::InputTag>             ("srcPFRho");
  mSrcPU             = cfg.getUntrackedParameter<edm::InputTag>    ("srcPU",edm::InputTag("addPileupInfo"));
  mGenJetsName       = cfg.getUntrackedParameter<edm::InputTag>    ("genjets",edm::InputTag(""));
  mPrintTriggerMenu  = cfg.getUntrackedParameter<bool>             ("printTriggerMenu",false);
  mIsMCarlo          = false;
  mUseGenInfo        = cfg.getUntrackedParameter<bool>             ("useGenInfo",false);
  mMinGenPt          = cfg.getUntrackedParameter<double>           ("minGenPt",30);
  processName_       = cfg.getParameter<std::string>               ("processName");
  triggerNames_      = cfg.getParameter<std::vector<std::string> > ("triggerName");
  triggerResultsTag_ = cfg.getParameter<edm::InputTag>             ("triggerResults");
  triggerEventTag_   = cfg.getParameter<edm::InputTag>             ("triggerEvent");
  mPFJECUncSrc       = cfg.getParameter<std::string>               ("jecUncSrc");
  mPFJECUncSrcNames  = cfg.getParameter<std::vector<std::string> > ("jecUncSrcNames");
  mXsec              = cfg.getUntrackedParameter<double>           ("Xsec",0.);
}
//////////////////////////////////////////////////////////////////////////////////////////
void ProcessedTreeProducer::beginJob()
{
  mTree = fs->make<TTree>("ProcessedTree","ProcessedTree");
  mEvent = new QCDEvent();
  mTree->Branch("events","QCDEvent",&mEvent);
  mTriggerNamesHisto = fs->make<TH1F>("TriggerNames","TriggerNames",1,0,1);
  mTriggerNamesHisto->SetBit(TH1::kCanRebin);
  for(unsigned i=0;i<triggerNames_.size();i++)
    mTriggerNamesHisto->Fill(triggerNames_[i].c_str(),1);
  mTriggerPassHisto = fs->make<TH1F>("TriggerPass","TriggerPass",1,0,1);
  mTriggerPassHisto->SetBit(TH1::kCanRebin);
  isPFJecUncSet_ = false;
  isCaloJecUncSet_ = false;
}
//////////////////////////////////////////////////////////////////////////////////////////
void ProcessedTreeProducer::endJob()
{
}
//////////////////////////////////////////////////////////////////////////////////////////
void ProcessedTreeProducer::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{
  bool changed(true);
  if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
    if (changed) {
      // check if trigger names in (new) config
      cout<<"New trigger menu found !!!"<<endl;
      triggerIndex_.clear();
      const unsigned int n(hltConfig_.size());
      for(unsigned itrig=0;itrig<triggerNames_.size();itrig++) {
        triggerIndex_.push_back(hltConfig_.triggerIndex(triggerNames_[itrig]));
        cout<<triggerNames_[itrig]<<" "<<triggerIndex_[itrig]<<" ";
        if (triggerIndex_[itrig] >= n)
          cout<<"does not exist in the current menu"<<endl;
        else
          cout<<"exists"<<endl;
      }
      cout << "Available TriggerNames are: " << endl;
      if (mPrintTriggerMenu)
        hltConfig_.dump("Triggers");
    }
  }
  else {
    cout << "ProcessedTreeProducer::analyze:"
         << " config extraction failure with process name "
         << processName_ << endl;
  }
}

//////////////////////////////////////////////////////////////////////////////////////////
void ProcessedTreeProducer::analyze(edm::Event const& event, edm::EventSetup const& iSetup)
{
  mIsMCarlo = !event.isRealData();

  vector<QCDCaloJet>    mCaloJets;
  vector<QCDPFJet>      mPFJets;
  vector<QCDJet>        mPFFatJets;
  vector<QCDPFJet>      tmpPFJets;
  vector<LorentzVector> mGenJets;
  QCDEventHdr mEvtHdr;
  QCDMET mCaloMet,mPFMet;
  //-------------- Basic Event Info ------------------------------
  mEvtHdr.setRun(event.id().run());
  mEvtHdr.setEvt(event.id().event());
  mEvtHdr.setLumi(event.luminosityBlock());
  mEvtHdr.setBunch(event.bunchCrossing());

  //-------------- Beam Spot --------------------------------------
  edm::Handle<reco::BeamSpot> beamSpot;
  event.getByLabel("offlineBeamSpot", beamSpot);
  if (beamSpot.isValid())
    mEvtHdr.setBS(beamSpot->x0(),beamSpot->y0(),beamSpot->z0());
  else
    mEvtHdr.setBS(-999,-999,-999);
  //-------------- HCAL Noise Summary -----------------------------
  edm::Handle<bool> noiseSummary;
  if (!mIsMCarlo) {
    event.getByLabel(edm::InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResult"), noiseSummary);
    mEvtHdr.setHCALNoise(*noiseSummary);
  }
  else
    mEvtHdr.setHCALNoise(true);
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
  vector<int> L1Prescales,HLTPrescales,Fired;
  vector<vector<LorentzVector> > mL1Objects,mHLTObjects;
  // sanity check
  assert(triggerResultsHandle_->size() == hltConfig_.size());
  //------ loop over all trigger names ---------
  for(unsigned itrig=0;itrig<triggerNames_.size() && !mIsMCarlo;itrig++) {
    bool accept(false);
    int preL1(-1);
    int preHLT(-1);
    int tmpFired(-1);
    vector<LorentzVector> vvL1,vvHLT;
    if (triggerIndex_[itrig] < hltConfig_.size()) {
      accept = triggerResultsHandle_->accept(triggerIndex_[itrig]);
      const std::pair<int,int> prescales(hltConfig_.prescaleValues(event,iSetup,triggerNames_[itrig]));
      preL1    = prescales.first;
      preHLT   = prescales.second;
      if (!accept)
        tmpFired = 0;
      else {
        mTriggerPassHisto->Fill(triggerNames_[itrig].c_str(),1);
        tmpFired = 1;
      }
      //--------- modules on this trigger path--------------
      const vector<string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex_[itrig]));
      const unsigned int moduleIndex(triggerResultsHandle_->index(triggerIndex_[itrig]));
      bool foundL1(false);

      for(unsigned int j=0; j<=moduleIndex; ++j) {
        const string& moduleLabel(moduleLabels[j]);
        const string  moduleType(hltConfig_.moduleType(moduleLabel));
        //--------check whether the module is packed up in TriggerEvent product
        const unsigned int filterIndex(triggerEventHandle_->filterIndex(edm::InputTag(moduleLabel,"",processName_)));
        if (filterIndex<triggerEventHandle_->sizeFilters()) {
          const Vids& VIDS (triggerEventHandle_->filterIds(filterIndex));
          const Keys& KEYS(triggerEventHandle_->filterKeys(filterIndex));
          const size_type nI(VIDS.size());
          const size_type nK(KEYS.size());
          assert(nI==nK);
          const size_type n(max(nI,nK));
          const TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());
          if (foundL1) {
            for(size_type i=0; i!=n; ++i) {
              const TriggerObject& TO(TOC[KEYS[i]]);
              TLorentzVector P4;
              P4.SetPtEtaPhiM(TO.pt(),TO.eta(),TO.phi(),TO.mass());
              LorentzVector qcdhltobj(P4.Px(),P4.Py(),P4.Pz(),P4.E());
              vvHLT.push_back(qcdhltobj);
              //cout<<TO.pt()<<endl;
            }
          }
          else {
            for(size_type i=0; i!=n; ++i) {
              const TriggerObject& TO(TOC[KEYS[i]]);
              TLorentzVector P4;
              P4.SetPtEtaPhiM(TO.pt(),TO.eta(),TO.phi(),TO.mass());
              LorentzVector qcdl1obj(P4.Px(),P4.Py(),P4.Pz(),P4.E());
              vvL1.push_back(qcdl1obj);
              //cout<<TO.pt()<<endl;
            }
            foundL1 = true;
          }
        }
      }// loop over modules
    }// if the trigger exists in the menu
    //cout<<triggerNames_[itrig]<<" "<<triggerIndex_[itrig]<<" "<<accept<<" "<<tmpFired<<endl;
    Fired.push_back(tmpFired);
    L1Prescales.push_back(preL1);
    HLTPrescales.push_back(preHLT);
    mL1Objects.push_back(vvL1);
    mHLTObjects.push_back(vvHLT);
  }// loop over trigger names
  mEvent->setTrigDecision(Fired);
  mEvent->setPrescales(L1Prescales,HLTPrescales);
  mEvent->setL1Obj(mL1Objects);
  mEvent->setHLTObj(mHLTObjects);
  //-------------- Vertex Info -----------------------------------
  edm::Handle<reco::VertexCollection> recVtxs;
  event.getByLabel(mOfflineVertices,recVtxs);
  //------------- reject events without reco vertices ------------
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
  //-------------- Rho ------------------------------------------------
  edm::Handle<double> rhoCalo;
  event.getByLabel(mSrcCaloRho,rhoCalo);
  edm::Handle<double> rhoPF;
  event.getByLabel(mSrcPFRho,rhoPF);
  mEvtHdr.setRho(*rhoCalo,*rhoPF);
  //-------------- Generator Info -------------------------------------
  edm::Handle<GenEventInfoProduct> hEventInfo;
  edm::Handle<std::vector<reco::GenParticle> > genParticles;
  //-------------- Simulated PU Info ----------------------------------
  edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
  if (mIsMCarlo && mUseGenInfo) {
    event.getByLabel("generator", hEventInfo);
    mEvtHdr.setPthat(hEventInfo->binningValues()[0]);
    mEvtHdr.setWeight(hEventInfo->weight());
    event.getByLabel(mSrcPU, PupInfo);
    event.getByLabel("genParticles",genParticles);
    std::vector<PileupSummaryInfo>::const_iterator PUI;
    int nbx = PupInfo->size();
    int ootpuEarly(0),ootpuLate(0),intpu(0);
    float Tnpv = -1.; // new variable for computing pileup weight factor for the event
    for(PUI = PupInfo->begin(); PUI != PupInfo->end(); ++PUI) {
      if (PUI->getBunchCrossing() < 0)
        ootpuEarly += PUI->getPU_NumInteractions();
      else if (PUI->getBunchCrossing() > 0)
        ootpuLate += PUI->getPU_NumInteractions();
      else {
        intpu += PUI->getPU_NumInteractions();
        Tnpv = PUI->getTrueNumInteractions();
       }
    }

    float parton_id_initial1(0), parton_id_initial2(0);

    parton_id_initial1 = (*genParticles)[4].pdgId();
    parton_id_initial2 = (*genParticles)[5].pdgId();

    mEvtHdr.setPU(nbx,ootpuEarly,ootpuLate,intpu);
    mEvtHdr.setTrPu(Tnpv);
    mEvtHdr.setXsec(mXsec);
    mEvtHdr.setInitialPartons(parton_id_initial1, parton_id_initial2);
  }
  else {
    mEvtHdr.setPthat(0);
    mEvtHdr.setWeight(0);
    mEvtHdr.setPU(0,0,0,0);
    mEvtHdr.setTrPu(0);
    mEvtHdr.setXsec(0.);
    mEvtHdr.setInitialPartons(0., 0.);
  }
  //---------------- Jets ---------------------------------------------
  edm::ESHandle<JetCorrectorParametersCollection> PFJetCorParColl;
  if (mPFPayloadName != "" && !isPFJecUncSet_){
    iSetup.get<JetCorrectionsRecord>().get(mPFPayloadName,PFJetCorParColl);
    JetCorrectorParameters const& PFJetCorPar = (*PFJetCorParColl)["Uncertainty"];
    mPFUnc = new JetCorrectionUncertainty(PFJetCorPar);
    if (mPFJECUncSrc != "") {
      for(unsigned isrc=0;isrc<mPFJECUncSrcNames.size();isrc++) {
        JetCorrectorParameters *par = new JetCorrectorParameters(mPFJECUncSrc,mPFJECUncSrcNames[isrc]);
        JetCorrectionUncertainty *tmpUnc = new JetCorrectionUncertainty(*par);
        mPFUncSrc.push_back(tmpUnc);
      }
    }
    isPFJecUncSet_ = true;
  }
  edm::ESHandle<JetCorrectorParametersCollection> CaloJetCorParColl;
  if (mCaloPayloadName != "" && !isCaloJecUncSet_){
    iSetup.get<JetCorrectionsRecord>().get(mCaloPayloadName,CaloJetCorParColl);
    JetCorrectorParameters const& CaloJetCorPar = (*CaloJetCorParColl)["Uncertainty"];
    mCALOUnc = new JetCorrectionUncertainty(CaloJetCorPar);
    isCaloJecUncSet_ = true;
  }
  edm::Handle<GenJetCollection>  genjets;
  edm::Handle<std::vector<pat::Jet> > pfjets;
  edm::Handle<std::vector<pat::Jet> > calojets;
  edm::Handle<JetExtendedAssociation::Container> calojetExtender;

  event.getByLabel(mPFJetsName,pfjets);
  event.getByLabel(mCaloJetsName,calojets);
  event.getByLabel(mCaloJetExtender,calojetExtender);
  if (mIsMCarlo) {
    event.getByLabel(mGenJetsName,genjets);
    event.getByLabel("genParticles",genParticles);
    for(GenJetCollection::const_iterator i_gen = genjets->begin(); i_gen != genjets->end(); i_gen++) {
      if (i_gen->pt() > mMinGenPt && fabs(i_gen->y()) < mMaxY) {
        mGenJets.push_back(i_gen->p4());
      }
    }
  }
  int njets(0);
  //----------- PFJets -------------------------
  for(std::vector<pat::Jet>::const_iterator i_pfjet = pfjets->begin(); i_pfjet != pfjets->end(); i_pfjet++) {
    QCDPFJet qcdpfjet;

    double scale = i_pfjet->pt()/i_pfjet->correctedJet("Uncorrected").pt();
    //---- preselection -----------------
    if (fabs(i_pfjet->y()) > mMaxY) continue;
    //---- vertex association -----------
    //---- get the vector of tracks -----
    reco::PFJet const * pfJet = dynamic_cast<reco::PFJet const *>( i_pfjet->originalObject() );
    reco::TrackRefVector vTrks(pfJet->getTrackRefs());
    float sumTrkPt(0.0),sumTrkPtBeta(0.0),sumTrkPtBetaStar(0.0),beta(0.0),betaStar(0.0);
    //---- loop over the tracks of the jet ----
    for(reco::TrackRefVector::const_iterator i_trk = vTrks.begin(); i_trk != vTrks.end(); i_trk++) {
      if (recVtxs->size() == 0) break;
      sumTrkPt += (*i_trk)->pt();
      //---- loop over all vertices ----------------------------
      for(unsigned ivtx = 0;ivtx < recVtxs->size();ivtx++) {
        //---- loop over the tracks associated with the vertex ---
        if (!((*recVtxs)[ivtx].isFake()) && (*recVtxs)[ivtx].ndof() >= mGoodVtxNdof && fabs((*recVtxs)[ivtx].z()) <= mGoodVtxZ) {
          for(reco::Vertex::trackRef_iterator i_vtxTrk = (*recVtxs)[ivtx].tracks_begin(); i_vtxTrk != (*recVtxs)[ivtx].tracks_end(); ++i_vtxTrk) {
            //---- match the jet track to the track from the vertex ----
            reco::TrackRef trkRef(i_vtxTrk->castTo<reco::TrackRef>());
            //---- check if the tracks match -------------------------
            if (trkRef == (*i_trk)) {
              if (ivtx == 0) {
                sumTrkPtBeta += (*i_trk)->pt();
              }
              else {
                sumTrkPtBetaStar += (*i_trk)->pt();
              }
              break;
            }
          }
        }
      }
    }
    if (sumTrkPt > 0) {
      beta     = sumTrkPtBeta/sumTrkPt;
      betaStar = sumTrkPtBetaStar/sumTrkPt;
    }
    qcdpfjet.setBeta(beta);
    qcdpfjet.setBetaStar(betaStar);
    //---- jec uncertainty --------------
    double unc(0.0);
    vector<float> uncSrc(0);
    if (mPFPayloadName != "") {
      mPFUnc->setJetEta(i_pfjet->eta());
      mPFUnc->setJetPt(i_pfjet->pt());
      unc = mPFUnc->getUncertainty(true);
    }
    if (mPFJECUncSrc != "") {
      for(unsigned isrc=0;isrc<mPFJECUncSrcNames.size();isrc++) {
        mPFUncSrc[isrc]->setJetEta(i_pfjet->eta());
        mPFUncSrc[isrc]->setJetPt(i_pfjet->pt());
        float unc1 = mPFUncSrc[isrc]->getUncertainty(true);
        uncSrc.push_back(unc1);
      }
    }
    qcdpfjet.setP4(i_pfjet->correctedJet("Uncorrected").p4());
    qcdpfjet.setCor(scale);
    qcdpfjet.setUnc(unc);
    qcdpfjet.setUncSrc(uncSrc);
    qcdpfjet.setArea(i_pfjet->jetArea());
    double chf   = i_pfjet->chargedHadronEnergyFraction();
    double nhf   = (i_pfjet->neutralHadronEnergy() + i_pfjet->HFHadronEnergy())/i_pfjet->energy();
    double phf   = i_pfjet->photonEnergyFraction();
    double elf   = i_pfjet->electronEnergyFraction();
    double muf   = i_pfjet->muonEnergyFraction();
    double hf_hf = i_pfjet->HFHadronEnergyFraction();
    double hf_phf= i_pfjet->HFEMEnergyFraction();
    int hf_hm    = i_pfjet->HFHadronMultiplicity();
    int hf_phm   = i_pfjet->HFEMMultiplicity();
    int chm      = i_pfjet->chargedHadronMultiplicity();
    int nhm      = i_pfjet->neutralHadronMultiplicity();
    int phm      = i_pfjet->photonMultiplicity();
    int elm      = i_pfjet->electronMultiplicity();
    int mum      = i_pfjet->muonMultiplicity();
    int npr      = i_pfjet->chargedMultiplicity() + i_pfjet->neutralMultiplicity();
    retpf.set(false);
    bool looseID  = pfJetIDLoose( *i_pfjet, retpf );
    retpf.set(false);
    bool tightID  = pfJetIDTight( *i_pfjet, retpf );
    qcdpfjet.setLooseID(looseID);
    qcdpfjet.setTightID(tightID);
    qcdpfjet.setFrac(chf,nhf,phf,elf,muf);
    qcdpfjet.setMulti(npr,chm,nhm,phm,elm,mum);
    qcdpfjet.setHFFrac(hf_hf,hf_phf);
    qcdpfjet.setHFMulti(hf_hm,hf_phm);
    if (mIsMCarlo) {
      //
      float bquark_3 = 0.0;
      float bquark_2 = 0.0;
      float parton_id = 0.0;
      for (reco::GenParticleCollection::const_iterator igen_par = genParticles->begin(); igen_par != genParticles->end(); igen_par++) {
        double deltaR2 = reco::deltaR(*i_pfjet,*igen_par);
        int pdgid = igen_par->pdgId();
        int status = igen_par->status();
        if(deltaR2 < 0.35 && status == 3 && abs(pdgid) == 5) bquark_3 = 1.0;
        if(deltaR2 < 0.35 && status == 2 && abs(pdgid) == 5) bquark_2 = 1.0;
      }

       if (njets>2) parton_id = 0.;
       else if (njets<=2 && fabs((*genParticles)[6].pdgId()) >= 32) parton_id = (*genParticles)[7+njets].pdgId();
       else parton_id = (*genParticles)[6+njets].pdgId();

//        cout<<"njets: "<<njets<<"\t parton: "<<parton_id<<endl;
      njets = njets + 1;
     //
      if (i_pfjet->genJet() == 0) {
        LorentzVector tmpP4(0.0,0.0,0.0,0.0);
        qcdpfjet.setGen(tmpP4,0);
      }
      else
        qcdpfjet.setGen(i_pfjet->genJet()->p4(),reco::deltaR(*i_pfjet,*(i_pfjet->genJet())));

      qcdpfjet.setFlavor(i_pfjet->partonFlavour());
      qcdpfjet.setBstatus(bquark_3, bquark_2);
      qcdpfjet.setPartonId(parton_id);
    }
    else {
      LorentzVector tmpP4(0.0,0.0,0.0,0.0);
      qcdpfjet.setGen(tmpP4,0);
      qcdpfjet.setFlavor(-99.0);
      qcdpfjet.setBstatus(-99.0, -99.0);
      qcdpfjet.setPartonId(0.);
    }
    //
    qcdpfjet.setBtag_tche( i_pfjet->bDiscriminator("trackCountingHighEffBJetTags") );
    qcdpfjet.setBtag_tchp( i_pfjet->bDiscriminator("trackCountingHighPurBJetTags") );
    qcdpfjet.setBtag_csv( i_pfjet->bDiscriminator("combinedSecondaryVertexBJetTags") );
    qcdpfjet.setBtag_ssvhe( i_pfjet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags") );
    qcdpfjet.setBtag_ssvhp( i_pfjet->bDiscriminator("simpleSecondaryVertexHighPurBJetTags") );
    qcdpfjet.setBtag_jp( i_pfjet->bDiscriminator("jetProbabilityBJetTags") );
    //
    if (qcdpfjet.ptCor() >= mMinPFPt)
      mPFJets.push_back(qcdpfjet);
    if (qcdpfjet.ptCor() >= mMinPFFatPt && fabs(qcdpfjet.eta()) < mMaxPFFatEta && qcdpfjet.looseID())
      tmpPFJets.push_back(qcdpfjet);
  }
  //----------- PFFatJets ----------------------
  sort(tmpPFJets.begin(),tmpPFJets.end(),sort_pfjets);
  if (tmpPFJets.size()>1) {
    LorentzVector lead[2], fat[2];
    float sumPt[2],sumPtUnc[2];
    for(unsigned i = 0; i<2; i++) {
      lead[i]     = tmpPFJets[i].p4()*tmpPFJets[i].cor();
      fat[i]      = tmpPFJets[i].p4()*tmpPFJets[i].cor();
      sumPt[i]    = tmpPFJets[i].ptCor();
      sumPtUnc[i] = tmpPFJets[i].ptCor() * tmpPFJets[i].unc();
    }
    double rmax = 1.1;
    for(unsigned i = 2; i<tmpPFJets.size(); i++) {
      LorentzVector cand = tmpPFJets[i].p4();
      double dR1 = deltaR(lead[0],cand);
      double dR2 = deltaR(lead[1],cand);
      int index(-1);
      if (dR1 < dR2 && dR1 < rmax)
        index = 0;
      if (dR1 > dR2 && dR2 < rmax)
        index = 1;
      if (index > -1) {
        fat[index]      += cand * tmpPFJets[i].cor();
        sumPt[index]    += tmpPFJets[i].ptCor();
        sumPtUnc[index] += tmpPFJets[i].ptCor()*tmpPFJets[i].unc();
      }
    }
    QCDJet fatJet[2];
    vector<float> uncSrc(0);
    for(unsigned i = 0; i<2; i++) {
      fatJet[i].setP4(fat[i]);
      fatJet[i].setLooseID(tmpPFJets[i].looseID());
      fatJet[i].setTightID(tmpPFJets[i].tightID());
      fatJet[i].setCor(1.0);
      fatJet[i].setArea(0.0);
      fatJet[i].setUncSrc(uncSrc);
      //
      fatJet[i].setBtag_tche(tmpPFJets[i].btag_tche());
      fatJet[i].setBtag_tchp(tmpPFJets[i].btag_tchp());
      fatJet[i].setBtag_csv(tmpPFJets[i].btag_csv());
      fatJet[i].setBtag_ssvhe(tmpPFJets[i].btag_ssvhe());
      fatJet[i].setBtag_ssvhp(tmpPFJets[i].btag_ssvhp());
      fatJet[i].setBtag_jp(tmpPFJets[i].btag_jp());
      fatJet[i].setFlavor(tmpPFJets[i].flavor());
      fatJet[i].setBstatus(tmpPFJets[i].bstatus3(), tmpPFJets[i].bstatus2());
      fatJet[i].setPartonId(tmpPFJets[i].PartonId());
      //
      if (sumPt[i] > 0)
        fatJet[i].setUnc(sumPtUnc[i]/sumPt[i]);
      else
        fatJet[i].setUnc(0.0);
      fatJet[i].setGen(tmpPFJets[i].genp4(),tmpPFJets[i].genR());
    }
    if (fatJet[0].pt()>fatJet[1].pt()) {
      mPFFatJets.push_back(fatJet[0]);
      mPFFatJets.push_back(fatJet[1]);
    }
    else {
      mPFFatJets.push_back(fatJet[1]);
      mPFFatJets.push_back(fatJet[0]);
    }
  }
  //----------- CaloJets -----------------------
  for(std::vector<pat::Jet>::const_iterator i_calojet = calojets->begin(); i_calojet != calojets->end(); i_calojet++) {
    reco::Jet const * caloJet = dynamic_cast<reco::Jet const *>( i_calojet->originalObject() );

    double scale = i_calojet->pt()/i_calojet->correctedJet("Uncorrected").pt();
    //---- preselection -----------------
    if (fabs(i_calojet->y()) > mMaxY) continue;
    double unc(0.0);
    vector<float> uncSrc(0);
    if (mCaloPayloadName != "") {
      mCALOUnc->setJetEta(i_calojet->eta());
      mCALOUnc->setJetPt(i_calojet->pt());
      unc = mCALOUnc->getUncertainty(true);
    }
    QCDCaloJet qcdcalojet;
    qcdcalojet.setP4(i_calojet->correctedJet("Uncorrected").p4());
    qcdcalojet.setCor(scale);
    qcdcalojet.setUnc(unc);
    qcdcalojet.setUncSrc(uncSrc);
    qcdcalojet.setArea(i_calojet->jetArea());
    double emf    = i_calojet->emEnergyFraction();
    int n90hits   = i_calojet->jetID().n90Hits;
    double fHPD   = i_calojet->jetID().fHPD;
    double fRBX   = i_calojet->jetID().fRBX;
    int nTrkVtx   = JetExtendedAssociation::tracksAtVertexNumber(*calojetExtender, *caloJet);
    int nTrkCalo  = JetExtendedAssociation::tracksAtCaloNumber(*calojetExtender, *caloJet);
    ret.set(false);
    bool looseID  = jetIDLoose( *i_calojet, ret );
    ret.set(false);
    bool tightID  = jetIDTight( *i_calojet, ret );
    qcdcalojet.setVar(emf,fHPD,fRBX,n90hits,nTrkCalo,nTrkVtx);
    qcdcalojet.setLooseID(looseID);
    qcdcalojet.setTightID(tightID);
    if (mIsMCarlo) {
      float bquark_3 = 0.0;
      float bquark_2 = 0.0;
      for (reco::GenParticleCollection::const_iterator igen_par = genParticles->begin(); igen_par != genParticles->end(); igen_par++) {
        double deltaR2 = reco::deltaR(*i_calojet,*igen_par);
        int pdgid = igen_par->pdgId();
        int status = igen_par->status();
        if(deltaR2 < 0.35 && status == 3 && abs(pdgid) == 5) bquark_3 = 1.0;
        if(deltaR2 < 0.35 && status == 2 && abs(pdgid) == 5) bquark_2 = 1.0;
      }

      if (i_calojet->genJet() == 0) {
        LorentzVector tmpP4(0.0,0.0,0.0,0.0);
        qcdcalojet.setGen(tmpP4,0);
      }
      else
        qcdcalojet.setGen(i_calojet->genJet()->p4(),reco::deltaR(*i_calojet,*(i_calojet->genJet())));

      qcdcalojet.setFlavor(i_calojet->partonFlavour());
      qcdcalojet.setBstatus(bquark_3, bquark_2);
    }
    else {
      LorentzVector tmpP4(0.0,0.0,0.0,0.0);
      qcdcalojet.setGen(tmpP4,0);
      qcdcalojet.setFlavor(-99.0);
      qcdcalojet.setBstatus(-99, -99);
    }

    //
    qcdcalojet.setBtag_tche( i_calojet->bDiscriminator("trackCountingHighEffBJetTags") );
    qcdcalojet.setBtag_tchp( i_calojet->bDiscriminator("trackCountingHighPurBJetTags") );
    qcdcalojet.setBtag_csv( i_calojet->bDiscriminator("combinedSecondaryVertexBJetTags") );
    qcdcalojet.setBtag_ssvhe( i_calojet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags") );
    qcdcalojet.setBtag_ssvhp( i_calojet->bDiscriminator("simpleSecondaryVertexHighPurBJetTags") );
    qcdcalojet.setBtag_jp( i_calojet->bDiscriminator("jetProbabilityBJetTags") );
    //

    if (qcdcalojet.ptCor() >= mMinCaloPt)
      mCaloJets.push_back(qcdcalojet);
  }

  //---------------- met ---------------------------------------------
  edm::Handle<PFMETCollection> pfmet;
  edm::Handle<CaloMETCollection> calomet;
  event.getByLabel(mPFMETName,pfmet);
  event.getByLabel(mCaloMETName,calomet);
  mPFMet.setVar((*pfmet)[0].et(),(*pfmet)[0].sumEt(),(*pfmet)[0].phi());
  mCaloMet.setVar((*calomet)[0].et(),(*calomet)[0].sumEt(),(*calomet)[0].phi());

  //-------------- fill the tree -------------------------------------
  sort(mCaloJets.begin(),mCaloJets.end(),sort_calojets);
  sort(mPFJets.begin(),mPFJets.end(),sort_pfjets);
  mEvent->setEvtHdr(mEvtHdr);
  mEvent->setCaloJets(mCaloJets);
  mEvent->setPFJets(mPFJets);
  mEvent->setFatJets(mPFFatJets);
  mEvent->setGenJets(mGenJets);
  mEvent->setCaloMET(mCaloMet);
  mEvent->setPFMET(mPFMet);
  mEvent->setL1Obj(mL1Objects);
  mEvent->setHLTObj(mHLTObjects);
  if ((mEvent->nPFJets() >= (unsigned)mMinNPFJets) && (mEvent->nCaloJets() >= (unsigned)mMinNCaloJets)) {
    if ((mEvent->pfmjjcor(0) >= mMinJJMass) || (mEvent->calomjjcor(0) >= mMinJJMass) || (mEvent->fatmjjcor(0) >= mMinJJMass)) {
      mTree->Fill();
    }
  }
  //if (mPFPayloadName != "") {
    //delete mPFUnc;
    //delete mPFUncSrc;
 //}
  //if (mCaloPayloadName != "")
    //delete mCALOUnc;
}
//////////////////////////////////////////////////////////////////////////////////////////
ProcessedTreeProducer::~ProcessedTreeProducer()
{
}

DEFINE_FWK_MODULE(ProcessedTreeProducer);
