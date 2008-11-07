// -*- C++ -*-
//
// Package:    BitPlotting
// Class:      BitPlotting
// 
/**\class BitPlotting BitPlotting.cc UserCode/BitPlotting/src/BitPlotting.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jean-Roch Vlimant
//         Created:  Thu Jul 17 12:57:41 CEST 2008
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//
// class decleration
//

#include "FWCore/Framework/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Run.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h" 
#include "TH1F.h"
#include "TProfile.h"
#include "TH2F.h"

#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"


class BitPlotting : public edm::EDAnalyzer {
   public:
      explicit BitPlotting(const edm::ParameterSet&);
      ~BitPlotting();

private:
  virtual void beginJob(const edm::EventSetup&) {}
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
  virtual void beginRun(const edm::Run  & r, const edm::EventSetup  &);
  //  virtual void endRun(const edm::Run &, const edm::EventSetup &);

  edm::InputTag inputTag_;
  edm::TriggerNames triggerNames_;

  std::vector<std::string > HLTPathsByName_;
  std::vector<unsigned int> HLTPathsByIndex_;

  std::vector<unsigned int> count_;
  //  std::map<std::string, unsigned int> count_;

  unsigned int total_;

  std::string out;

  MonitorElement * h1_;
  MonitorElement * h2_;
  MonitorElement * pf_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
BitPlotting::BitPlotting(const edm::ParameterSet& iConfig) :
  inputTag_ (iConfig.getParameter<edm::InputTag> ("TriggerResultsTag")),
  HLTPathsByName_(iConfig.getParameter<std::vector<std::string > >("HLTPaths")),
  total_(0)
{
  count_.resize(HLTPathsByName_.size());
  HLTPathsByIndex_.resize(HLTPathsByName_.size());
  out = iConfig.getUntrackedParameter<std::string>("out");
}


BitPlotting::~BitPlotting(){}

//
// member functions
//

void BitPlotting::beginRun(const edm::Run  & r, const edm::EventSetup  &){
  
  std::stringstream rNs;
  rNs<<r.run();
  std::string rN = rNs.str();
  
  float min = -0.5;
  float max = HLTPathsByName_.size()-0.5;
  uint nBin = HLTPathsByName_.size();
  
  edm::LogError("BitPlotting")<<"this is the beginning of a NEW run: "<< r.run();

  //  h1_ = edm::Service<TFileService>()->make<TH1F>(("passingBits_"+rN).c_str(),("passing bits in run "+rN).c_str(),nBin,min,max);
  //  h2_ = edm::Service<TFileService>()->make<TH2F>(("passingBits_2_"+rN).c_str(),("correlation between bits in run "+rN).c_str(),nBin,min,max, nBin,min,max);
  //  pf_ = edm::Service<TFileService>()->make<TProfile>(("fraction_"+rN).c_str(),("fraction of passing bits in run "+rN).c_str(),nBin,min,max, 0, 1);

  edm::Service<DQMStore>()->setCurrentFolder("BitPlotting");

  h1_ = edm::Service<DQMStore>()->book1D(std::string("passingBits_"+rN),std::string("passing bits in run "+rN),nBin,min,max);
  h2_ = edm::Service<DQMStore>()->book2D(std::string("passingBits_2_"+rN),std::string("correlation between bits in run "+rN),nBin,min,max, nBin,min,max);
  pf_ = edm::Service<DQMStore>()->bookProfile(std::string("fraction_"+rN),std::string("fraction of passing bits in run "+rN),nBin,min,max, 1000, 0.0, 1.0);
  
  for (uint i=0; i!=nBin; ++i){
    h1_->getTH1F()->GetXaxis()->SetBinLabel(i+1,HLTPathsByName_[i].c_str());
    h2_->getTH2F()->GetXaxis()->SetBinLabel(i+1,HLTPathsByName_[i].c_str());
    h2_->getTH2F()->GetYaxis()->SetBinLabel(i+1,HLTPathsByName_[i].c_str());
    pf_->getTProfile()->GetXaxis()->SetBinLabel(i+1,HLTPathsByName_[i].c_str());
  }

}


// ------------ method called to for each event  ------------
void
BitPlotting::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   using namespace edm;
   using namespace std;

   const string invalid("@@invalid@@");

   // get hold of TriggerResults Object
   Handle<TriggerResults> trh;
   iEvent.getByLabel(inputTag_,trh);

   if (trh.failedToGet()) {
     edm::LogError("BitPlotting")<<" could not get: "<<inputTag_;
	 return;
   }

   // get hold of trigger names - based on TriggerResults object!
   triggerNames_.init(*trh);

   unsigned int n = HLTPathsByName_.size();
   for (unsigned int i=0; i!=n; i++) {
     HLTPathsByIndex_[i]=triggerNames_.triggerIndex(HLTPathsByName_[i]);
   }
   

   std::vector<bool> validity;
   validity.resize(n);
    
   for (unsigned int i=0; i!=n; i++) {
     validity[i]=( (HLTPathsByIndex_[i]<trh->size()) && (HLTPathsByName_[i]!=invalid) );
     
     /*
     LogTrace("") << " " << HLTPathsByIndex_[i]
		  << " " << HLTPathsByName_[i]
		  << " " << validity;

     if(!validity){
       edm::LogError("BitPlotting")<<HLTPathsByName_[i]<<" is not in the HLT menu.";
       return ;
     }
     */
     //     if (!validity) throw cms::Exception("Configuration");
       //       << " HLTHighLevel [instance: " << *moduleLabel()
       //       << " - path: " << *pathName()
       //       << "] configured with unknown HLT path name "
       //       << i << " " << HLTPathsByName_[i] <<"\n";
   }

   std::stringstream report;
   std::string sep=" ";
   bool atLeasOne=false;
   for (unsigned int i=0; i!=n; i++) {
     if (!validity[i]) continue;
     if (HLTPathsByIndex_[i]<trh->size()) {
       if (trh->accept(HLTPathsByIndex_[i])) {
	 report<<sep<<HLTPathsByName_[i];
	 count_[i]++;
	 sep=", ";
	 atLeasOne=true;
	 //	 edm::LogError("BitPlotting")<<"filling: "<<i<<" to: "<< h1_->getName()<<std::endl;
	 h1_->Fill(i);
	 pf_->Fill(i,1);
	 
	 for (unsigned int j=0; j!=n; j++) {
	   if (!validity[j]) continue;
	   if (HLTPathsByIndex_[j]<trh->size()) {
	     if (trh->accept(HLTPathsByIndex_[j])) {
	       h2_->Fill(i,j);
	     }
	   }
	 }
       }
       else{
	 pf_->Fill(i,0);
       }
     }
   }
   if (atLeasOne){
     edm::LogError("BitPlotting")<<report.str();
   }

   total_++;

   //   edm::LogError("BitPlotting")<<"# entries:"<<h1_->getTH1F()->GetEntries();

}


/*
void BitPlotting::endRun(const edm::Run &, const edm::EventSetup &){

  std::stringstream report;
  report<<" out of: "<<total_<<" events.\n";
  for (uint i=0; i!=HLTPathsByName_.size();i++){
    report<<HLTPathsByName_[i]<<" passed: "<<count_[i]<<" times.\n";
    count_[i]=0;
  }

  edm::LogError("BitPlotting")<<report.str();
  total_=0;
}
*/
// ------------ method called once each job just after ending the event loop  ------------
void 
BitPlotting::endJob() {

  std::stringstream report;
  report<<" out of: "<<total_<<" events.\n";
  for (uint i=0; i!=HLTPathsByName_.size();i++){
    report<<HLTPathsByName_[i]<<" passed: "<<count_[i]<<" times.\n";
    count_[i]=0;
  }

  edm::LogError("BitPlotting")<<report.str();
  total_=0;
  if(out) edm::Service<DQMStore>()->save(out);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BitPlotting);
