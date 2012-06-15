// -*- C++ -*-
//
// Package:    GluonTagLikelihood
// Class:      GluonTagLikelihood
// 

//
// Original Author:  Francesco FIORI
//         Created:  Wed May 16 11:03:08 CEST 2012
// $Id: GluonTagLikelihood.cc,v 1.1 2012/06/11 14:04:34 paoloa Exp $
//

#include <memory>
#include <TROOT.h>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "KKousour/MultiJetAnalysis/plugins/Classifier_Class.C"

using namespace std;
using namespace reco;
using namespace edm;
using namespace pat;

class GluonTagLikelihood : public edm::EDProducer {
   public:
      explicit GluonTagLikelihood(const edm::ParameterSet&);
      ~GluonTagLikelihood();

   private:
      virtual void beginJob();
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob();
      virtual double CorrectAndCompute(std::vector<double>&, double, double, double);
      virtual double Interpolate(double, int, int, double, double);

      // ----------member data ---------------------------
      edm::InputTag srcJets_,srcRho_;
      std::vector<std::string> vars_;  
      double corr_[6][3][5];
};

////////////////////////////////////////////////////////////////////////////////////////////////////
GluonTagLikelihood::GluonTagLikelihood(const edm::ParameterSet& iConfig)
{
  srcJets_ = iConfig.getParameter<edm::InputTag>("jets");
  srcRho_  = iConfig.getParameter<edm::InputTag>("rho");

  produces<edm::ValueMap<float> >().setBranchAlias("GluonTagLikelihood");
}
////////////////////////////////////////////////////////////////////////////////////////////////////
GluonTagLikelihood::~GluonTagLikelihood()
{
}
////////////////////////////////////////////////////////////////////////////////////////////////////
void GluonTagLikelihood::beginJob()
{
  const char* inVars[] = {"axis1","axis2","Mult","JetR","JetPull"};
  for (int i=0; i<5; i++) {
    vars_.push_back(inVars[i]);
  }
  double a[6][3][5] = { 
    {
      // a30c
      {0.000816, 0.000764, 0.000072, -0.012881, -0.002864},
      // a30t
      {0.002178, 0.001773, 0.000156, 0.285514, -0.005477},
      // a30f
      {0.001291, 0.001259, 0.000094, 0.203733, -0.004888},
    },
    {                             
      // a50c
      {0.000524, 0.000488, 0.000037, -0.005547, -0.001459},
      // a50t
      {0.001664, 0.001322, 0.000122, 0.256516, -0.004148},
      // a50f
      {0.001064, 0.000978, 0.000081, 0.202807, -0.004651},
    },
    {                                                   
      // a80c
      {0.000362, 0.000332, 0.000021, -0.002348, -0.000901},
      // a80t 
      {0.001140, 0.000923, 0.000075, 0.246119, -0.003209},
      // a80f[5]
      {0.000741, 0.000702, 0.000049, 0.202268, -0.003967},
    },
    {                                                             
      // a120c
      {0.000254, 0.000233, 0.000014, -0.001701, -0.000605},
      // a120t 
      {0.000777, 0.000646, 0.000044, 0.234861, -0.002532},
      // a120f 
      {0.000676, 0.000606, 0.000044, 0.212068, -0.003914},
    },
    {                                                   
      // a170c
      {0.000177, 0.000157, 0.000008, -0.005231, -0.000424},
      // a170t 
      {0.000555, 0.000457, 0.000026, 0.233618, -0.002029},
      // a170f 
      {0.000333, 0.000365, 0.000015, 0.208433, -0.002963},
    },
    {                                                     
      //a300 
      {0.000100, 0.000088, 0.000004, -0.006157, -0.000345},
      // a300t 
      {0.000387, 0.000298, 0.000015, 0.244468, -0.001433},
      // a300f 
      {-0.002678, -0.001875, -0.000106, 0.157607, 0.010459},
    }
  };
  for(int ipt=0;ipt<6;ipt++) {
    for(int ieta=0;ieta<3;ieta++) {
      for(int ipar=0;ipar<5;ipar++) {
        corr_[ipt][ieta][ipar] = a[ipt][ieta][ipar];
      }
    }
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////
void GluonTagLikelihood::endJob()
{
}
////////////////////////////////////////////////////////////////////////////////////////////////////
double GluonTagLikelihood::Interpolate(double pt, int ptlow, int pthigh, double mvalow, double mvahigh) 
{
  double interpolation;
  interpolation = (mvahigh-mvalow)/(pthigh-ptlow)*(pt-ptlow)+mvalow;
  return interpolation;
}
////////////////////////////////////////////////////////////////////////////////////////////////////
void GluonTagLikelihood::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::auto_ptr<std::vector<float> > mva(new std::vector<float>());

  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByLabel(srcJets_,jets);
  edm::View<pat::Jet> pat_jets = *jets;

  edm::Handle<double> rho;
  iEvent.getByLabel(srcRho_,rho);

  for(edm::View<pat::Jet>::const_iterator jet = pat_jets.begin(); jet != pat_jets.end(); ++jet) {
    vector<double> values;
    if (fabs(jet->eta()) < 2) {
      values.push_back(jet->userFloat("axis1_QC"));
      values.push_back(jet->userFloat("axis2_QC"));
      values.push_back(jet->userFloat("nChg_QC"));
      values.push_back(jet->userFloat("jetRchg_QC"));
      values.push_back(jet->userFloat("pull_QC"));
    }
    else {
      values.push_back(jet->userFloat("axis1"));
      values.push_back(jet->userFloat("axis2"));
      values.push_back(jet->userFloat("nChg_ptCut")+jet->userFloat("nNeutral_ptCut"));
      values.push_back(jet->userFloat("jetR"));
      values.push_back(jet->userFloat("pull"));
    }
    double mvaVal = CorrectAndCompute(values,jet->eta(),jet->pt(),*rho);
    mva->push_back(mvaVal);
  }
  std::auto_ptr<edm::ValueMap<float> > out(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler(*out);
  filler.insert(jets, mva->begin(), mva->end());
  filler.fill();
  // put value map into event
  iEvent.put(out);
}
////////////////////////////////////////////////////////////////////////////////////////////////////
double GluonTagLikelihood::CorrectAndCompute(std::vector<double>& vec, double eta, double pt, double rho) {

  if (vec.size() != 5) return -999;
  std::vector<double> vec1;
  vec1 = vec; //this is for the linear interpolation
  //---- align indices --------------
  int index[5] = {0,1,3,4,2};
  int ieta,ipt;
  if (fabs(eta) < 2) 
    ieta = 0;
  else if (fabs(eta) >= 2 && fabs(eta) < 3) 
    ieta = 1;
  else
    ieta = 2;
  if (pt >= 20 && pt < 30)
    ipt = 0;
  else if (pt >= 30 && pt < 50)
    ipt = 1;
  else if (pt >= 50 && pt < 80)
    ipt = 2;
  else if (pt >= 80 && pt < 120)
    ipt = 3;
  else if (pt >= 120 && pt < 170)
    ipt = 4; 
  else if (pt >= 170 && pt < 300)
    ipt = 5;
  else
    ipt = 6;

  for(int i=0;i<5;i++) {
    if (ieta == 0 && i == 2) continue; //in the central region the multiplicity is uncorrected
    if (ipt < 6) {
      vec[i] -= corr_[ipt][ieta][index[i]]*rho;
      vec1[i]-= corr_[ipt+1][ieta][index[i]]*rho;
    }
    else {
      vec[i] -= corr_[ipt-1][ieta][index[i]]*rho;
    }
  }
  double mvaval(-999);
  if (ieta == 0) {
    if (ipt == 0) {
      ReadLikelihood30C *MVAEvaluator = new ReadLikelihood30C(vars_);
      mvaval = MVAEvaluator->GetMvaValue(vec);
      delete MVAEvaluator;
    }
    else if (ipt == 1) {
      ReadLikelihood30C* MVAEvaluator = new ReadLikelihood30C(vars_);
      double mvaval1 = MVAEvaluator->GetMvaValue(vec);
      delete MVAEvaluator;

      ReadLikelihood50C* MVAEvaluator2 = new ReadLikelihood50C(vars_);
      double mvaval2 = MVAEvaluator2->GetMvaValue(vec1);
      delete MVAEvaluator2;

      mvaval = Interpolate(pt,30,50,mvaval1,mvaval2);
    }
    else if (ipt == 2) {
      ReadLikelihood50C* MVAEvaluator = new ReadLikelihood50C(vars_);
      double mvaval1 = MVAEvaluator->GetMvaValue(vec);
      delete MVAEvaluator;
    
      ReadLikelihood80C* MVAEvaluator2 = new ReadLikelihood80C(vars_);
      double mvaval2 = MVAEvaluator2->GetMvaValue(vec1);
      delete MVAEvaluator2;

      mvaval = Interpolate(pt,50,80,mvaval1,mvaval2);
    }
    else if (ipt == 3) {
      ReadLikelihood80C* MVAEvaluator = new ReadLikelihood80C(vars_);
      double mvaval1 = MVAEvaluator->GetMvaValue(vec);
      delete MVAEvaluator;

      ReadLikelihood120C* MVAEvaluator2 = new ReadLikelihood120C(vars_);
      double mvaval2 = MVAEvaluator2->GetMvaValue(vec1);
      delete MVAEvaluator2;

      mvaval = Interpolate(pt,80,120,mvaval1,mvaval2);
    }
    else if (ipt == 4) {
      ReadLikelihood120C* MVAEvaluator = new ReadLikelihood120C(vars_);
      double mvaval1 = MVAEvaluator->GetMvaValue(vec);
      delete MVAEvaluator;

      ReadLikelihood170C* MVAEvaluator2 = new ReadLikelihood170C(vars_);
      double mvaval2 = MVAEvaluator2->GetMvaValue(vec1);
      delete MVAEvaluator2;

      mvaval = Interpolate(pt,120,170,mvaval1,mvaval2);
    }
    else if (ipt == 5) {
      ReadLikelihood170C* MVAEvaluator = new ReadLikelihood170C(vars_);
      double mvaval1 = MVAEvaluator->GetMvaValue(vec);
      delete MVAEvaluator;

      ReadLikelihood300C* MVAEvaluator2 = new ReadLikelihood300C(vars_);
      double mvaval2 = MVAEvaluator2->GetMvaValue(vec1);
      delete MVAEvaluator2;

      mvaval = Interpolate(pt,170,300,mvaval1,mvaval2);
    }
    else {
      ReadLikelihood300C* MVAEvaluator = new ReadLikelihood300C(vars_);
      mvaval = MVAEvaluator->GetMvaValue(vec);
      delete MVAEvaluator;
    }
  }
  // transition
  if (ieta == 1) {
    if (ipt == 0) {
      ReadLikelihood30T* MVAEvaluator = new ReadLikelihood30T(vars_);
      mvaval = MVAEvaluator->GetMvaValue(vec);
      delete MVAEvaluator;
    }
    else if (ipt == 1) {
      ReadLikelihood30T* MVAEvaluator = new ReadLikelihood30T(vars_);
      double mvaval1 = MVAEvaluator->GetMvaValue(vec);
      delete MVAEvaluator;

      ReadLikelihood50T* MVAEvaluator2 = new ReadLikelihood50T(vars_);
      double mvaval2 = MVAEvaluator2->GetMvaValue(vec1);
      delete MVAEvaluator2;

      mvaval = Interpolate(pt,30,50,mvaval1,mvaval2);
    }
    else if (ipt == 2) {
      ReadLikelihood50T* MVAEvaluator = new ReadLikelihood50T(vars_);
      double mvaval1 = MVAEvaluator->GetMvaValue(vec);
      delete MVAEvaluator;

      ReadLikelihood80T* MVAEvaluator2 = new ReadLikelihood80T(vars_);
      double mvaval2 = MVAEvaluator2->GetMvaValue(vec1);
      delete MVAEvaluator2;

      mvaval = Interpolate(pt,50,80,mvaval1,mvaval2);
    }
    else if (ipt == 3) {
      ReadLikelihood80T* MVAEvaluator = new ReadLikelihood80T(vars_);
      double mvaval1=MVAEvaluator->GetMvaValue(vec);
      delete MVAEvaluator;
   
      ReadLikelihood120T* MVAEvaluator2 = new ReadLikelihood120T(vars_);
      double mvaval2 = MVAEvaluator2->GetMvaValue(vec1);
      delete MVAEvaluator2;
      mvaval = Interpolate(pt,80,120,mvaval1,mvaval2);
    }
    else if (ipt == 4) {
      ReadLikelihood120T* MVAEvaluator = new ReadLikelihood120T(vars_);
      double mvaval1 = MVAEvaluator->GetMvaValue(vec);
      delete MVAEvaluator;

      ReadLikelihood170T* MVAEvaluator2 = new ReadLikelihood170T(vars_);
      double mvaval2 = MVAEvaluator2->GetMvaValue(vec1);
      delete MVAEvaluator2;

      mvaval = Interpolate(pt,120,170,mvaval1,mvaval2);
    }
    else if (ipt == 5) {
      ReadLikelihood170T* MVAEvaluator = new ReadLikelihood170T(vars_);
      double mvaval1 = MVAEvaluator->GetMvaValue(vec);
      delete MVAEvaluator;

      ReadLikelihood300T* MVAEvaluator2 = new ReadLikelihood300T(vars_);
      double mvaval2 = MVAEvaluator2->GetMvaValue(vec1);
      delete MVAEvaluator2;

      mvaval = Interpolate(pt,170,300,mvaval1,mvaval2);
    }
    else {
      ReadLikelihood300T* MVAEvaluator = new ReadLikelihood300T(vars_);
      mvaval = MVAEvaluator->GetMvaValue(vec);
      delete MVAEvaluator;
    }
  }
  //forward region
  if (ieta == 2) {
    if (ipt == 0) {
      ReadLikelihood30F* MVAEvaluator = new ReadLikelihood30F(vars_);
      mvaval = MVAEvaluator->GetMvaValue(vec);
      delete MVAEvaluator;
    }
    else if (ipt == 1) {
      ReadLikelihood30F* MVAEvaluator = new ReadLikelihood30F(vars_);
      double mvaval1 = MVAEvaluator->GetMvaValue(vec);
      delete MVAEvaluator;

      ReadLikelihood50F* MVAEvaluator2 = new ReadLikelihood50F(vars_);
      double mvaval2 = MVAEvaluator2->GetMvaValue(vec1);
      delete MVAEvaluator2;

      mvaval = Interpolate(pt,30,50,mvaval1,mvaval2);
    }
    else if (ipt == 2) {
      ReadLikelihood50F* MVAEvaluator = new ReadLikelihood50F(vars_);
      double mvaval1 = MVAEvaluator->GetMvaValue(vec);
      delete MVAEvaluator;

      ReadLikelihood80F* MVAEvaluator2 = new ReadLikelihood80F(vars_);
      double mvaval2 = MVAEvaluator2->GetMvaValue(vec1);
      delete MVAEvaluator2;

      mvaval = Interpolate(pt,50,80,mvaval1,mvaval2);
    }
    else if (ipt == 3) {
      ReadLikelihood80F* MVAEvaluator = new ReadLikelihood80F(vars_);
      double mvaval1 = MVAEvaluator->GetMvaValue(vec);
      delete MVAEvaluator;

      ReadLikelihood120F* MVAEvaluator2 = new ReadLikelihood120F(vars_);
      double mvaval2 = MVAEvaluator2->GetMvaValue(vec1);
      delete MVAEvaluator2;

      mvaval = Interpolate(pt,80,120,mvaval1,mvaval2);
    }
    else if (ipt == 4){
      ReadLikelihood120F* MVAEvaluator = new ReadLikelihood120F(vars_);
      double mvaval1 = MVAEvaluator->GetMvaValue(vec);
      delete MVAEvaluator;

      ReadLikelihood170F* MVAEvaluator2 = new ReadLikelihood170F(vars_);
      double mvaval2 = MVAEvaluator2->GetMvaValue(vec1);
      delete MVAEvaluator2;

      mvaval=Interpolate(pt,120,170,mvaval1,mvaval2);
    }
    else if (ipt == 5){
      ReadLikelihood170F* MVAEvaluator = new ReadLikelihood170F(vars_);
      double mvaval1 = MVAEvaluator->GetMvaValue(vec);
      delete MVAEvaluator;

      ReadLikelihood300F* MVAEvaluator2 = new ReadLikelihood300F(vars_);
      double mvaval2 = MVAEvaluator2->GetMvaValue(vec1);
      delete MVAEvaluator2;

      mvaval = Interpolate(pt,170,300,mvaval1,mvaval2);
    }
    else {
      mvaval = -999; //kinematically forbidden
    }
  }
  return mvaval;
}

//define this as a plug-in
DEFINE_FWK_MODULE(GluonTagLikelihood);
