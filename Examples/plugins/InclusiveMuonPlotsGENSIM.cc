/** \class InclusiveMuonPlotsGENSIM
 *  Make inclusive muon plots
 *
 *  Based on r 1.4 of InclusiveMuonPlots by G.P.  This version 
 *  adds figures that only I am interested in . . . 
 *
 *  \author G. Petrucciani - UCSD (Giovanni.Petrucciani@cern.ch), ...
 */

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "MuonAnalysis/Examples/interface/muonStations.h"

// for "luminosity"
#include "DataFormats/Common/interface/MergeableCounter.h"

// for selection cut
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "PhysicsTools/HepMCCandAlgos/interface/GenParticlesHelper.h"

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TObjString.h>
#include <TDirectory.h>

#include <map>
#include <string>

#include "boost/lexical_cast.hpp"

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

using boost::lexical_cast;

class InclusiveMuonPlotsGENSIM: public edm::EDAnalyzer {
    public:
        /// Constructor
        InclusiveMuonPlotsGENSIM(const edm::ParameterSet& pset) ;

        /// Destructor
        virtual ~InclusiveMuonPlotsGENSIM() ;

        // Operations
        void analyze(const edm::Event & event, const edm::EventSetup& eventSetup) ;

        void endLuminosityBlock(const edm::LuminosityBlock &, const edm::EventSetup&);

        void book(const TFileService &fs, const edm::ParameterSet &pset, const std::string &name, const std::string &basename) ;
        void book2d(const TFileService &fs, const edm::ParameterSet &pset, const std::string &name, const std::string &basename) ;
        void book(const TFileService &fs, const edm::ParameterSet &pset, const std::string &name) { book(fs,pset,name,name); }

        void bookProf(const TFileService &fs, const edm::ParameterSet &pset, const std::string &name, const std::string &basename) ;
        void bookProf(const TFileService &fs, const edm::ParameterSet &pset, const std::string &name) { bookProf(fs,pset,name,name); }

    private:
        edm::InputTag muons_;
        edm::InputTag particleSrc_;
  //edm::InputTag dilepton_src_;
        StringCutObjectSelector<reco::GenParticle> selector_;
        StringCutObjectSelector<pat::Muon> selectorReco_;

        edm::InputTag primaryVertices_;
        edm::InputTag normalization_;

        // we don't care too much about performance
        std::map<std::string, TH1*>      plots;
        std::map<std::string, TProfile*> profiles;

        TH1D *luminosity;
};

/// Constructor
InclusiveMuonPlotsGENSIM::InclusiveMuonPlotsGENSIM(const edm::ParameterSet& pset):

    muons_(pset.getParameter<edm::InputTag>("muons")),
    particleSrc_(pset.getParameter<edm::InputTag>("particleSrc")),
    //dilepton_src_(pset.getParameter<edm::InputTag>("dilepton_src")),
    selector_(pset.getParameter<std::string>("selection")),
    selectorReco_(pset.getParameter<std::string>("selectionReco")),
    primaryVertices_(pset.getParameter<edm::InputTag>("primaryVertices")),
    luminosity(0) // by default, we don't have luminosity info
{

    edm::Service<TFileService> fs;

    TFileDirectory md = fs->mkdir("metadata");
    TDirectory *md_dir = md.cd();
    md_dir->WriteTObject(new TObjString(muons_.encode().c_str()), "muons");
    md_dir->WriteTObject(new TObjString(pset.getParameter<std::string>("selection").c_str()), "selection");
    md_dir->WriteTObject(new TObjString(pset.getParameter<std::string>("selectionReco").c_str()), "selectionReco");

    book(*fs, pset, "p"); 
    book(*fs, pset, "pt"); 
    book(*fs, pset, "eta"); 
    book(*fs, pset, "phi"); 
    book(*fs, pset, "charge"); 
    book(*fs, pset, "pdg"); 
    book(*fs, pset, "mass");

    book(*fs, pset, "prodz", "z"); 
    book(*fs, pset, "prodr", "r"); 
    book(*fs, pset, "prodd", "r"); 
    book2d(*fs, pset, "prodrz", "rz"); 
    book2d(*fs, pset, "pteta", "pteta");
    book2d(*fs, pset, "peta", "peta");

    book(*fs, pset, "p_reco","p"); 
    book(*fs, pset, "pt_reco","pt"); 
    book(*fs, pset, "eta_reco","eta"); 
    book(*fs, pset, "phi_reco","phi"); 
    book2d(*fs, pset, "pteta_reco", "pteta");
    book2d(*fs, pset, "peta_reco", "peta");

    book(*fs, pset, "p_","p"); 
    book(*fs, pset, "pt_","pt"); 
    book(*fs, pset, "eta_","eta"); 
    book(*fs, pset, "phi_","phi"); 
    book2d(*fs, pset, "pteta_", "pteta");
    book2d(*fs, pset, "peta_", "peta");

    if (pset.existsAs<edm::InputTag>("normalization")) {
        normalization_ = pset.getParameter<edm::InputTag>("normalization");
        luminosity = fs->make<TH1D>("normalization", "normalization", 1, 0, 1);
        luminosity->Sumw2();
    }
}

/// Destructor
InclusiveMuonPlotsGENSIM::~InclusiveMuonPlotsGENSIM()
{

}

void InclusiveMuonPlotsGENSIM::book(const TFileService &fs, const edm::ParameterSet &pset, const std::string &name, const std::string &basename) 
{
    typedef std::vector<double> vdouble;
    if (pset.existsAs<vdouble>(basename+"Bins")) {
        vdouble bins = pset.getParameter<vdouble>(basename+"Bins");
        plots[name] = fs.make<TH1D>(name.c_str(), name.c_str(), bins.size()-1, &bins[0]);
    } else {
        uint32_t nbins = pset.getParameter<uint32_t>(basename+"Bins");
        vdouble  range = pset.getParameter<vdouble>(basename+"Range");
        if (range.size() != 2) throw cms::Exception("Configuration") << "parameter '" << basename << "Range' is not of the form (min, max).\n";
        plots[name] = fs.make<TH1D>(name.c_str(), name.c_str(), nbins, range[0], range[1]);
    }
}

void InclusiveMuonPlotsGENSIM::book2d(const TFileService &fs, const edm::ParameterSet &pset, const std::string &name, const std::string &basename) 
{
    typedef std::vector<double> vdouble;
    //if (pset.existsAs<vdouble>(basename+"Bins")) {
    //vdouble bins = pset.getParameter<vdouble>(basename+"Bins");
    //plots[name] = fs.make<TH1D>(name.c_str(), name.c_str(), bins.size()-1, &bins[0]);
    //} else {
    uint32_t nbinsx = pset.getParameter<uint32_t>(basename+"XBins");
    vdouble  rangex = pset.getParameter<vdouble>(basename+"XRange");
    uint32_t nbinsy = pset.getParameter<uint32_t>(basename+"YBins");
    vdouble  rangey = pset.getParameter<vdouble>(basename+"YRange");
    if (rangex.size() != 2 || rangey.size() != 2) throw cms::Exception("Configuration") << "parameter '" << basename <<"Range' is not of the form (min, max).\n";
    plots[name] = fs.make<TH2D>(name.c_str(), name.c_str(), nbinsx, rangex[0], rangex[1],nbinsy,rangey[0],rangey[1]);
    //}
}

void InclusiveMuonPlotsGENSIM::bookProf(const TFileService &fs, const edm::ParameterSet &pset, const std::string &name, const std::string &basename) 
{
    typedef std::vector<double> vdouble;
    if (pset.existsAs<vdouble>(basename+"Bins")) {
        vdouble bins = pset.getParameter<vdouble>(basename+"Bins");
        profiles[name] = fs.make<TProfile>(name.c_str(), name.c_str(), bins.size()-1, &bins[0]);
    } else {
        uint32_t nbins = pset.getParameter<uint32_t>(basename+"Bins");
        vdouble  range = pset.getParameter<vdouble>(basename+"Range");
        if (range.size() != 2) throw cms::Exception("Configuration") << "parameter '" << basename << "Range' is not of the form (min, max).\n";
        profiles[name] = fs.make<TProfile>(name.c_str(), name.c_str(), nbins, range[0], range[1]);
    }
}


void InclusiveMuonPlotsGENSIM::analyze(const edm::Event & event, const edm::EventSetup& eventSetup){
    using namespace edm;
    using namespace std;

    //std::cout << "Event " << event.id() << std::endl;

    // get generator particle collection
    Handle<reco::GenParticleCollection> particles;
    event.getByLabel(particleSrc_, particles);

    Handle<View<reco::Muon> > muons;
    event.getByLabel(muons_, muons);
    
    //Handle<reco::CompositeCandidateView> dileptons;
    //event.getByLabel(dilepton_src_, dileptons);

    Handle<vector<reco::Vertex> > vertices;
    event.getByLabel(primaryVertices_, vertices);
    int k =0;

    for(reco::GenParticleCollection::const_iterator part=particles->begin(); 
	part!=particles->end(); ++part){
      
      if (!selector_(*part)) continue;

      //if(part->mother().isValid()) plots["mass"]->Fill(part->mother()->mass());      
      plots["p"  ]->Fill(part->p());
      plots["pt" ]->Fill(part->pt());
      plots["eta"]->Fill(part->eta());
      plots["phi"]->Fill(part->phi());
      plots["charge"]->Fill(part->charge());
      plots["pdg"]->Fill(part->pdgId());
      
      plots["prodz"]->Fill(part->vz());
      float r = sqrt(part->vx()*part->vx() + part->vy()*part->vy());
      float d = sqrt(part->vx()*part->vx() + part->vy()*part->vy() + part->vz()*part->vz());
      plots["prodr"]->Fill(r);
      plots["prodd"]->Fill(d);
      plots["prodrz"]->Fill(part->vz(),r);
      plots["pteta"]->Fill(part->eta(),part->pt());
      plots["peta"]->Fill(part->eta(),part->p());
    }
    
    ////////////start
    reco::GenParticleRefVector allStatus3Zs; 
    using namespace GenParticlesHelper;
    findParticles( *particles,allStatus3Zs, 23, 3);
    //std::cout << " nZs " << allStatus3Zs.size() << std::endl;
    for( IGR iZ=allStatus3Zs.begin(); iZ!=allStatus3Zs.end(); ++iZ) {
      // look for all status 1 (stable) descendents 
      reco::GenParticleRefVector descendents;
      findDescendents( *iZ, descendents, 3,13);
      //std::cout << "  daughters " << descendents.size() << std::endl;
      
      //////
      for(IGR igr = descendents.begin(); 
	  igr!= descendents.end(); ++igr ) {
	plots["p_reco"  ]->Fill((*igr)->p());
	plots["pt_reco" ]->Fill((*igr)->pt());
	plots["eta_reco"]->Fill((*igr)->eta());
	plots["phi_reco"]->Fill((*igr)->phi());
	plots["pteta_reco"]->Fill((*igr)->eta(),(*igr)->pt());
	plots["peta_reco"]->Fill((*igr)->eta(),(*igr)->p());
      }
      //////
      if(descendents.size() >= 2) {
	if( ! ( (abs(descendents[0]->eta()) < 2.1 && abs(descendents[1]->eta()) < 2.1)
		&&  ( (descendents[0]->pt()>10. && descendents[1]->pt()>6.) || (descendents[0]->pt()>6. && descendents[1]->pt()>10.) ) ) ) continue;
	//std::cout << "     " <<  descendents[0]->pt() << " " << descendents[1]->pt() << std::endl;
	for(IGR igr = descendents.begin(); 
	    igr!= descendents.end(); ++igr ) {
	  plots["p_"  ]->Fill((*igr)->p());
	  plots["pt_" ]->Fill((*igr)->pt());
	  plots["eta_"]->Fill((*igr)->eta());
	  plots["phi_"]->Fill((*igr)->phi());
	  plots["pteta_"]->Fill((*igr)->eta(),(*igr)->pt());
	  plots["peta_"]->Fill((*igr)->eta(),(*igr)->p());
	}
      }
    }
    //////////////////end
    
    foreach (const reco::Muon &recomu, *muons) {
      // we want to make a pat::Muon so that we can access directly muonID in the cuts
      const pat::Muon &mu = (typeid(recomu) == typeid(pat::Muon) ? static_cast<const pat::Muon &>(recomu) : pat::Muon(recomu));
      
      if (! selectorReco_(mu)) continue;
      if (! mu.genParticleRef().isAvailable()) continue;
      //if (! selector_(*mu.genParticleRef())) continue;
      /*
      plots["p_"  ]->Fill(mu.genParticleRef()->p());
      plots["pt_" ]->Fill(mu.genParticleRef()->pt());
      plots["eta_"]->Fill(mu.genParticleRef()->eta());
      plots["phi_"]->Fill(mu.genParticleRef()->phi());
      plots["pteta_"]->Fill(mu.genParticleRef()->eta(),mu.genParticleRef()->pt());
      plots["peta_"]->Fill(mu.genParticleRef()->eta(),mu.genParticleRef()->p());	
      */
      /*
      plots["p_reco"  ]->Fill(mu.p());
      plots["pt_reco" ]->Fill(mu.pt());
      plots["eta_reco"]->Fill(mu.eta());
      plots["phi_reco"]->Fill(mu.phi());
      plots["pteta_reco"]->Fill(mu.eta(),mu.pt());
      plots["peta_reco"]->Fill(mu.eta(),mu.p());	
      */
      }
    

}

void InclusiveMuonPlotsGENSIM::endLuminosityBlock(const edm::LuminosityBlock & iLumi, const edm::EventSetup & iSetup) 
{
    if (luminosity != 0) {
        edm::Handle<edm::MergeableCounter> mc;
        iLumi.getByLabel(normalization_, mc);
        luminosity->Fill(0.5, double(mc->value));
        // set the correct uncertainty from counting statistics
        luminosity->SetBinError(1, sqrt(luminosity->GetBinContent(1)));
    }
}

DEFINE_FWK_MODULE(InclusiveMuonPlotsGENSIM);







