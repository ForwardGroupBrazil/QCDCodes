/** \class STAMuonAnalyzer
 *  Analyzer of the StandAlone muon tracks
 *
 *  $Date: 2006/09/01 14:35:48 $
 *  $Revision: 1.4 $
 *  \author R. Bellan - INFN Torino <riccardo.bellan@cern.ch>
 */

#include "UserCode/AEverett/test/STAMuonAnalyzer.h"

// Collaborating Class Header
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "RecoMuon/TrackingTools/interface/MuonPatternRecoDumper.h"

#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace std;
using namespace edm;

/// Constructor
STAMuonAnalyzer::STAMuonAnalyzer(const ParameterSet& pset){
  theSTAMuonLabel = pset.getUntrackedParameter<string>("StandAloneTrackCollectionLabel");
  theSeedCollectionLabel = pset.getUntrackedParameter<string>("MuonSeedCollectionLabel");

  theRootFileName = pset.getUntrackedParameter<string>("rootFileName");

  theDataType = pset.getUntrackedParameter<string>("DataType");
  
  if(theDataType != "RealData" && theDataType != "SimData")
    cout<<"Error in Data Type!!"<<endl;

  numberOfSimTracks=0;
  numberOfRecTracks=0;
}

/// Destructor
STAMuonAnalyzer::~STAMuonAnalyzer(){
}

void STAMuonAnalyzer::beginJob(const EventSetup& eventSetup){
  // Create the root file
  theFile = new TFile(theRootFileName.c_str(), "RECREATE");
  theFile->cd();

  hPtRec = new TH1F("pTRec","p_{T}^{rec}",250,0,120);
  hPtSim = new TH1F("pTSim","p_{T}^{gen} ",250,0,120);

  hPTDiff = new TH1F("pTDiff","p_{T}^{rec} - p_{T}^{gen} ",250,-120,120);
  hPTDiff2 = new TH1F("pTDiff2","p_{T}^{rec} - p_{T}^{gen} ",250,-120,120);

  hPTDiffvsEta = new TH2F("PTDiffvsEta","p_{T}^{rec} - p_{T}^{gen} VS #eta",100,-2.5,2.5,250,-120,120);
  hPTDiffvsPhi = new TH2F("PTDiffvsPhi","p_{T}^{rec} - p_{T}^{gen} VS #phi",100,-6,6,250,-120,120);

  hPres = new TH1F("pTRes","pT Resolution",100,-2,2);
  h1_Pres = new TH1F("invPTRes","1/pT Resolution",100,-2,2);

// my histos...
 hNsim = new TH1F("hNsim","Nmu sim",20,0,20);
 hNrec = new TH1F("hNrec","Nmu rec",20,0,20);
 hNmatch = new TH1F("hNmatch","Nmu rec matched",20,0,20);
 hetasim = new TH1F("hetasim","eta sim",50,-2.5,2.5);
 hptsim  = new TH1F("hptsim","pt rec",80,0.,160.);
 hdR  = new TH1F("hdR","dR ",50,0.,2.5);
 hetarec = new TH1F("hetarec","eta rec",50,-2.5,2.5);
 hetarec1 = new TH1F("hetarec1","eta rec",50,-2.5,2.5);
 hetarec_all = new TH1F("hetarec_all","eta rec (all)",50,-2.5,2.5);
 hptrec  = new TH1F("hptrec","pt rec",80,0.,160.);
 hptrec1  = new TH1F("hptrec1","pt rec",80,0.,160.);
 hdptrec  = new TH1F("hdptrec","pt rec",50,-1.,1.);
 hdptrecvsEta = new TH2F("hdptrecvsEta","(pTrec-pTgen)/pTgen vs eta",50,-2.5,2.5,80,-2.,2.);
 hptrecvsEta = new TH2F("hptrecvsEta","pTrec unmatched vs eta",50,-2.5,2.5,80,0.,200.);
 hminv  = new TH1F("hminv","Minv",60,0.,120.);
 hminvrec  = new TH1F("hminvrec","Minvrec",70,0.,140.);
 hminvrecZ = new TH1F("hminvrecZ","Minvrec",70,0.,140.);

 hetaeff = new TH1F("hetaeff","Efficiency",50,-2.5,2.5);
 hetapur = new TH1F("hetapur","Purity",50,-2.5,2.5);

}

void STAMuonAnalyzer::endJob(){
  if(theDataType == "SimData"){
    cout << endl << endl << "Number of Sim tracks: " << numberOfSimTracks << endl;
  }

  cout << "Number of Reco tracks: " << numberOfRecTracks << endl << endl;
    
  hetaeff->Divide(hetarec,hetasim,1,1);
  hetapur->Divide(hetarec1,hetarec_all,1,1);

  // Write the histos to file
  theFile->cd();
  hPtRec->Write();
  hPtSim->Write();
  hPres->Write();
  h1_Pres->Write();
  hPTDiff->Write();
  hPTDiff2->Write();
  hPTDiffvsEta->Write();
  hPTDiffvsPhi->Write();

  hNsim->Write();
  hNrec->Write();
  hNmatch->Write();
  hetasim->Write();
  hptsim->Write();
  hdR->Write();
  hetarec->Write();
  hetarec1->Write();
  hetarec_all->Write();
  hptrec->Write();
  hptrec1->Write();
  hdptrec->Write(); 
  hdptrecvsEta->Write();
  hptrecvsEta->Write(); 
  hminv->Write();  
  hminvrec->Write(); 
  hminvrecZ->Write();   

  hetaeff->Write();
  hetapur->Write();

  theFile->Close();
}
 

void STAMuonAnalyzer::analyze(const Event & event, const EventSetup& eventSetup){
  
  cout << "Run: " << event.id().run() << " Event: " << event.id().event() << endl;
  MuonPatternRecoDumper debug;
  
  // Get the RecTrack collection from the event
  Handle<reco::TrackCollection> staTracks;
  event.getByLabel(theSTAMuonLabel, staTracks);

  ESHandle<MagneticField> theMGField;
  eventSetup.get<IdealMagneticFieldRecord>().get(theMGField);

  ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
  eventSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);
  
  double recPt=0.;
  double simPt=0.;
  float pxsim[20], pysim[20],pzsim[20], etasim[20], phisim[20], ptsim[20], qsim[20];  // sim muons
  float pxrec[20], pyrec[20],pzrec[20], etarec[20], phirec[20], ptrec[20], qrec[20];  // rec muons
  int nmusim = 0;
  int nmurec = 0;
  
  // Get the SimTrack collection from the event
  if(theDataType == "SimData"){
    Handle<SimTrackContainer> simTracks;
    event.getByLabel("g4SimHits",simTracks);
    
    numberOfRecTracks += staTracks->size();

    SimTrackContainer::const_iterator simTrack;

// loop on simulated  tracks.....-----------------
    cout<<"Simulated tracks: "<<endl;
    for (simTrack = simTracks->begin(); simTrack != simTracks->end(); ++simTrack){
      if (abs((*simTrack).type()) == 13) {
	cout<<"Sim pT: "<<(*simTrack).momentum().perp()<<endl;
	simPt=(*simTrack).momentum().perp();
        float simEta=(*simTrack).momentum().eta();
	if (simPt > 4. && abs(simEta) < 2.4 )  {
 	  nmusim++;
	  pxsim[nmusim] = (*simTrack).momentum().x();
	  pysim[nmusim] = (*simTrack).momentum().y();
	  pzsim[nmusim] = (*simTrack).momentum().z();
          etasim[nmusim] = (*simTrack).momentum().eta();
	  phisim[nmusim] = (*simTrack).momentum().phi();
          ptsim[nmusim] = (*simTrack).momentum().perp();
          qsim[nmusim] = -1.;
	  if( (*simTrack).type() == -13 ) qsim[nmusim] = 1.;
	} 
	cout<<"Sim Eta: "<<(*simTrack).momentum().eta()<<endl;
        cout<<"Sim Phi: "<<(*simTrack).momentum().phi()<<endl;
	numberOfSimTracks++;
      }    
    } //-------- next simulated track -----------------
    cout << endl;
  } 
  
  reco::TrackCollection::const_iterator staTrack;
  
  cout<<"Reconstructed tracks: " << staTracks->size() << endl;

// loop on reco tracks --------------------------
  for (staTrack = staTracks->begin(); staTrack != staTracks->end(); ++staTrack){
    reco::TransientTrack track(*staTrack,&*theMGField,theTrackingGeometry); 
    
    cout << debug.dumpFTS(track.impactPointTSCP().theState());
    
    recPt = track.impactPointTSCP().momentum().perp();    
    cout<<" p: "<<track.impactPointTSCP().momentum().mag()<< " pT: "<<recPt<<endl;
    cout<<" chi2: "<<track.chi2()<<endl;
    
 // clean ghosts...(it will be fixed in 1_0_2 )
    int ighost = 0;
    if (nmurec > 0 ) { 
      for ( int imu = 1; imu != nmurec+1; ++imu ) {
        float pT1= sqrt(pxrec[imu]*pxrec[imu] + pyrec[imu]*pyrec[imu]);
        float dpT = abs(track.impactPointTSCP().momentum().perp()- pT1 );
	if (dpT < 0.001) ighost = 1;
	cout << " pT1 = " << pT1 <<", dpT = " << dpT << endl;
      }
    }   
    if (ighost == 1) cout << " ***  ghost found *** " << endl;
    if ( ighost == 0)  {
      hPtRec->Fill(recPt);
      nmurec++;
      pxrec[nmurec] = track.impactPointTSCP().momentum().x();
      pyrec[nmurec] = track.impactPointTSCP().momentum().y();
      pzrec[nmurec] = track.impactPointTSCP().momentum().z();    
      etarec[nmurec] = track.impactPointTSCP().momentum().eta(); 
      phirec[nmurec] = track.impactPointTSCP().momentum().phi(); 
      ptrec[nmurec] = track.impactPointTSCP().momentum().perp(); 
      qrec[nmurec] = track.impactPointTSCP().charge(); 
    }  
    TrajectoryStateOnSurface innerTSOS = track.innermostMeasurementState();
    cout << "Inner TSOS:"<<endl;
    cout << debug.dumpTSOS(innerTSOS);
    cout<<" p: "<<innerTSOS.globalMomentum().mag()<< " pT: "<<innerTSOS.globalMomentum().perp()<<endl;

    if(recPt && theDataType == "SimData"){  
      hPres->Fill( (recPt-simPt)/simPt);
      hPtSim->Fill(simPt);
      hPTDiff->Fill(recPt-simPt);
      //      hPTDiff2->Fill(track.innermostMeasurementState().globalMomentum().perp()-simPt);
      hPTDiffvsEta->Fill(track.impactPointTSCP().position().eta(),recPt-simPt);
      hPTDiffvsPhi->Fill(track.impactPointTSCP().position().phi(),recPt-simPt);
      if( ((recPt-simPt)/simPt) <= -0.4)
	cout<<"Out of Res: "<<(recPt-simPt)/simPt<<endl;
      h1_Pres->Fill( ( 1/recPt - 1/simPt)/ (1/simPt));
    }
 cout << " --------- end track ----------------" << endl;    
  }  // next STAmuon reco track ------------------------------


  // analysis ===================================
  hNsim->Fill(nmusim);
  //ADAM
  if(nmusim > 0) hNrec->Fill(nmurec);
  int imuZ1 = 0, imuZ2 = 0;
  
  // simul. invariant Z mass    
  if(nmusim > 1 ) {
    for ( int imu1 = 1; imu1 != nmusim; ++imu1 ) {
      float e1 = sqrt(ptsim[imu1]*ptsim[imu1]+pzsim[imu1]*pzsim[imu1]);
      for ( int imu2 = imu1+1; imu2 != nmusim+1; ++imu2 ) {
	if( qsim[imu1]*qsim[imu2] < 0 ) {
	  float e2 = sqrt(ptsim[imu2]*ptsim[imu2]+pzsim[imu2]*pzsim[imu2]);
	  float pxtot = pxsim[imu1] + pxsim[imu2];
	  float pytot = pysim[imu1] + pysim[imu2];
	  float pztot = pzsim[imu1] + pzsim[imu2];
	  float minv = sqrt( (e1+e2)*(e1+e2) - pxtot*pxtot -  pytot*pytot -  pztot*pztot );
	  hminv->Fill(minv);
	  cout << " minv " << minv << " pt1 " <<  ptsim[imu1] << " pt2 " << ptsim[imu2] << endl;
	  if (abs(minv-91.2) < 5.) {
	    imuZ1 = imu1;
	    imuZ2 = imu2;
	  }
	} // endif opposite charge muons 
      }	// next sim track
    } // next sim track
  } // endif nsim > 1
  
  //  matched tracks...

  //start ADAM
  for ( int imur = 1; imur != nmurec+1; ++imur ) {
    hetarec_all-> Fill(etarec[imur]);
  }
  //end ADAM
  float pxgood[20], pygood[20], pzgood[20],  ptgood[20], qgood[20];
  int isimlink[20], ireclink[20];
  int ngood = 0;  // matched tracks with simulation
  for ( int imu = 1; imu != nmusim+1; ++imu ) {
    hptsim-> Fill(ptsim[imu]);
    if( ptsim[imu] > 5. ) hetasim-> Fill(etasim[imu]);
    float dRmin = 9999.;
    int jbest = 0;  // best reco candidate
    for ( int imur = 1; imur != nmurec+1; ++imur ) {
      float deta = etarec[imur] - etasim[imu];
      float dphi = phirec[imur] - phisim[imu];
      float dR = sqrt(deta*deta + dphi*dphi);
      if(dR < dRmin) {
	dRmin=dR;
	jbest = imur;
      }
    }  // next rec muon
    if( jbest > 0) {
      hdR->Fill(dRmin);
      //ADAM
      if (true || dRmin < 0.2) {
	hptrec-> Fill(ptsim[imu]);
	hptrec1-> Fill(ptrec[jbest]);
	//ADAM
	//float dptrec = (ptrec[jbest]-ptsim[imu])/ptsim[imu];
	float dptrec = (qrec[jbest]/ptrec[jbest]-qsim[imu]/ptsim[imu])/(qsim[imu]/ptsim[imu]);
	hdptrec-> Fill(dptrec);
	hdptrecvsEta-> Fill(etasim[imu],dptrec); 
	if( ptsim[imu] > 5. ) {
	  hetarec-> Fill(etasim[imu]);
	  hetarec1-> Fill(etarec[jbest]);
	}	
	ngood++;  // matched tracks
	pxgood[ngood] = pxrec[jbest];
	pygood[ngood] = pyrec[jbest];
	pzgood[ngood] = pzrec[jbest];
	ptgood[ngood] = ptrec[jbest];
	qgood[ngood] = qrec[jbest];
	isimlink[ngood] = imu;
	ireclink[ngood] = jbest;
      }  // endif dR < 0.3
    } // endif jbest > 0 (found good reco track matched to this sim track 
  } // next sim muon
  
  // plot ghosts (unmatched tracks)...
  if (nmurec > 0 ) {
    for ( int imur = 1; imur != nmurec+1; ++imur ) {
      int imatch = 0;
      if ( ngood > 0 ) {
	for ( int imug = 1; imug != ngood+1; ++imug ) {
	  if( imur == ireclink[imug]) imatch = 1;
	}
      }
      if (imatch == 0)  hptrecvsEta-> Fill(etarec[imur],ptrec[imur]); 
    }  // next rec muon
  } // endif nmurec > 0
  
  //  reconstructed mZ...
  hNmatch->Fill(ngood);
  if(ngood > 1 ) {
    for ( int imu1 = 1; imu1 != ngood; ++imu1 ) {
      float e1 = sqrt(ptgood[imu1]*ptgood[imu1]+pzgood[imu1]*pzgood[imu1]);
      for ( int imu2 = imu1+1; imu2 != ngood+1; ++imu2 ) {
	if( qgood[imu1]*qgood[imu2] < 0 )       {
	  float e2 = sqrt(ptgood[imu2]*ptgood[imu2]+pzgood[imu2]*pzgood[imu2]);
	  float pxtot = pxgood[imu1] + pxgood[imu2];
	  float pytot = pygood[imu1] + pygood[imu2];
	  float pztot = pzgood[imu1] + pzgood[imu2];
	  float minv = sqrt( (e1+e2)*(e1+e2) - pxtot*pxtot -  pytot*pytot -  pztot*pztot );
	  cout << " minvrec " << minv << " pt1 " <<  ptgood[imu1] << " pt2 " << ptgood[imu2] << endl;
	  hminvrec->Fill(minv);
	  if(isimlink[imu1] == imuZ1 || isimlink[imu1] == imuZ2 ) {
	    if(isimlink[imu2] == imuZ1 || isimlink[imu2] == imuZ2 ) hminvrecZ->Fill(minv);
	  } // endif reco tk matched with tk from Z decay
	} // endif opposite charge muons 
      }	// next rec track
    } // next rec track
  } // endif ngood > 1
  
  cout<<"====== end event ================"<<endl;  
}

DEFINE_FWK_MODULE(STAMuonAnalyzer)
