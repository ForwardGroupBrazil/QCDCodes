#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
//#include "Math/LorentzVector.h"
#include <vector.h>
#include "RecoTree.h"
#include <map.h>
#include "TMath.h"
#include <sstream>

template <class T>
inline std::string to_string (const T& t)
{
std::stringstream ss;
ss << t;
return ss.str();
}

int oneKink ( TTree* tree, char *fileName="kink.root", int event = 0) {
  using namespace std;
  
  TFile *histFileAlgo = new TFile(fileName,"RECREATE");

  Init(tree);
  
  TH3F* glb_poorMan[10];
  TH2F* glb_kinkHit[10];

  for (unsigned int ww=0;ww<10;ww++){
    std::string ss = to_string(ww);
    TString name("mu_");
    TString name2("poorManmu_");
    name+=ss;
    name2+=ss;

    glb_kinkHit[ww] = new TH2F(name,"Kink value per hit",76,-0.5,75.5,100,0.,100.);
    glb_poorMan[ww] = new TH3F(name2,"test for hits", 200, -800.0, 800.0, 200, -800.0, 800.0, 250, -1045.0, 1045.0); 
  }


  TH1F *recHitsTest = new TH1F("recHitsTest","test to fill x position of recHits",1600,-800,800); 
  //  TH3F *poorMansDisplay = new TH3F("poorMansDisplay","test for hits", 200, -800.0, 800.0, 200, -800.0, 800.0, 250, -1045.0, 1045.0); 
  //  TH2F *kinkPerHit = new TH2F("kinkPerHit","Kink value per hit",76,-0.5,75.5,100,0.,10.);

  tree->GetEntry(event);

  cout << "nMu " << nMu << endl;

  for (int iMu=0; iMu < nMu; iMu++) {
    
    for (map<int,std::vector<double> >::const_iterator hitsX = (*l3RecHitsX).begin(); hitsX != (*l3RecHitsX).end(); hitsX++) {        
      for (map<int,std::vector<double> >::const_iterator hitsY = (*l3RecHitsY).begin(); hitsY != (*l3RecHitsY).end(); hitsY++) {          
	for (map<int,std::vector<double> >::const_iterator hitsZ = (*l3RecHitsZ).begin(); hitsZ != (*l3RecHitsZ).end(); hitsZ++) {
	  if (hitsX->first == hitsY->first && hitsX->first == hitsZ->first) { //same L2 muon              
	    for (int i = 0; i < hitsX->second.size(); i++) {
	      if(iMu<10) glb_poorMan[iMu]->Fill(hitsX->second.at(i),hitsY->second.at(i),hitsZ->second.at(i));
	    }            
	  }
	} 
      }
    }
        
    double resultGlb = 0.;
    for (map<int,std::vector<double> >::const_iterator ii = (*l3RecHitsPhiTM).begin(); ii != (*l3RecHitsPhiTM).end(); ii++) {
      for (map<int,std::vector<double> >::const_iterator jj = (*l3RecHitsPhiTSOS).begin(); jj != (*l3RecHitsPhiTSOS).end(); jj++) {
	for (map<int,std::vector<double> >::const_iterator kk = (*l3RecHitsErrorTM).begin(); kk != (*l3RecHitsErrorTM).end(); kk++) {
	  if (ii->first == iMu  && 
	      ii->first == jj->first && 
	      ii->first == kk->first) {
	    for (int i = 0; i < ii->second.size(); i++) {
	      double diff = fabs(jj->second.at(i) - ii->second.at(i));
	      if ( diff > TMath::Pi() ) diff = 2*TMath::Pi() - diff;
	      double error = kk->second.at(i);
	      double s = ( error > 0. ) ? (diff*diff)/error : (diff*diff);
	      if(iMu<10 && error>0) glb_kinkHit[iMu]->Fill(i,s);
	      resultGlb += s;
	    }
	  }
	}
      } 
    }
    float kinkVal = (*muGlbKink).at(iMu);
    cout << "iMu " << iMu << " Bit  " << (*l3AssociationMyBit).at(iMu) << " Saved val " << kinkVal << " and new val " << resultGlb << endl;
  }
 
  histFileAlgo->Write("",TObject::kOverwrite);
  histFileAlgo->Close();

  return 0;
 
}
