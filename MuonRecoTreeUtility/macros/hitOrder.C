#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
//#include "Math/LorentzVector.h"
#include <vector.h>
#include "RecoTree.h"
//#include "CMS1.h"
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

int hitOrder ( TTree* tree, char *fileName="hitOrder.root", int event = 0) {
  using namespace std;
  
  TFile *histFileAlgo = new TFile(fileName,"RECREATE");

  Init(tree);
  
  TH2F* glb_rHit[10];
  TH2F* glb_zHit[10];
  TH2F* glb_RHit[10];

  for (unsigned int ww=0;ww<10;ww++){
    std::string ss = to_string(ww);
    TString namer("r_");
    TString namez("z_");
    TString nameR("R_");
    TString name2("poorManmu_");
    namer+=ss;
    namez+=ss;
    nameR+=ss;

    glb_rHit[ww] = new TH2F(namer,"r value per hit",76,-0.5,75.5,1000,0.,1000.);
    glb_zHit[ww] = new TH2F(namez,"z value per hit",76,-0.5,75.5,2000,-1000.,1000.);
    glb_RHit[ww] = new TH2F(nameR,"R value per hit",76,-0.5,75.5,5000,0.,5000.);

  }


  tree->GetEntry(event);

  cout << "Event " << EventNumber << " nMu " << nMu << endl;

  for (int iMu=0; iMu < nMu; iMu++) {
    
    for (map<int,std::vector<double> >::const_iterator hitsX = (*l3RecHitsX).begin(); hitsX != (*l3RecHitsX).end(); hitsX++) {        
      for (map<int,std::vector<double> >::const_iterator hitsY = (*l3RecHitsY).begin(); hitsY != (*l3RecHitsY).end(); hitsY++) {          
	for (map<int,std::vector<double> >::const_iterator hitsZ = (*l3RecHitsZ).begin(); hitsZ != (*l3RecHitsZ).end(); hitsZ++) {
	  if (hitsX->first == iMu && hitsX->first == hitsY->first && hitsX->first == hitsZ->first) { //same L2 muon              
	    for (int i = 0; i < hitsX->second.size(); i++) {

	      float r = sqrt(hitsX->second.at(i)*hitsX->second.at(i) + hitsY->second.at(i)*hitsY->second.at(i));
	      float R = sqrt(hitsX->second.at(i)*hitsX->second.at(i) + hitsY->second.at(i)*hitsY->second.at(i) + hitsZ->second.at(i)*hitsZ->second.at(i));

	      if(iMu<10) glb_rHit[iMu]->Fill(i,r);
	      if(iMu<10) glb_zHit[iMu]->Fill(i,hitsZ->second.at(i));
	      if(iMu<10) glb_RHit[iMu]->Fill(i,R);
	    }            
	  }
	} 
      }
    }
        
  }
 
  histFileAlgo->Write("",TObject::kOverwrite);
  histFileAlgo->Close();

  return 0;
 
}
