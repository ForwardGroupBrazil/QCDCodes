static Int_t nE(0), nJ(0);
#include <stdio.h>
#include <time.h>

 
TChain* GetTX(const char* title, int skip){

  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();
  
  nE = skip;
  ifstream file;
  file.open(title);

  Int_t newline=0;
  char c;
  char line[1000];  

  TChain *chain = new TChain("Events");

  while( file.getline(line, 1000)){
    newline++;
    if(newline<skip) continue;
    cout << "read line: " << line << endl;
    TFile *f;
    TString name(line);
    Each(chain,f,name,newline);
    
  }
  
  cout << " =============================== " << endl;
  cout << "number of events = " << nE <<endl;
  cout << " =============================== " << endl; 

  return chain;
  
}

void Each(TChain* chain, TFile* g, TString& s,int newline){
  
  ofstream file("thefile",ios::app);
  file <<  s.Data();
  file.close(); 
  nE++;  
  cout << nE << " " << s << endl;

  //TFileCollection fc("dum");
  //TChain chain("Events");

  try {
    time_t start;
    time (&start);
    // DCACHE:  "root://dcache-00.rcac.purdue.edu/pnfs/rcac.purdue.edu/data/";
    // CASTOR : "rfio:/castor/cern.ch/cms/";
    //    TString prestuff = "root://dcache-00.rcac.purdue.edu/pnfs/rcac.purdue.edu/data/";
    TString prestuff = "dcap://dcache.rcac.purdue.edu:22125/pnfs/rcac.purdue.edu/data";
    //cout << prestuff+s << endl;
    
    TXNetFile::AsyncOpen(prestuff+s);

    Bool_t too_long=false;
      

    int j=0;
    while (TXNetFile::GetAsyncOpenStatus(prestuff+s) == TXNetFile::kAOSInProgress){
      time_t here;
      time (&here);

      if( difftime(here,start)>j ){ 
	j++;
	cout << j << "...";
      }
 
     if(difftime(here,start) > 25){
	too_long=true;
	ofstream file("thefile",ios::app);
	file << "      TOOLONG" << '\n';
	file.close();
	break;
     }
    }
    
    //cout << " seconds, too long?  " << too_long << endl;
   


    if(!too_long){
      //fc.Add(new TFileInfo(prestuff+s));
	  chain->Add(prestuff+s);
    }
    
  }
  //return fc;
}
