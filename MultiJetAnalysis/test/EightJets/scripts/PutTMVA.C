#include "TMVA/Reader.h"
using namespace TMVA;
void PutTMVA(TString FileName, TString SIGNAL)
{
  TMVA::Tools::Instance();
  // --- Create the Reader object
  TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent");
  // Create a set of variables and declare them to the reader
  // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
  Float_t var[20];
  
  reader->AddVariable("pt0 := pt[0]",&var[0]);
  reader->AddVariable("pt1 := pt[1]",&var[1]);
  reader->AddVariable("pt2 := pt[2]",&var[2]);
  reader->AddVariable("pt3 := pt[3]",&var[3]);
  reader->AddVariable("pt4 := pt[4]",&var[4]);
  reader->AddVariable("pt5 := pt[5]",&var[5]);
  reader->AddVariable("pt6 := pt[6]",&var[6]);
  reader->AddVariable("pt7 := pt[7]",&var[7]);
  reader->AddVariable("b4j := m4jBalance",&var[8]);
  reader->AddVariable("b2j := m2jSig/m2jAve",&var[9]);
  reader->AddVariable("m2j := m2jAve",&var[10]);  
  reader->AddVariable("m4j := m4jAve",&var[11]); 
  reader->AddVariable("m8j := m8j",&var[12]);
  reader->AddVariable("ht  := ht" ,&var[13]);

  // --- Book the MVA methods
  reader->BookMVA("FIS","weights/factory_"+SIGNAL+"_Fisher.weights.xml");
  reader->BookMVA("LIK","weights/factory_"+SIGNAL+"_LikelihoodD.weights.xml");
  reader->BookMVA("BDT","weights/factory_"+SIGNAL+"_BDT.weights.xml");
  reader->BookMVA("MLP","weights/factory_"+SIGNAL+"_MLP_ANN.weights.xml");
  
  float pt[8],ht4j[2],m4j[2],m8j,ht,m2jAve,m4jAve,m2jSig,m4jBalance;

  TFile *inf   = TFile::Open(FileName+".root");
  TFile *outf  = TFile::Open(FileName+"_tmva.root","RECREATE");
  TTree *trIN  = (TTree*)inf->Get("multijets/events");
  TTree *trOUT = (TTree*)trIN->CloneTree();

  float BDT,MLP,FIS,LIK;
  TBranch *brFIS = trOUT->Branch("Fisher",&FIS,"FIS/F");
  TBranch *brLIK = trOUT->Branch("LikelihoodD",&LIK,"LIK/F");
  TBranch *brBDT = trOUT->Branch("BDT",&BDT,"BDT/F");
  TBranch *brMLP = trOUT->Branch("MLP",&MLP,"MLP/F");
  for(int i=0;i<trOUT->GetEntries();i++) {
    trOUT->GetEntry(i);
    trOUT->SetBranchAddress("pt",&pt);
    trOUT->SetBranchAddress("ht",&ht);
    trOUT->SetBranchAddress("ht4j",&ht4j);
    trOUT->SetBranchAddress("m8j",&m8j);
    trOUT->SetBranchAddress("m4j",&m4j);
    trOUT->SetBranchAddress("m2jAve",&m2jAve);
    trOUT->SetBranchAddress("m4jAve",&m4jAve);
    trOUT->SetBranchAddress("m2jSig",&m2jSig);
    trOUT->SetBranchAddress("m4jBalance",&m4jBalance);
    
    var[0]  = pt[0];
    var[1]  = pt[1];
    var[2]  = pt[2];
    var[3]  = pt[3];
    var[4]  = pt[4];
    var[5]  = pt[5];
    var[6]  = pt[6];
    var[7]  = pt[7];
    var[8]  = m4jBalance;
    var[9]  = m2jSig/m2jAve;
    var[10] = m2jAve;
    var[11] = m4jAve;
    var[12] = m8j;
    var[13] = ht;
    
    FIS = reader->EvaluateMVA("FIS");
    LIK = reader->EvaluateMVA("LIK");
    BDT = reader->EvaluateMVA("BDT");
    MLP = reader->EvaluateMVA("MLP");

    brFIS->Fill();
    brLIK->Fill();
    brBDT->Fill();
    brMLP->Fill();
  }
  TDirectoryFile *dir = (TDirectoryFile*)outf->mkdir("multijets");
  dir->cd();
  trOUT->Write("events");
  outf->Close();
  inf->Close();
}
