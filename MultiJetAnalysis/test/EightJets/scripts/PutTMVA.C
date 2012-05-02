#include "TMVA/Reader.h"
using namespace TMVA;
void PutTMVA(TString FileName, TString TRAIN)
{
  TMVA::Tools::Instance();
  // --- Create the Reader object
  TMVA::Reader *reader      = new TMVA::Reader("!Color:!Silent");
  TMVA::Reader *readerJESup = new TMVA::Reader("!Color:!Silent");
  TMVA::Reader *readerJESlo = new TMVA::Reader("!Color:!Silent");
  // Create a set of variables and declare them to the reader
  // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
  Float_t var[20],varJESup[20],varJESlo[20];
  const float JES = 0.03;
  // nominal JES
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
  // JES down
  readerJESlo->AddVariable("pt0 := pt[0]",&varJESlo[0]);
  readerJESlo->AddVariable("pt1 := pt[1]",&varJESlo[1]);
  readerJESlo->AddVariable("pt2 := pt[2]",&varJESlo[2]);
  readerJESlo->AddVariable("pt3 := pt[3]",&varJESlo[3]);
  readerJESlo->AddVariable("pt4 := pt[4]",&varJESlo[4]);
  readerJESlo->AddVariable("pt5 := pt[5]",&varJESlo[5]);
  readerJESlo->AddVariable("pt6 := pt[6]",&varJESlo[6]);
  readerJESlo->AddVariable("pt7 := pt[7]",&varJESlo[7]);
  readerJESlo->AddVariable("b4j := m4jBalance",&varJESlo[8]);
  readerJESlo->AddVariable("b2j := m2jSig/m2jAve",&varJESlo[9]);
  readerJESlo->AddVariable("m2j := m2jAve",&varJESlo[10]);  
  readerJESlo->AddVariable("m4j := m4jAve",&varJESlo[11]); 
  readerJESlo->AddVariable("m8j := m8j",&varJESlo[12]);
  readerJESlo->AddVariable("ht  := ht" ,&varJESlo[13]);
  // JES up
  readerJESup->AddVariable("pt0 := pt[0]",&varJESup[0]);
  readerJESup->AddVariable("pt1 := pt[1]",&varJESup[1]);
  readerJESup->AddVariable("pt2 := pt[2]",&varJESup[2]);
  readerJESup->AddVariable("pt3 := pt[3]",&varJESup[3]);
  readerJESup->AddVariable("pt4 := pt[4]",&varJESup[4]);
  readerJESup->AddVariable("pt5 := pt[5]",&varJESup[5]);
  readerJESup->AddVariable("pt6 := pt[6]",&varJESup[6]);
  readerJESup->AddVariable("pt7 := pt[7]",&varJESup[7]);
  readerJESup->AddVariable("b4j := m4jBalance",&varJESup[8]);
  readerJESup->AddVariable("b2j := m2jSig/m2jAve",&varJESup[9]);
  readerJESup->AddVariable("m2j := m2jAve",&varJESup[10]);  
  readerJESup->AddVariable("m4j := m4jAve",&varJESup[11]); 
  readerJESup->AddVariable("m8j := m8j",&varJESup[12]);
  readerJESup->AddVariable("ht  := ht" ,&varJESup[13]);

  // --- Book the MVA methods
  //reader->BookMVA("LIK","weights/factory_"+TRAIN+"_LikelihoodD.weights.xml");
  reader->BookMVA("FIS","weights/factory_"+TRAIN+"_Fisher.weights.xml");
  reader->BookMVA("BDT","weights/factory_"+TRAIN+"_BDT.weights.xml");
  reader->BookMVA("MLP","weights/factory_"+TRAIN+"_MLP_ANN.weights.xml");
  readerJESlo->BookMVA("FIS_JESlo","weights/factory_"+TRAIN+"_Fisher.weights.xml");
  readerJESlo->BookMVA("BDT_JESlo","weights/factory_"+TRAIN+"_BDT.weights.xml");
  readerJESlo->BookMVA("MLP_JESlo","weights/factory_"+TRAIN+"_MLP_ANN.weights.xml");
  readerJESup->BookMVA("FIS_JESup","weights/factory_"+TRAIN+"_Fisher.weights.xml");
  readerJESup->BookMVA("BDT_JESup","weights/factory_"+TRAIN+"_BDT.weights.xml");
  readerJESup->BookMVA("MLP_JESup","weights/factory_"+TRAIN+"_MLP_ANN.weights.xml");
  float pt[8],ht4j[2],m4j[2],m8j,ht,m2jAve,m4jAve,m2jSig,m4jBalance;

  TFile *inf   = TFile::Open(FileName+".root");
  TFile *outf  = TFile::Open(FileName+"_tmva"+TRAIN+".root","RECREATE");
  TTree *trIN  = (TTree*)inf->Get("multijets/events");
  TTree *trOUT = (TTree*)trIN->CloneTree();

  float BDT,MLP,FIS,LIK,BDT_JESlo,MLP_JESlo,FIS_JESlo,LIK_JESlo,BDT_JESup,MLP_JESup,FIS_JESup,LIK_JESup;
  //TBranch *brLIK = trOUT->Branch("LikelihoodD",&LIK,"LIK/F");
  TBranch *brFIS = trOUT->Branch("Fisher",&FIS,"FIS/F");
  TBranch *brBDT = trOUT->Branch("BDT",&BDT,"BDT/F");
  TBranch *brMLP = trOUT->Branch("MLP",&MLP,"MLP/F");
  TBranch *brFIS_JESup = trOUT->Branch("Fisher_JESup",&FIS_JESup,"FIS_JESup/F");
  TBranch *brBDT_JESup = trOUT->Branch("BDT_JESup",&BDT_JESup,"BDT_JESup/F");
  TBranch *brMLP_JESup = trOUT->Branch("MLP_JESup",&MLP_JESup,"MLP_JESup/F");
  TBranch *brFIS_JESlo = trOUT->Branch("Fisher_JESlo",&FIS_JESlo,"FIS_JESlo/F");
  TBranch *brBDT_JESlo = trOUT->Branch("BDT_JESlo",&BDT_JESlo,"BDT_JESlo/F");
  TBranch *brMLP_JESlo = trOUT->Branch("MLP_JESlo",&MLP_JESlo,"MLP_JESlo/F");
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
    
    varJESlo[0]  = (1-JES)*pt[0];
    varJESlo[1]  = (1-JES)*pt[1];
    varJESlo[2]  = (1-JES)*pt[2];
    varJESlo[3]  = (1-JES)*pt[3];
    varJESlo[4]  = (1-JES)*pt[4];
    varJESlo[5]  = (1-JES)*pt[5];
    varJESlo[6]  = (1-JES)*pt[6];
    varJESlo[7]  = (1-JES)*pt[7];
    varJESlo[8]  = m4jBalance;
    varJESlo[9]  = m2jSig/m2jAve;
    varJESlo[10] = (1-JES)*m2jAve;
    varJESlo[11] = (1-JES)*m4jAve;
    varJESlo[12] = (1-JES)*m8j;
    varJESlo[13] = (1-JES)*ht;
    
    varJESup[0]  = (1+JES)*pt[0];
    varJESup[1]  = (1+JES)*pt[1];
    varJESup[2]  = (1+JES)*pt[2];
    varJESup[3]  = (1+JES)*pt[3];
    varJESup[4]  = (1+JES)*pt[4];
    varJESup[5]  = (1+JES)*pt[5];
    varJESup[6]  = (1+JES)*pt[6];
    varJESup[7]  = (1+JES)*pt[7];
    varJESup[8]  = m4jBalance;
    varJESup[9]  = m2jSig/m2jAve;
    varJESup[10] = (1+JES)*m2jAve;
    varJESup[11] = (1+JES)*m4jAve;
    varJESup[12] = (1+JES)*m8j;
    varJESup[13] = (1+JES)*ht;
    
    //LIK = reader->EvaluateMVA("LIK");
    FIS = reader->EvaluateMVA("FIS");
    BDT = reader->EvaluateMVA("BDT");
    MLP = reader->EvaluateMVA("MLP");
    FIS_JESlo = readerJESlo->EvaluateMVA("FIS_JESlo");
    BDT_JESlo = readerJESlo->EvaluateMVA("BDT_JESlo");
    MLP_JESlo = readerJESlo->EvaluateMVA("MLP_JESlo");
    FIS_JESup = readerJESup->EvaluateMVA("FIS_JESup");
    BDT_JESup = readerJESup->EvaluateMVA("BDT_JESup");
    MLP_JESup = readerJESup->EvaluateMVA("MLP_JESup");

    //brLIK->Fill();
    brFIS->Fill();
    brBDT->Fill();
    brMLP->Fill();
    brFIS_JESlo->Fill();
    brBDT_JESlo->Fill();
    brMLP_JESlo->Fill();
    brFIS_JESup->Fill();
    brBDT_JESup->Fill();
    brMLP_JESup->Fill();
  }
  TDirectoryFile *dir = (TDirectoryFile*)outf->mkdir("multijets");
  dir->cd();
  trOUT->Write("events");
  outf->Close();
  inf->Close();
}
