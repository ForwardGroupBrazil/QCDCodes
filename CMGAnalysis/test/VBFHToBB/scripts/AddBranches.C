#include "TMVA/Reader.h"
#include "TFile.h"
#include "TKey.h"
#include "TCollection.h"
#include "TLorentzVector.h"
#include "Settings.h"
using namespace TMVA;
using namespace TMath;
void AddBranches()
{
  float varCLA[50],varREG[50],varQGL[5];
  TMVA::Reader *readerQGL[7][3];
  cout<<"Booking QGL readers"<<endl;
  for(int ipt=0;ipt<7;ipt++) {
    for(int ieta=0;ieta<3;ieta++) {
      readerQGL[ipt][ieta] = new TMVA::Reader("!Color:!Silent");
      readerQGL[ipt][ieta]->AddVariable("axis1"  ,&varQGL[0]);
      readerQGL[ipt][ieta]->AddVariable("axis2"  ,&varQGL[1]);
      readerQGL[ipt][ieta]->AddVariable("Mult"   ,&varQGL[2]);
      readerQGL[ipt][ieta]->AddVariable("JetR"   ,&varQGL[3]);
      readerQGL[ipt][ieta]->AddVariable("JetPull",&varQGL[4]);
      TString ss = PT_CAT[ipt]+"_"+ETA_CAT[ieta]+"_Likelihood.xml";
      readerQGL[ipt][ieta]->BookMVA("LIK_"+PT_CAT[ipt]+"_"+ETA_CAT[ieta],"/afs/cern.ch/work/k/kkousour/private/data/vbfhbb/xml_files/"+ss);
    }
  }
  
  TMVA::Reader *readerCLA    = new TMVA::Reader("!Color:!Silent");
  // ---- classification TMVA -----------------------
  readerCLA->AddVariable("dEtaqq"           ,&varCLA[0]);
  readerCLA->AddVariable("dEtaqqEta-dEtaqq" ,&varCLA[1]);
  readerCLA->AddVariable("mqq"              ,&varCLA[2]);
  readerCLA->AddVariable("etaBoostqq"       ,&varCLA[3]);
  readerCLA->AddVariable("jetBtag[0]"       ,&varCLA[4]);
  readerCLA->AddVariable("jetBtag[1]"       ,&varCLA[5]);
  readerCLA->AddVariable("jetQGLnew[2]"     ,&varCLA[6]);
  readerCLA->AddVariable("jetQGLnew[3]"     ,&varCLA[7]);
  readerCLA->AddVariable("softHt"           ,&varCLA[8]);
  readerCLA->AddVariable("softMulti"        ,&varCLA[9]);
  readerCLA->AddVariable("jetEtaBtag[2]"    ,&varCLA[10]);
  readerCLA->AddVariable("jetPt[4]/mqq"     ,&varCLA[11]);
  // spectator variables: not used for the training but recoreded
  readerCLA->AddSpectator("mbb"             ,&varCLA[12]);
  readerCLA->AddSpectator("dPhibb"          ,&varCLA[13]);
  readerCLA->AddSpectator("rho"             ,&varREG[14]);

  // --- Book the MVA methods
  readerCLA->BookMVA("BDT_CLA","/afs/cern.ch/work/k/kkousour/private/data/vbfhbb/xml_files/factory_hybrid21_125_BDT.weights.xml");
  readerCLA->BookMVA("MLP_CLA","/afs/cern.ch/work/k/kkousour/private/data/vbfhbb/xml_files/factory_hybrid21_125_MLP_ANN.weights.xml");
  // ---- regression TMVA ---------------------------
  TMVA::Reader *readerREG = new TMVA::Reader("!Color:!Silent");
  
  readerREG->AddVariable("jetBtag"   ,&varREG[0]);
  readerREG->AddVariable("jetPt"     ,&varREG[1]);
  readerREG->AddVariable("jetEta"    ,&varREG[2]);
  readerREG->AddVariable("jetMetPhi" ,&varREG[3]);
  readerREG->AddVariable("jetChf"    ,&varREG[4]);
  readerREG->AddVariable("jetPhf"    ,&varREG[5]);
  readerREG->AddVariable("jetNhf"    ,&varREG[6]);
  readerREG->AddVariable("jetElf"    ,&varREG[7]);
  readerREG->AddVariable("jetMuf"    ,&varREG[8]);
  readerREG->AddVariable("jetPtD"    ,&varREG[9]);
  readerREG->AddVariable("jetVtxPt"  ,&varREG[10]);
  readerREG->AddVariable("jetVtx3dL" ,&varREG[11]);
  readerREG->AddVariable("jetVtx3deL",&varREG[12]);
  readerREG->AddVariable("met"       ,&varREG[13]);
  readerREG->AddVariable("rho"       ,&varREG[14]);
    
  readerREG->BookMVA("BDT_REG","/afs/cern.ch/work/k/kkousour/private/data/vbfhbb/xml_files/factoryJetReg_BDT.weights.xml");
  
  float dEtaqq,dEtaqqEta,dPhibb,mqq,mbb,softHt,etaBoostqq;
  int btagIdx[5],jetPart[5],softMulti;
  float jetQGLnew[5],jetBtag[5],jetPt[5],jetPtD[5],jetEta[5],jetPhi[5],jetMass[5],jetChf[5],jetNhf[5],jetPhf[5],jetElf[5],jetMuf[5],jetVtxPt[5],jetVtx3dL[5],jetVtx3deL[5],jetAxis[2][5],jetAxis_QC[2][5],jetPull[5],jetPull_QC[5],jetR[5],jetnChg_QC[5],jetnChg_ptCut[5],jetnNeutral_ptCut[5];
  float met,metPhi,rho,mbbCor;
  float BDT_CLA,MLP_CLA;
  
  TString FileName[19] = {
    "flatTree_GluGlu-Powheg125_preselect_hard",
    "flatTree_GluGlu-Madgraph125_preselect_hard",
    "flatTree_VBF-Powheg115_preselect_hard",
    "flatTree_VBF-Powheg120_preselect_hard",
    "flatTree_VBF-Powheg125_preselect_hard",
    "flatTree_VBF-Powheg130_preselect_hard",
    "flatTree_VBF-Powheg135_preselect_hard",
    "flatTree_QCD-HT100_preselect_hard",
    "flatTree_QCD-HT250_preselect_hard",
    "flatTree_QCD-HT500_preselect_hard",
    "flatTree_QCD-HT1000_preselect_hard",
    "flatTree_ZJets_preselect_hard",
    "flatTree_T_preselect_hard",
    "flatTree_Tbar_preselect_hard",
    "flatTree_TTJets_preselect_hard",
    "flatTree_MultiJet-Run2012A_preselect_hard",
    "flatTree_BJetPlusX-Run2012B_preselect_hard",
    "flatTree_BJetPlusX-Run2012C_preselect_hard",
    "flatTree_BJetPlusX-Run2012D_preselect_hard"
  };
  TKey *key;
  TString PATH("/afs/cern.ch/work/k/kkousour/private/data/vbfhbb/");
  for(int iFile=16;iFile<19;iFile++) {
    cout<<"Opening file: "<<FileName[iFile]<<endl;
    TFile *inf   = TFile::Open(PATH+FileName[iFile]+".root");
    TFile *outf  = TFile::Open(PATH+FileName[iFile]+"_tmva.root","RECREATE");
    // find the directories
    std::vector<TString> DIR_NAME;
    TIter nextkey(inf->GetListOfKeys());
    while ((key = (TKey*)nextkey())) {
      if (TString(key->GetClassName()).CompareTo("TDirectoryFile") == 0) {
        DIR_NAME.push_back(TString(key->GetName()));
      }  
    }
    for(unsigned int idir=0;idir<DIR_NAME.size();idir++) {
      TTree *trIN  = (TTree*)inf->Get(DIR_NAME[idir]+"/events");
      cout<<"Cloning tree "<<DIR_NAME[idir]<<endl;
      TTree *trOUT = (TTree*)trIN->CloneTree(-1,"fast");
      cout<<"Filling new tree"<<endl;  
      TBranch *brBDT    = trOUT->Branch("BDT"      ,&BDT_CLA    ,"BDT_CLA/F");
      TBranch *brMLP    = trOUT->Branch("MLP"      ,&MLP_CLA    ,"MLP_CLA/F");
      TBranch *brMbbCOR = trOUT->Branch("mbbCor"   ,&mbbCor     ,"mbbCor/F");
      TBranch *brQGLnew = trOUT->Branch("jetQGLnew",&jetQGLnew  ,"jetQGLnew[5]/F");
    
      trOUT->SetBranchAddress("btagIdx"          ,&btagIdx);
      trOUT->SetBranchAddress("mqq"              ,&mqq);
      trOUT->SetBranchAddress("mbb"              ,&mbb);
      trOUT->SetBranchAddress("dEtaqq"           ,&dEtaqq);
      trOUT->SetBranchAddress("dEtaqqEta"        ,&dEtaqqEta);
      trOUT->SetBranchAddress("dPhibb"           ,&dPhibb);
      trOUT->SetBranchAddress("etaBoostqq"       ,&etaBoostqq);
      trOUT->SetBranchAddress("nSoftTrackJets"   ,&softMulti);
      trOUT->SetBranchAddress("softHt"           ,&softHt);
      trOUT->SetBranchAddress("jetPt"            ,&jetPt);
      trOUT->SetBranchAddress("jetPtD"           ,&jetPtD);
      trOUT->SetBranchAddress("jetPart"          ,&jetPart);
      trOUT->SetBranchAddress("jetEta"           ,&jetEta);
      trOUT->SetBranchAddress("jetPhi"           ,&jetPhi);
      trOUT->SetBranchAddress("jetBtag"          ,&jetBtag);
      trOUT->SetBranchAddress("jetMass"          ,&jetMass);
      trOUT->SetBranchAddress("jetChf"           ,&jetChf);
      trOUT->SetBranchAddress("jetNhf"           ,&jetNhf);
      trOUT->SetBranchAddress("jetPhf"           ,&jetPhf);
      trOUT->SetBranchAddress("jetElf"           ,&jetElf);
      trOUT->SetBranchAddress("jetMuf"           ,&jetMuf);
      trOUT->SetBranchAddress("jetVtxPt"         ,&jetVtxPt); 
      trOUT->SetBranchAddress("jetVtx3dL"        ,&jetVtx3dL);
      trOUT->SetBranchAddress("jetVtx3deL"       ,&jetVtx3deL); 
      trOUT->SetBranchAddress("jetAxis"          ,&jetAxis);
      trOUT->SetBranchAddress("jetAxis_QC"       ,&jetAxis_QC);
      trOUT->SetBranchAddress("jetPull"          ,&jetPull);
      trOUT->SetBranchAddress("jetPull_QC"       ,&jetPull_QC);
      trOUT->SetBranchAddress("jetR"             ,&jetR);
      trOUT->SetBranchAddress("jetnChg_QC"       ,&jetnChg_QC);
      trOUT->SetBranchAddress("jetnChg_ptCut"    ,&jetnChg_ptCut);
      trOUT->SetBranchAddress("jetnNeutral_ptCut",&jetnNeutral_ptCut);
      trOUT->SetBranchAddress("rho"              ,&rho);
      trOUT->SetBranchAddress("met"              ,&met);
      trOUT->SetBranchAddress("metPhi"           ,&metPhi);
  
      int decade(0);
      int NN = trOUT->GetEntries();
      cout<<"Reading "<<NN<<" entries"<<endl;
      for(int i=0;i<trOUT->GetEntries();i++) {
        double progress = 10.0*i/(1.0*NN);
        int k = TMath::FloorNint(progress); 
        if (k > decade) 
          cout<<10*k<<" %"<<endl;
        decade = k;
        trOUT->GetEntry(i);

        for(int j=0;j<5;j++) { 
          float pt  = jetPt[j];
          float eta = fabs(jetEta[j]);
          int ipt  = TMath::Max(0,FindIndex(8,PT_QGL,pt));
          int ieta = TMath::Max(0,FindIndex(4,ETA_QGL,eta));
          if (ieta == 0) {
            varQGL[0] = jetAxis_QC[0][j]-cor_rho[ipt][ieta][0]*rho;
            varQGL[1] = jetAxis_QC[1][j]-cor_rho[ipt][ieta][1]*rho;
            varQGL[2] = jetnChg_QC[j]-cor_rho[ipt][ieta][3]*rho;
            varQGL[3] = jetR[j]-cor_rho[ipt][ieta][4]*rho;
            varQGL[4] = jetPull_QC[j]-cor_rho[ipt][ieta][2]*rho;
          }
          else {
            varQGL[0] = jetAxis[0][j]-cor_rho[ipt][ieta][0]*rho;
            varQGL[1] = jetAxis[1][j]-cor_rho[ipt][ieta][1]*rho;
            varQGL[2] = jetnChg_ptCut[j]+jetnNeutral_ptCut[j]-cor_rho[ipt][ieta][3]*rho;
            varQGL[3] = jetR[j]-cor_rho[ipt][ieta][4]*rho;
            varQGL[4] = jetPull[j]-cor_rho[ipt][ieta][2]*rho;
          }
          jetQGLnew[j] = readerQGL[ipt][ieta]->EvaluateMVA("LIK_"+PT_CAT[ipt]+"_"+ETA_CAT[ieta]);
        }

        varCLA[0]  = dEtaqq;
        varCLA[1]  = dEtaqqEta-dEtaqq;
        varCLA[2]  = mqq;
        varCLA[3]  = etaBoostqq;
        varCLA[4]  = jetBtag[btagIdx[0]];
        varCLA[5]  = jetBtag[btagIdx[1]];
        varCLA[6]  = jetQGLnew[btagIdx[2]];
        varCLA[7]  = jetQGLnew[btagIdx[3]];
        varCLA[8]  = softHt;
        varCLA[9]  = softMulti;
        varCLA[10] = jetEta[btagIdx[2]];
        varCLA[11] = jetPt[4]/mqq;
        //--- spectator vatiables -----
        varCLA[12] = mbb;
        varCLA[13] = dPhibb;
  
        BDT_CLA    = readerCLA->EvaluateMVA("BDT_CLA");
        MLP_CLA    = readerCLA->EvaluateMVA("MLP_CLA");

        brBDT->Fill();
        brMLP->Fill();
        brQGLnew->Fill();
    
        TLorentzVector P4[2];
        for(int j=0;j<2;j++) {
          int ib = btagIdx[j];
          varREG[0]  = jetBtag[ib];
          varREG[1]  = jetPt[ib]; 
          varREG[2]  = jetEta[ib];
          float tmpPhi = fabs(jetPhi[ib]-metPhi);
          if (tmpPhi > 3.14159) {
            tmpPhi -= 3.14159;
          }
          varREG[3]  = tmpPhi;
          varREG[4]  = jetChf[ib];
          varREG[5]  = jetPhf[ib];
          varREG[6]  = jetNhf[ib];
          varREG[7]  = jetElf[ib];
          varREG[8]  = jetMuf[ib];
          varREG[9]  = jetPtD[ib];
          varREG[10] = jetVtxPt[ib];
          varREG[11] = jetVtx3dL[ib];
          varREG[12] = jetVtx3deL[ib];
          varREG[13] = met;
          varREG[14] = rho;
          float tmvaREG = readerREG->EvaluateRegression("BDT_REG")[0];
          float cor = tmvaREG/jetPt[ib];
          P4[j].SetPtEtaPhiM(cor*jetPt[ib],jetEta[ib],jetPhi[ib],cor*jetMass[ib]);
        }
        TLorentzVector corP4(P4[0]+P4[1]);
        mbbCor = corP4.M();
        brMbbCOR->Fill();
      }
      TDirectoryFile *dir = (TDirectoryFile*)outf->mkdir(DIR_NAME[idir]);
      dir->cd();
      trOUT->Write("events");
    }// end of directory loop  
    outf->Close();
    inf->Close();
  }// end of file loop  
}
