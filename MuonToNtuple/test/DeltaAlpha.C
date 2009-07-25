#include "TMath.h"

void DeltaAlpha( void ) 
{
    gStyle->SetOptStat(0);
    gROOT->ProcessLine(".L setTDRStyle.C");
    gROOT->LoadMacro("adamFunctions.C");
    setTDRStyle();
    gROOT->SetStyle("tdrStyle");
    gROOT->ForceStyle(true);

    TCanvas *c1 = new TCanvas("dAlpha", "dAlpha", 800, 600);

    TString _fsig = "test";
    TString _fbg2 = "test";
    TString _fbg3 = ""; //"QCD_BCtoMu_Pt30to50";
    TString _fbg4 = ""; //"InclusivePPmuX";
    TString _fbg5 = "";
    TH1D *h1, *h2, *h3, *h4, *h5;
    h1 = GetDistance("../rootfiles/"+_fsig+".root"); 
    c1->Update();
    h2 = GetDistance("../rootfiles/"+_fbg2+".root"); 
    c1->Update();
    if( _fbg3 != "" ) {
	h3 = GetDistance("../rootfiles/"+_fbg3+".root"); 
        c1->Update();
    }
    if( _fbg4 != "" ) {
	h4 = GetDistance("../rootfiles/"+_fbg4+".root"); 
        c1->Update();
    }
    if( _fbg5 != "" ) {
	h5 = GetDistance("../rootfiles/"+_fbg5+".root"); 
        c1->Update();
    }

    c1->SetLogy();
    c1->cd();
    h1->SetTitle("#Delta(#alpha) between two muons");
    h1->GetXaxis()->SetTitle("#Delta(#eta)+#Delta(#phi)-2#pi");
    h1->DrawNormalized();
    h2->SetLineColor(2);
    h2->DrawNormalized("same");
    if( _fbg3 != "" ) {
        h3->SetLineColor(4);
        h3->DrawNormalized("same");
    }
    if( _fbg4 != "" ) {
        h4->SetLineColor(6);
        h4->DrawNormalized("same");
    }
    if( _fbg5 != "" ) {
        h5->SetLineColor(6);
        h5->DrawNormalized("same");
    }

    TLegend *legend = new TLegend(0.6, 0.7, 0.95, 0.95);
    legend->AddEntry(h1, _fsig, "l");
    legend->AddEntry(h2, _fbg2, "l");
    if( _fbg3 != "" ) legend->AddEntry(h3, _fbg3, "l");
    if( _fbg4 != "" ) legend->AddEntry(h4, _fbg4, "l");
    if( _fbg5 != "" ) legend->AddEntry(h5, _fbg5, "l");
    legend->Draw();
    printCanvasesType(".eps");//c1->Print("./plots/Dist_"+_fbg2+".eps");
    //c1->Print("./plots/Dist_"+_fbg2+".pdf");
    //c1->Print("./plots/Dist_"+_fbg2+".gif");
}

TH1D* GetDistance( TString fname )
{
    const double par_mass = 0.105658;
    const int num = 500;

    TFile *f = new TFile(fname); 
    TTree *t1 = (TTree*)f->Get("DiMuonTree");

    int npair;
    int hlt_mu7;
    int hlt_l1mu;
    int sign[num];
    double vtxChi2[num];
    double reco_mass[num];
    int muon1_type[num];
    int muon2_type[num];
    double muon1_pt[num];
    double muon2_pt[num];
    double muon1_eta[num];
    double muon2_eta[num];
    double muon1_phi[num];
    double muon2_phi[num];
    double muon1_trkiso[num];
    double muon2_trkiso[num];

    double muon1_vx[num];
    double muon2_vx[num];
    double muon1_vy[num];
    double muon2_vy[num];
    double muon1_vz[num];
    double muon2_vz[num];

    t1->SetBranchAddress("nPair", &npair);
    t1->SetBranchAddress("InvMass", &reco_mass);
    //t1->SetBranchAddress("vtxChi2", &vtxChi2);
    t1->SetBranchAddress("HLT_L1Mu", &hlt_l1mu);
    t1->SetBranchAddress("HLT_Mu7", &hlt_mu7);
    t1->SetBranchAddress("isOppSign", &sign);
    t1->SetBranchAddress("Muon1_muonType", &muon1_type);
    t1->SetBranchAddress("Muon2_muonType", &muon2_type);
    t1->SetBranchAddress("Muon1_pT", &muon1_pt);
    t1->SetBranchAddress("Muon2_pT", &muon2_pt);
    t1->SetBranchAddress("Muon1_eta", &muon1_eta);
    t1->SetBranchAddress("Muon2_eta", &muon2_eta);
    t1->SetBranchAddress("Muon1_phi", &muon1_phi);
    t1->SetBranchAddress("Muon2_phi", &muon2_phi);
    t1->SetBranchAddress("Muon1_trkiso", &muon1_trkiso);
    t1->SetBranchAddress("Muon2_trkiso", &muon2_trkiso);

    t1->SetBranchAddress("Muon1_vx", &muon1_vx);
    t1->SetBranchAddress("Muon2_vx", &muon2_vx);
    t1->SetBranchAddress("Muon1_vy", &muon1_vy);
    t1->SetBranchAddress("Muon2_vy", &muon2_vy);
    t1->SetBranchAddress("Muon1_vz", &muon1_vz);
    t1->SetBranchAddress("Muon2_vz", &muon2_vz);

    TH1D* _h1 = new TH1D("_h1"+fname, "_h1"+fname, 80, -4, 4);

    cout << "entry = " << t1->GetEntries() << endl;
    for( int i = 0; i < t1->GetEntries(); i++ ) {
	t1->GetEntry(i);
	if( i == 10000 ) cout << "i = " << i << endl;
	
	const int _npair = npair;
	int _index = -1;
	double _best_val = -99999;
	for( int j = 0; j < npair; j++ ) {
	    if( !hlt_l1mu && !hlt_mu7 ) continue;
	    if( muon1_type[j] != 1 ) continue;
	    if( muon2_type[j] != 1 ) continue;
	    if( sign[j] == 0 ) continue;
	    if( muon1_pt[j] < 7 ) continue;
	    if( muon2_pt[j] < 7 ) continue;
	    if( fabs(muon1_eta[j]) > 2.4 ) continue;
	    if( fabs(muon2_eta[j]) > 2.4 ) continue;
	    if( muon1_trkiso[j] > 3.0 ) continue;
	    if( muon2_trkiso[j] > 3.0 ) continue;

	    if( reco_mass[j] > _best_val ) {
		_best_val = reco_mass[j];
		_index = j;
	    }
	}
	for( int j = 0; j < npair; j++ ) {
	    if( j != _index ) continue;
	    
	    double dA = deltaA(muon1_eta[j],muon1_phi[j],
			       muon2_eta[j],muon2_phi[j]);
	    _h1->Fill(dA);
	}
    }
    _h1->Draw();
    return _h1;
}

