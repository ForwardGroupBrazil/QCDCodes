{
  _file0->cd();
  globalMatchingAnalyser->cd("matchAnalyzer");
  pos_tkAll->Draw("AP") ;
  pos_tkAll->SetLineColor(kBlue);
  pos_tkAll->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  pos_tkAll->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  sta_muon->Draw("P same");
  sta_muon->SetLineColor(kRed);
  region_fixed->Draw("P same") ;
  region_fixed->SetLineColor(kYellow);
  ////region_dynamic->Draw("P same") ;
  //pos_tkCandFixed->SetMarkerStyle(30);
  //pos_tkCandFixed->Draw("P same");
  pos_selectedTkCandFixed->SetMarkerStyle(4);
  pos_selectedTkCandFixed->Draw("P same");
  tk_muon->SetMarkerStyle(25);
  tk_muon->Draw("P same");
  glb_combined->SetLineColor(kGreen);
  glb_combined->Draw("P same");
}
