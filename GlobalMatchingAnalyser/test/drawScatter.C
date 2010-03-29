{
  _file0->cd();
  globalMatchingAnalyser->cd("matchAnalyzer");
  pos_tkAll->Draw("AP") ;
  pos_tkAll->GetHistogram()->GetXaxis()->SetRangeUser(-5.,5.);
  pos_tkAll->GetHistogram()->GetYaxis()->SetRangeUser(-5.,5.);
  sta_muon->Draw("P same");
  region_fixed->Draw("P same") ;
  //region_dynamic->Draw("P same") ;
  pos_tkCandFixed->Draw("P same");
  pos_selectedTkCandFixed->Draw("P same");
  tk_muon->Draw("P same"); 

}
