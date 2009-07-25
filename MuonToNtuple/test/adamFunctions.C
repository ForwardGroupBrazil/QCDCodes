double deltaPhi(double phi1, double phi2) { 
  double result = phi1 - phi2;
  while (result > TMath::Pi()) result -= 2*TMath::Pi();
  while (result <= -TMath::Pi) result += 2*TMath::Pi();
  return result;
}

double deltaPhi(float phi1, double phi2) {
  return deltaPhi(static_cast<double>(phi1), phi2);
}

double deltaPhi(double phi1, float phi2) {
  return deltaPhi(phi1, static_cast<double>(phi2));
}

double deltaPhi(float phi1, float phi2) {
  return deltaPhi(static_cast<double>(phi1),
		  static_cast<double>(phi2));
} 

double deltaR2(double eta1, double phi1, double eta2, double phi2) {
  double deta = eta1 - eta2;
  double dphi = deltaPhi(phi1, phi2);
  return deta*deta + dphi*dphi;
}

double deltaR(double eta1, double phi1, double eta2, double phi2) {
  return sqrt(deltaR2 (eta1, phi1, eta2, phi2));
}

double deltaA(double eta1, double phi1, double eta2, double phi2) {
  return ((eta1-eta2)+(deltaPhi(phi1,phi2)) - 2*TMath::Pi());
}

// Print all canvases in separate EPS files
void printCanvasesType(TString type=".eps"){
  TIter iter(gROOT->GetListOfCanvases());
  TCanvas *c;
  while( (c = (TCanvas *)iter()) ) {
    c->cd();
    c->Print(0,type);
  }
}
