#include "TGraph.h"
#include "TGraphErrors.h"
#include <vector>

TGraph* fit(const TH2* histo, int fitType = 2){
  vector<const TGraphErrors*> results;

  static int i = 0;
  i++;

  vector<double> width;
  vector<double> mass;

  // Errors
  vector<double> widthError;
  vector<double> massError;

  // Chi2
  vector<double> chi2;

  //X bin centre and width
  vector<double> binCentre;
  vector<double> binWidth;

  // Fit slices projected along Y from bins in X 
  double cont_min = 100;    //Minimum number of entries
  Int_t binx =  histo->GetXaxis()->GetNbins();

  //histo->RebinY(4);
  //histo->RebinX(2);

  for (int i = 1; i <= binx ; i++) {
    TString iString(i);
    TH1 *histoY =  histo->ProjectionY(" ", i, i);
    double cont = histoY->GetEntries();
    //cout << "Entries: " << cont << endl;
    if (cont >= cont_min) {
      float minfit = histoY->GetMean() - histoY->GetRMS();
      float maxfit = histoY->GetMean() + histoY->GetRMS();
      
      TF1 *fitFcn = new TF1(TString("g")+histo->GetName()+iString,"gaus",minfit,maxfit);
      double x1,x2;
      fitFcn->GetRange(x1,x2);
      //cout << "Range: " << x1 << " " << x2 << endl; 

      histoY->Fit(fitFcn,"QR0","",x1,x2);

//      histoY->Fit(fitFcn->GetName(),"RME");
      double *par = fitFcn->GetParameters();
      double *err = fitFcn->GetParErrors();

      // Values
      mass.push_back(par[1]);
      width.push_back(par[2]);

      // Errors
      massError.push_back(err[1]);
      widthError.push_back(err[2]);

      chi2.push_back(fitFcn->GetChisquare());
      
      double xx= histo->GetXaxis()->GetBinCenter(i);
      binCentre.push_back(xx);
      
      double ex = 0; //FIXME: you can use the bin width
      binWidth.push_back(ex); 
      if (fitFcn) delete fitFcn;
      if (histoY) delete histoY;
    }
    else {
      if (histoY) delete histoY;
      continue;
    }
  }

  // Put the fit results in arrays for TGraphErrors
  const int nn= mass.size();
  double *x = new double[nn];
  double *ym = new double[nn];
  double *e = new double[nn];
  double *eym = new double[nn];
  double *yw = new double[nn];
  double *eyw = new double[nn];
  double *yc = new double[nn];
  
  for (int j=0;j<nn;j++){
    x[j]=binCentre[j];
    ym[j]=mass[j];
    eym[j]=massError[j];
    yw[j]=width[j];
    eyw[j]=widthError[j];

    /*
    if(j==nn-1)
      if(eyw[j] > (eyw[j-1]+0.6)) eyw[j]-=0.6;
    if(j==0){
      if(eyw[0] > (eyw[1]+0.6)) eyw[0]-=0.6;}
    else
      if(eyw[j] > (eyw[j-1]+0.6)) eyw[j]-=0.6;
    if(eyw[j] > 1.) eyw[j]=1;
    */

    yc[j]=chi2[j];
    e[j]=binWidth[j];
  }

  //Create TGraphErrors
  TString name = histo->GetName();
  TGraphErrors *grM = new TGraphErrors(nn,x,ym,e,eym);
  grM->SetTitle(name+"_Mean");
  grM->SetName(name+"_Mean");
  TGraphErrors *grW = new TGraphErrors(nn,x,yw,e,eyw);
  grW->SetTitle(name+"_Width");
  grW->SetName(name+"_Width");
  TGraphErrors *grC = new TGraphErrors(nn,x,yc,e,e);
  grC->SetTitle(name+"_chi2");
  grC->SetName(name+"_chi2");

  if(x) delete x;
  if(ym) delete ym;
  if(e) delete e;
  if(eym) delete eym;
  if(yw) delete yw;
  if(eyw) delete eyw;
  if(yc) delete yc;

  if(fitType==1) return grM;
  if(fitType==2) return grW;
  if(fitType==3) return grC;

  return grW;
}

std::pair<double,double>
fit(const TH1* histo, int fitType = 2)
{
  static int i = 0;
  i++;

  vector<double> width;
  vector<double> mass;

  // Errors
  vector<double> widthError;
  vector<double> massError;
  
  // Chi2
  vector<double> chi2;
  
  
  double cont_min = 1; 
  Int_t binx =  histo->GetXaxis()->GetNbins();
  
  //histo->RebinY(4);
  //histo->RebinX(2);
  
  double cont = histo->GetEntries();
  
  if (cont >= cont_min) {
    
    float minfit = histo->GetMean() - histo->GetRMS();
    float maxfit = histo->GetMean() + histo->GetRMS();
    
    
    TString iString(i);   
    
    
    TF1 *fitFcn = new TF1(TString("g")+histo->GetName()+iString,"gaus",minfit,maxfit);
    double x1,x2;
    fitFcn->GetRange(x1,x2);
    histo->Fit(fitFcn,"QR0","",x1,x2);
    double *par = fitFcn->GetParameters();
    double *err = fitFcn->GetParErrors();
    
    // Values
    mass.push_back(par[1]);
    width.push_back(par[2]);
    
    // Errors
    massError.push_back(err[1]);
    widthError.push_back(err[2]);
    
    chi2.push_back(fitFcn->GetChisquare());
    
    
    if (fitFcn) delete fitFcn;
  }
  
  const int nn= mass.size();
  double *x = new double[nn];
  double *ym = new double[nn];
  double *e = new double[nn];
  double *eym = new double[nn];
  double *yw = new double[nn];
  double *eyw = new double[nn];
  double *yc = new double[nn];
  
  for (int j=0;j<nn;j++){
    //x[j]=binCentre[j];
    ym[j]=mass[j];
    eym[j]=massError[j];
    yw[j]=width[j];
    eyw[j]=widthError[j];
    
    yc[j]=chi2[j];
    e[j]=0.;
    
  }
  
  if(nn>1) cout << "ERROR in fit macro: nn > 1 nn=" << nn << endl;
  
  std::pair<double,double> returnVal;
  if(fitType==1) {returnVal.first = ym[0]; returnVal.second = eym[0]; }
  if(fitType==2) {returnVal.first = yw[0]; returnVal.second = eyw[0]; }
  if(fitType==3) {returnVal.first = yc[0]; returnVal.seconsd = e[0]; }
  
  if(x) delete x;
  if(ym) delete ym;
  if(e) delete e;
  if(eym) delete eym;
  if(yw) delete yw;
  if(eyw) delete eyw;
  if(yc) delete yc;
  
  return returnVal;
}
