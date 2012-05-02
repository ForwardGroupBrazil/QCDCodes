#include "PutTMVA.C"
#include "FillHistos.C"
void PutFillAll(TString TRAIN)
{
  cout<<"Adding MVA branches to signal"<<endl;
  PutTMVA("data/flatTree_"+TRAIN+"_weights",TRAIN);
  cout<<"Adding MVA branches to qcd"<<endl;
  PutTMVA("data/flatTree_qcd_weights",TRAIN);
  cout<<"Adding MVA branches to data"<<endl;
  PutTMVA("data/flatTree_data",TRAIN);
  cout<<"Filling Histograms"<<endl;
  FillHistos("flatTree_"+TRAIN+"_weights_tmva"+TRAIN+".root",true);
  FillHistos("flatTree_qcd_weights_tmva"+TRAIN+".root",true);
  FillHistos("flatTree_data_tmva"+TRAIN+".root",false);
}
