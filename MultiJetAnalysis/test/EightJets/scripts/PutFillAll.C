#include "PutTMVA.C"
#include "FillHistos.C"
void PutFillAll()
{
  cout<<"Adding MVA branches to signal"<<endl;
  PutTMVA("flatTree_Col800Sig200_weights","Col800Sig200");
  cout<<"Adding MVA branches to qcd"<<endl;
  PutTMVA("flatTree_qcd_weights","Col800Sig200");
  cout<<"Adding MVA branches to data"<<endl;
  PutTMVA("flatTree_data","Col800Sig200");
  
  cout<<"Filling Histograms"<<endl;
  FillHistos("flatTree_Col800Sig200_weights_tmva.root",true);
  FillHistos("flatTree_qcd_weights_tmva.root",true);
  FillHistos("flatTree_data_tmva.root",false);
}
