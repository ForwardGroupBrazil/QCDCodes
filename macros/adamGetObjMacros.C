TDirectory * 
getDirectory( TString & fileName_ , TString & directoryName_ )
{
  TFile * file = new TFile(fileName_);
  return file->GetDirectory(directoryName_,true);
}

TDirectory * 
getDirectory( TFile * file_ , TString & directoryName_ )
{
  TDirectory * dir =  file_->GetDirectory(directoryName_);
  return dir;
}

TH1 *
getHistogram(TDirectory * dir_, TString & histoName_)
{
  TH1 * theHisto = (dir_) ? dir_->Get(histoName_) : 0;
  return theHisto;
}
