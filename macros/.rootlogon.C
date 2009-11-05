// This is the file rootlogon.C
{
   printf("\n Beginning new ROOT session\n");

   TStyle *myStyle  = new TStyle("MyStyle","My Root Styles");

   // from ROOT plain style

   myStyle->SetTitleColor(1);
   myStyle->SetStatColor(0);

   myStyle->SetLabelSize(0.03,"xyz"); // size of axis values

   myStyle->SetPalette(1);

   // default canvas positioning
   myStyle->SetCanvasBorderMode(0);
   myStyle->SetCanvasColor(0);
   myStyle->SetCanvasDefX(900);
   myStyle->SetCanvasDefY(20);
   myStyle->SetCanvasDefH(550);
   myStyle->SetCanvasDefW(540);

   myStyle->SetPadBorderMode(0);
   myStyle->SetPadColor(0);
   myStyle->SetPadGridX(true);
   myStyle->SetPadGridY(true);
   myStyle->SetGridColor(0);
   myStyle->SetGridStyle(3);
   myStyle->SetGridWidth(1);
   myStyle->SetPadTickX(1);
   myStyle->SetPadTickY(1);

// For the frame:
  myStyle->SetFrameBorderMode(0);
  myStyle->SetFrameBorderSize(1);
  myStyle->SetFrameFillColor(0);
  myStyle->SetFrameFillStyle(0);
  myStyle->SetFrameLineColor(1);
  myStyle->SetFrameLineStyle(1);
  myStyle->SetFrameLineWidth(1);

   myStyle->SetPadBottomMargin(0.1);
   myStyle->SetPadTopMargin(0.1);
   myStyle->SetPadLeftMargin(0.1);
   myStyle->SetPadRightMargin(0.1);

   /*
  myStyle->SetHistLineColor(1);
  myStyle->SetHistLineStyle(0);
  myStyle->SetHistLineWidth(1);
 myStyle->SetEndErrorSize(2);
 myStyle->SetErrorX(0.);
  myStyle->SetMarkerStyle(20);
 */


   // US letter
   myStyle->SetPaperSize(20, 24);


   // gROOT->SetStyle("MyStyle"); //uncomment to set this style

   bool foundIt=true;
   // see if CMSSW has been setup
   char *cmsbase=gSystem->Getenv("CMSSW_BASE");
   if (cmsbase==NULL) {
     cout << " CMSSW environment has not been setup -- "
	  << " FWLite libraries will not be loaded\n" << endl;
     foundIt=false;
   } else {
     cout << " CMSSW environment has been setup \n" << endl;

     char *search=gSystem->Getenv("LD_LIBRARY_PATH");
     string cms_path = search;
     
     TString FWLiteLib = "libFWCoreFWLite.so";
     const char* foundlib =gSystem->Which(search, FWLiteLib, 0);
     
     if (! foundlib) {
       FWLiteLib = "libPhysicsToolsFWLite.so";
       foundlib =gSystem->Which(search, FWLiteLib, 0);
       if (! foundlib) {
	 cout << "Could not find any FWLite libraries to load " << endl;       
	 foundIt=false;
       }
     }
   }
   if (foundIt){
     cout << "Loading: " << FWLiteLib << endl;
     gSystem->Load(FWLiteLib);
     AutoLibraryLoader::enable();
     gSystem->Load("libTreeLib.so");
   }
}

