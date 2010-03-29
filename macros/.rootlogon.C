// This is the file rootlogon.C
{
   printf("\n Beginning new ROOT session\n");

   int macroAvailable = gROOT->LoadMacro("adamStyles.C");
   printf("\n macroAvailable = %d \n",macroAvailable);
   if(macroAvailable==0){
     setTDRStyle();
     gROOT->SetStyle("tdrStyle");
     gROOT->ForceStyle();
   }

   // US letter
   //gStyle->SetPaperSize(20, 24);

   bool foundIt=true;
   // see if CMSSW has been setup
   char *cmsbase=gSystem->Getenv("CMSSW_BASE");
   if (cmsbase==NULL) {
     printf("\n CMSSW environment has not been setup -- \n FWLite libraries will not be loaded\n");
     foundIt=false;
   } else {
     printf(" CMSSW environment has been setup \n");

     char *search=gSystem->Getenv("LD_LIBRARY_PATH");
     string cms_path = search;
     
     TString FWLiteLib = "libFWCoreFWLite.so";
     const char* foundlib =gSystem->Which(search, FWLiteLib, 0);
     
     if (! foundlib) {
       FWLiteLib = "libPhysicsToolsFWLite.so";
       foundlib =gSystem->Which(search, FWLiteLib, 0);
       if (! foundlib) {
	 printf("Could not find any FWLite libraries to load \n");
	 foundIt=false;
       }
     }
   }
   if (foundIt){
     printf(" Loading: %s \n",FWLiteLib.Data());
     gSystem->Load(FWLiteLib);
     AutoLibraryLoader::enable();
     TString treeLib = "libTreeLib.so";
     const char* foundtreelib =gSystem->Which(search, treeLib, 0);
     if(foundtreelib) {
       printf(" Loading: %s \n",treeLib.Data());
       gSystem->Load(treeLib);
     } else {
       printf("\n No TreeLib found \n");
     }
     gSystem->Load("libSmatrix");
     printf("\n All loaded . . . \n");
   }
}

