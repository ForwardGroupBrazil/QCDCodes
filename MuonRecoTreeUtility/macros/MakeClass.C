/*Root macro to make a header and .C file for basic analysis
  of CMS1 ntuples. Usage:

  [kalavase@stau ~/rootmacros]$ root
  root [0] .L makeCMS1Header.C
  root [1] makeCMS1Header("tablemaker_Zmumu_ntuple.root")

*/




//makeCMS1Files(std::string fname) {
makeCMS1Files( TChain* ev ) {

  using namespace std;
  
  //  TFile *f = TFile::Open( fname.c_str() );
  //  if(f->IsZombie()) { 
  //    cout << "File is not a valid root file, or root file is corruped" << endl;
  //    cout << "Exiting..." << endl;
  //  }
  
  cout << "line 24" << endl;
  ofstream headerf;
  ofstream codef;
  headerf.open("test.h");
  codef.open("test.C");
  headerf << "#include <vector>" << endl;
  headerf << "#ifdef __MAKECINT__" << endl;
  headerf << "//#pragma link C++ class vector<vector<float> >+;" << endl;
  headerf << "//#pragma link C++ class vector<vector<int> >+;" << endl;
  headerf << "#endif" << endl << endl;
  headerf << "#include \"TTree.h\"" << endl;
  headerf << "#include \"map.h\"" << endl;
  headerf << "#include \"utility.h\"" << endl;
  headerf << "#include \"string.h\"" << endl;
   
  cout << "line 36" << endl;
  
  //  TTree *ev = (TTree*)f->Get("MuTrigMC");
  TObjArray *brancharray = ev->GetListOfBranches();  cout << "line 39" << endl;
  TObjArray *leafarray = ev->GetListOfLeaves();  cout << "line 40" << endl;
  
  for(Int_t i = 0 ; i < leafarray->GetSize(); i++) {
    TLeaf *leaf = (TLeaf*)leafarray->At(i);
    string leaftype(leaf->GetTypeName());
    string leafname(leaf->GetName());
     
    if((leaftype=="Int_t" || leaftype=="Float_t") && leafname!="_" ) 
      headerf << "\t" << leaftype << "\t" << leafname << ";" << endl;
  }


  
  for(Int_t i = 0; i< brancharray->GetSize(); i++) {

    //Class name is blank for a int of float
    TBranch *branch = (TBranch*)brancharray->At(i);
    string branchname( branch->GetName() );
    string classname( branch->GetClassName() );
    string classtitle( branch->GetTitle() );
    
    if(classname!="TH2F") {
      if(classname != "") 
	headerf << "\t" << classname << "\t*" << branchname << ";" << endl;
    }
  }

  headerf << "void Init(TTree *tree) {" << endl;
  

  //SetBranchAddresses
  for(Int_t i = 0 ; i < leafarray->GetSize(); i++) {
    TLeaf *leaf = (TLeaf*)leafarray->At(i);
    string leaftype(leaf->GetTypeName());
    string leafname(leaf->GetName());
    
    
    if( (leaftype=="Int_t" || leaftype=="Float_t") && leafname!="_" && !(leafname.find(".first")!=string::npos) && !(leafname.find(".second")!=string::npos) ) {
      //if( (leaftype=="Int_t" || leaftype=="Float_t") && leafname!="_" ) {
      headerf << "\ttree->SetBranchAddress(\"" << leafname
	      << "\", &" << leafname << ");"
	      << endl;
    }
  }

  for(Int_t i = 0; i< brancharray->GetSize(); i++) {
    TBranch *branch = (TBranch*)brancharray->At(i);
    string branchname( branch->GetName() );
    string classname( branch->GetClassName() );
    string classtitle( branch->GetTitle() );
    
    if(classname!="TH2F") {
      if(classname != "") 
	headerf << "\ttree->SetBranchAddress(\"" << branchname
		<< "\", &" << branchname << ");"
		<< endl;
    }
  }

  headerf << "}" << endl;



  //now make the source file
  codef << "#include \"TH1F.h\"" << endl;
  codef << "#include \"TH2F.h\"" << endl;
  codef << "//#include \"Math/LorentzVector.h\"" << endl;
  codef << "#include <vector>" << endl;
  codef << "#include \"test.h\"" << endl << endl;
  codef << "int ScanTree ( TTree* tree) {" << endl << endl;
  codef << "Init(tree);" << endl;
  codef << "  TH1F *l3TestPt = new TH1F("l3TestPt","histogram of L3 p_{T}",100,0,100);" << endl;
  codef << "  TH1F *l3deltaPt = new TH1F("l3deltaPt","histogram of L3 #Delta p_{T}",100,-10.,10.);" << endl;
  codef << endl << "\tint nEntries = tree->GetEntries();" << endl << endl;
  codef << "\t //Event Loop" << endl;
  codef << "\tfor( int i = 0; i < nEntries; i++) {" << endl;
  codef << "\t\ttree->GetEntry(i);" << endl << endl << endl;
  codef << "    for (int iL3 = 0; iL3 < nL3; iL3++) {" << endl;
  codef << "      l3TestPt->Fill((*l3Pt).at(iL3));" << endl;
  codef << "      l3deltaPt->Fill((*l2Pt).at(iL3)-(*tkTrackPt).at(iL3));" << endl;
  codef << "    }" << endl;
  codef << "\t\tcout << \"//N L3: \" << nL3 << endl;" << endl;
  
  codef << " \t}" << endl << "\treturn 0;" << endl << "}" << endl;
  
  
  
  headerf.close();
  codef.close();
  //  f->Close();
}
