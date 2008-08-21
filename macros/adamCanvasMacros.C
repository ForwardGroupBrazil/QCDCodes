/*
 * TCanvas(const char* name, const char* title = "", Int_t form = 1)
 * Create a new canvas with a predefined size form.
 * If form < 0  the menubar is not shown.
 *
 *  form = 1    700x500 at 10,10 (set by TStyle::SetCanvasDefH,W,X,Y)
 *  form = 2    500x500 at 20,20
 *  form = 3    500x500 at 30,30
 *  form = 4    500x500 at 40,40
 *  form = 5    500x500 at 50,50
 */

/*
 * TCanvas(const char *name, const char *title, Int_t ww, Int_t wh)
 *
 *  Create a new canvas at a random position.
 *
 *  ww is the canvas size in pixels along X
 *      (if ww < 0  the menubar is not shown)
 *  wh is the canvas size in pixels along Y
 */

//1
TCanvas * newCanvas(TString name="", TString title="",
                     Int_t xdiv=0, Int_t ydiv=0, Int_t form = 1, Int_t w=-1){
  static int i = 1;
  if (name == "") {
    name = TString("Canvas ") + i;
    i++;
  }
  if (title == "") title = name;
  if (w<0) {
    TCanvas * c = new TCanvas(name,title, form);
  } else {
    TCanvas * c = new TCanvas(name,title,form,w);
  }
  if (xdiv*ydiv!=0) c->Divide(xdiv,ydiv);
  c->cd(1);
  return c;
}

// Print all canvases in a single PS file
void printCanvasesPS(TString name){
  TPostScript * ps = new TPostScript(name,112);
  TIter iter(gROOT->GetListOfCanvases());
  TCanvas *c;
  while( (c = (TCanvas *)iter()) )
    {
      c->cd();
      cout << "Printing " << c->GetName() << endl;
      ps->NewPage();
      c->Draw();
    }
  cout << " File " << name << " was created" << endl;
  ps->Close();
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
