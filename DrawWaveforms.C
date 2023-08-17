#include <TApplication.h>
#include <TGClient.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TGLayout.h>
#include <TGWindow.h>
#include <TGLabel.h>
#include <TGNumberEntry.h>
#include <TString.h>
#include <iostream>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "ChannelEntry.h"
#include "TGFileDialog.h"
#include <TGMenu.h>
#include <RQ_OBJECT.h>
#include <TRootEmbeddedCanvas.h>
// #include "DECONVOLUTION.h"
#include <TGDockableFrame.h>

enum ETestCommandIdentifiers {
   M_FILE_OPEN,
   M_FILE_EXIT

};
const char *filetypes[] = { "All files",     "*",
                            "ROOT files",    "*.root",
                            "ROOT macros",   "*.C",
                            "Text files",    "*.[tT][xX][tT]",
                            0,               0 };


/////////////
class MyMainFrame : public TGMainFrame {

private:
   static const Int_t total_channels = 128;
   TRootEmbeddedCanvas  *fEcan;
   TGNumberEntry       *fNumber;
   Int_t entrynum;
   TTree *t1;
   ChannelEntry channel_info[total_channels];
   // DECONVOLUTION decon;
   TString fName;
   Int_t n;
   TGNumberEntry       *fChannel;
   Int_t ChNum;

   TGDockableFrame    *fMenuDock;
   TGMenuBar          *fMenuBar;
   TGPopupMenu        *fMenuFile, *fMenuView;
   TGLayoutHints      *fMenuBarLayout, *fMenuBarItemLayout, *fMenuBarHelpLayout;

public:
   MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h, TString s);
   virtual ~MyMainFrame();
   void DoExit();
   void DoDraw();
   //void EventInfo(Int_t event, Int_t px, Int_t py, TObject *selected);
   void SetEntry();
   void SetChannel();
   void ReadFile(TString s);

   void CloseWindow();
   void DoButton();
   void HandleMenu(Int_t id);
   void HandlePopup() { printf("menu popped up\n"); }
   void HandlePopdown() { printf("menu popped down\n"); }



   ClassDef(MyMainFrame, 0)
};


MyMainFrame::MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h, TString s) :
   TGMainFrame(p, w, h)
{
      // decon.SetResponse();

   
///////////////////
   fMenuDock = new TGDockableFrame(this);
   AddFrame(fMenuDock, new TGLayoutHints(kLHintsExpandX, 0, 0, 1, 0));
   fMenuDock->SetWindowName("GuiTest Menu");
   fMenuBarLayout = new TGLayoutHints(kLHintsTop | kLHintsExpandX);
   fMenuBarItemLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 4, 0, 0);
   fMenuFile = new TGPopupMenu(gClient->GetRoot());
   fMenuFile->AddEntry("&Open...", M_FILE_OPEN);
   fMenuFile->AddEntry("E&xit", M_FILE_EXIT);

   fMenuView = new TGPopupMenu(gClient->GetRoot());
   fMenuDock->Connect("Undocked()", "MyMainFrame", this, "HandleMenu(=M_VIEW_UNDOCK)");
   fMenuFile->Connect("Activated(Int_t)", "MyMainFrame", this,
                      "HandleMenu(Int_t)");
   fMenuFile->Connect("PoppedUp()", "MyMainFrame", this, "HandlePopup()");
   fMenuFile->Connect("PoppedDown()", "MyMainFrame", this, "HandlePopdown()");
   fMenuBar = new TGMenuBar(fMenuDock, 1, 1, kHorizontalFrame);
   fMenuBar->AddPopup("&File", fMenuFile, fMenuBarItemLayout);
   fMenuDock->AddFrame(fMenuBar, fMenuBarLayout);
//////////////////
   ChNum = 3500;
   // Create the embedded canvas
   fEcan = new TRootEmbeddedCanvas(0,this,1920,1080);
   Int_t wid = fEcan->GetCanvasWindowId();
   TCanvas *myc = new TCanvas("MyCanvas", 10,10,wid);
   fEcan->AdoptCanvas(myc);
   myc->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)","MyMainFrame",this,
               "EventInfo(Int_t,Int_t,Int_t,TObject*)");

   AddFrame(fEcan, new TGLayoutHints(kLHintsTop | kLHintsLeft |
                                     kLHintsExpandX  | kLHintsExpandY,0,0,1,1));
   TGHorizontalFrame *hframe = new TGHorizontalFrame(this, 200, 40);

////////chnum
   TGLabel* fLchannel = new TGLabel(this, "Channel Number");
   AddFrame(fLchannel,  new TGLayoutHints(kLHintsTop | kLHintsCenterX, 0 , 0 , 0 , 0));
   fChannel = new TGNumberEntry(this, 0, 9,999, TGNumberFormat::kNESInteger,
                                               TGNumberFormat::kNEANonNegative,
                                               TGNumberFormat::kNELLimitMinMax,
                                               0,128);
   fChannel->Connect("ValueSet(Long_t)", "MyMainFrame", this, "SetChannel()");
   AddFrame(fChannel, new TGLayoutHints(kLHintsTop | kLHintsCenterX, 0 , 0 , 0 , 0));
////////////evnum
   TGLabel* fLevent = new TGLabel(this, "Event Number");
   AddFrame(fLevent,  new TGLayoutHints(kLHintsTop | kLHintsCenterX, 0 , 0 , 0 , 0));
   fNumber = new TGNumberEntry(this, 0, 9,999, TGNumberFormat::kNESInteger,
                                               TGNumberFormat::kNEANonNegative,
                                               TGNumberFormat::kNELLimitMinMax,
                                               0,1000000000000);
   fNumber->Connect("ValueSet(Long_t)", "MyMainFrame", this, "SetEntry()");
   (fNumber->GetNumberEntry())->Connect("ReturnPressed()", "MyMainFrame", this,"SetEntry()");

   AddFrame(fNumber, new TGLayoutHints(kLHintsTop | kLHintsCenterX, 0 , 0 , 0 , 0));
   //AddFrame(fNumber, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 5, 5));


   TGTextButton *exit = new TGTextButton(hframe, "&Exit ");
   exit->Connect("Pressed()", "MyMainFrame", this, "DoExit()");
   hframe->AddFrame(exit, new TGLayoutHints(kLHintsRight, 5, 5, 3, 4));
   AddFrame(hframe, new TGLayoutHints(kLHintsRight, 2, 2, 2, 2));
   
   // Set a name to the main frame
   SetWindowName("Embedded Canvas Status Info");
   MapSubwindows();

   // Initialize the layout algorithm via Resize()
   Resize(GetDefaultSize());

   // Map main frame
   MapWindow();
   Print();
}

MyMainFrame::~MyMainFrame()
{
   // Clean up main frame...
   Cleanup();
   delete fEcan;
}

void MyMainFrame::DoDraw()
{
   n = channel_info[ChNum].wf_size;
   Printf("Slot DoDraw()");
   channel_info[ChNum].Set_Zero_Level_Area(300);
   int zl = channel_info[ChNum].Get_Zero_Level();
   int amp = channel_info[ChNum].Get_Amplitude();
   TCanvas *c1 = fEcan->GetCanvas();
   c1->SetFillColor(42);
   c1->SetGrid();
   Double_t x[2048] = {0.}; Double_t y[1024] = {0.}; Double_t z[2048] = {0.};
   cout << amp << endl;



   //for (int cc = 0; cc < 2; cc++) channel_info[ChNum].SplineWf();
   for (Int_t i=0;i<n;i++) {
        x[i] = i;

      z[i] = channel_info[ChNum].wf[i];
   }
   // decon.Reset();
   // decon.SetZl(channel_info[ChNum].Get_Zero_Level(50));
   // decon.SetSignal(channel_info[ChNum].wf);
   // decon.deconvolution();
   // decon.GetSignal(y);
                channel_info[ChNum].SplineWf();
                channel_info[ChNum].SplineWf();
                channel_info[ChNum].SplineWf();
   for (Int_t i=0;i<n;i++) {
      y[i] = channel_info[ChNum].wf[i];
   }


   ///////////
      TGraph *gr = new TGraph(n,x,z);
      gr->SetLineColor(2);
      gr->SetLineWidth(4);
      gr->SetMarkerColor(1);
      gr->SetMarkerStyle(21);
      gr->SetMarkerSize(1);

      gr->SetTitle("ADC Waveform");
      gr->GetXaxis()->SetTitle("Pulse time [ns]");
      gr->GetYaxis()->SetTitle("ADC channels");
      gr->Draw("APL");

      TGraph *gr1 = new TGraph(n,x,y);
      gr1->SetLineColor(4);
      gr1->SetLineWidth(2);
      gr1->SetMarkerColor(1);
      gr1->SetMarkerStyle(21);
      gr1->SetMarkerSize(1);
      gr1->Draw("same");
      c1->Update();
      c1->GetFrame()->SetFillColor(21);
      c1->GetFrame()->SetBorderSize(12);
      c1->Modified();
      c1->Update();      
   // TCanvas::Update() draws the frame, after which it can be changed



}

void MyMainFrame::DoExit()
{
   printf("Exit application...");
   gApplication->Terminate(0);
}

void MyMainFrame::ReadFile(TString s)
{
   TFile *f = new TFile(s);
   t1 = (TTree*)f->Get("adc64_data");
   for (Int_t channel = 0; channel < total_channels; channel++) 
   {
      t1->SetBranchAddress(Form("channel_%i",channel), &channel_info[channel]);
   }
}

void MyMainFrame::CloseWindow()
{
   gApplication->Terminate();
}

void MyMainFrame::HandleMenu(Int_t id)
{
   switch (id) 
   {
      case M_FILE_OPEN:
      {
         static TString dir(".");
         TGFileInfo fi;
         fi.fFileTypes = filetypes;
         fi.SetIniDir(dir);
         printf("fIniDir = %s\n", fi.fIniDir);
         new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
         printf("Open file: %s (dir: %s)\n", fi.fFilename, fi.fIniDir);
         dir = fi.fIniDir;
         fName = fi.fFilename;
         ReadFile(fName);
      }
      break;
   }
}



void MyMainFrame::SetEntry()
{
   entrynum = fNumber->GetNumberEntry()->GetIntNumber();
   t1->GetEntry(entrynum);
   DoDraw();
}
void MyMainFrame::SetChannel()
{
   ChNum = fChannel->GetNumberEntry()->GetIntNumber();
   cout << ChNum << endl;
    DoDraw();
}

void DrawWaveforms(  TString file_path = "/home/doc/SiPM_low_energy_detector/new_adc_data/Jp/diff_k/Co57",
                     TString file_name = "Co57_k_2.root")
{
   new MyMainFrame(gClient->GetRoot(), 200, 200, file_path+"/"+file_name);
}
