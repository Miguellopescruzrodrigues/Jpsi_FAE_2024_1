#include <stdio.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <iostream>
#include <math.h>
#include <TLatex.h>
#include <TLegend.h>
using namespace std;

 void Juntar() {

//string quantity = "Pt"; string title = "P_{t}";
//string quantity = "Eta";  string title = "#eta";
string quantity = "Phi"; string title = "#Phi";

string trackerFile_DATA = "results/efficiencies/efficiency/Jpsi_MC_2020/" + quantity + "_" + "trackerMuon.root";
string trackerFile_MC = "results/efficiencies/efficiency/Jpsi_Run_2011/" + quantity + "_" + "trackerMuon.root";
//string standaloneFile = "results/efficiencies/Z0_Run_2011/" + quantity + "_" + "standaloneMuon.root";
//string globalFile = "results/efficiencies/Z0_Run_2011/" + quantity + "_" + "globalMuon.root";

string Titulocanvas = "Efficiency of Tracker Muon Probe";

    TFile *f1 =new TFile(trackerFile_DATA.c_str());
    TFile *f2 =new TFile(trackerFile_MC.c_str());


    //TFile *f2 =new TFile(globalFile.c_str());
    //TFile *f3 =new TFile(standaloneFile.c_str());


  	// TCanvas *c0 = (TCanvas*)f1->Get("trackerMuon_Pt_Efficiency;1");
  	// TCanvas *c1 = (TCanvas*)f1->Get("globalMuon_Pt_Efficiency;1");

TCanvas* c2 = new TCanvas();
string tracker_DATA_EF = "trackerMuon_" + quantity + "_Efficiency;1";
string tracker_MC_EF = "trackerMuon_" + quantity + "_Efficiency;1";
//string standaloneGet = "standaloneMuon_" + quantity + "_Efficiency;1";
//string globalGet = "globalMuon_" + quantity + "_Efficiency;1";

  	TEfficiency* tracker_DATA = (TEfficiency*)f1->Get(tracker_DATA_EF.c_str());
  	 tracker_DATA->SetLineColor(kGreen - 2);
  	 tracker_DATA->SetMarkerColor(kGreen - 2);
     tracker_DATA->SetLineWidth(1);
     tracker_DATA->SetTitle(Titulocanvas.c_str());

    TEfficiency* tracker_MC = (TEfficiency*)f2->Get(tracker_MC_EF.c_str());
     tracker_MC->SetLineColor(kRed);
  	 tracker_MC->SetMarkerColor(kRed);
     tracker_MC->SetLineWidth(1);


    // TEfficiency* standalone = (TEfficiency*)f3->Get(standaloneGet.c_str());
    //  standalone->SetLineColor(kGreen);
  	//  standalone->SetMarkerColor(kGreen);
    //  standalone->SetLineWidth(2);

TH1F *green = new TH1F("h2","Ex",1,-10,10);
TH1F *red = new TH1F("h3","Ex",1,-10,10);
TH1F *blue = new TH1F("h5","Ex",1,-10,10);


green->SetLineColor(kGreen - 2);
red->SetLineColor(kRed);
blue->SetLineColor(kBlue);

TLegend *leg = new TLegend(0.3,0.6,0.1,0.9);
leg->AddEntry(tracker_DATA,"J/#psi real data","pl");
leg->AddEntry(tracker_MC,"Simulated data","pl");
//leg->AddEntry(standalone,"Standalone Muon","pl");
//leg->AddEntry(global,"Global Muon","pl");


 tracker_DATA->Draw();

gPad->Update(); 
auto graph = tracker_DATA->GetPaintedGraph(); 
graph->SetMinimum(0.7);
graph->SetMaximum(1.1); 
gPad->Update(); 

 //trackerFile_MC->Draw("same");
 tracker_MC->Draw("same");

     //CMS Open Data
    TLatex* txCOD = new TLatex();
    txCOD->SetTextSize(0.04);
    txCOD->SetTextAlign(12);
    txCOD->SetTextFont(42);
    txCOD->SetNDC(kTRUE);
    txCOD->DrawLatex(0.2, 0.8, "#bf{CMS Open Data}");


leg->Draw();
  }