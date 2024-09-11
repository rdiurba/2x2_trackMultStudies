	
#include "TH1.h"
#include "TGraph.h"
#include "TH2.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TVectorD.h"
#include "TTimeStamp.h"
#include <fstream>
#include "TMinuit.h"
#include "TString.h"
#include <vector>
#include <string.h>
#include "TLatex.h"
#include "TPaveStats.h"
#include "TDatime.h"
#include "TColor.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TPrincipal.h"
#include "TDecompChol.h"
#include "TEfficiency.h"
#include <iomanip>

void plotMinervaValidation()
{


//std::cout<<std::setprecision(2);
gROOT->LoadMacro("protoDUNEStyle.C");
gROOT->SetStyle("protoDUNEStyle");
gROOT->ForceStyle();
gStyle->SetTitleX(0.35);
gStyle->SetOptFit(111);
//gStyle->SetPadRightMargin(0.15);
//gStyle->SetPadLeftMargin(0.15);
TCanvas c1=TCanvas();

TLatex tL;
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
TFile fSim("minervaPlotterSim.root");
TH1D* startXMC=(TH1D*)fSim.Get("startX");
TH1D* startYMC=(TH1D*)fSim.Get("startY");
TH1D* startZMC=(TH1D*)fSim.Get("startZ");
TH1D* endXMC=(TH1D*)fSim.Get("endX");
TH1D* endYMC=(TH1D*)fSim.Get("endY");
TH1D* endZMC=(TH1D*)fSim.Get("endZ");
TH1D* lengthMC=(TH1D*)fSim.Get("length");
TH1D* dotProductMC=(TH1D*)fSim.Get("dotProduct");
TH1D* kalmanQualMC=(TH1D*)fSim.Get("kalmanQual");
startXMC->SetName("startXMC");
startYMC->SetName("startYMC");
startZMC->SetName("startZMC");
endXMC->SetName("endXMC");
endYMC->SetName("endYMC");
endZMC->SetName("endZMC");
lengthMC->SetName("lengthMC");
dotProductMC->SetName("dotProductMC");
kalmanQualMC->SetName("kalmanQual");





TFile fMinerva("minervaPlotterData.root");
TH1D* startXData=(TH1D*)fMinerva.Get("startX");
TH1D* startYData=(TH1D*)fMinerva.Get("startY");
TH1D* startZData=(TH1D*)fMinerva.Get("startZ");
TH1D* endXData=(TH1D*)fMinerva.Get("endX");
TH1D* endYData=(TH1D*)fMinerva.Get("endY");
TH1D* endZData=(TH1D*)fMinerva.Get("endZ");
TH1D* lengthData=(TH1D*)fMinerva.Get("length");
TH1D* dotProductData=(TH1D*)fMinerva.Get("dotProduct");
TH1D* kalmanQualData=(TH1D*)fMinerva.Get("kalmanQual");

startXMC->Scale(0.3);
startYMC->Scale(0.3);
startZMC->Scale(0.3);
dotProductMC->Scale(0.3);
kalmanQualMC->Scale(0.3);
lengthMC->Scale(0.3);

/*
startXData->Scale(startXMC->GetEntries()/startXData->GetEntries());
startYData->Scale(startXMC->GetEntries()/startXData->GetEntries());
startZData->Scale(startXMC->GetEntries()/startXData->GetEntries());
endXData->Scale(startXMC->GetEntries()/startXData->GetEntries());
endYData->Scale(startXMC->GetEntries()/startXData->GetEntries());
endZData->Scale(startXMC->GetEntries()/startXData->GetEntries());
dotProductData->Scale(startXMC->GetEntries()/startXData->GetEntries());
kalmanQualData->Scale(startXMC->GetEntries()/startXData->GetEntries());
*/
dotProductMC->GetXaxis()->CenterTitle();
dotProductMC->GetYaxis()->CenterTitle();
dotProductMC->SetLineColor(kRed);
dotProductMC->SetMarkerColor(kRed);


  TLegend *lVertex = new TLegend(0.3,0.6,0.5,0.85);
  lVertex->SetTextFont(133);
  lVertex->SetTextSize(25);
  lVertex->SetHeader("Scaled to 3E18 POT");
  lVertex->AddEntry(dotProductData,"Data","lp");
  lVertex->AddEntry(dotProductMC,"Simulation","l");


dotProductMC->GetYaxis()->SetRangeUser(0,dotProductMC->GetMaximum()*1.8);
dotProductMC->GetXaxis()->SetTitle("cos(#theta)");
dotProductMC->GetYaxis()->SetTitle("Number of Tracks");
dotProductMC->SetTitle("US MINERvA Tracks");
dotProductMC->Draw("HIST"); dotProductData->Draw("E0 P SAME");      tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}"); lVertex->Draw("SAME");
c1.Print("dotProductData.png"); c1.Print("dotProductData.pdf");

startXMC->GetYaxis()->SetRangeUser(0,startXMC->GetMaximum()*1.8);
startXMC->GetXaxis()->CenterTitle();
startXMC->GetYaxis()->CenterTitle();
startXMC->SetLineColor(kRed);
startXMC->SetMarkerColor(kRed);
startXMC->GetXaxis()->SetTitle("Starting Poisiton X [cm]");
startXMC->GetYaxis()->SetTitle("Number of Tracks");
startXMC->SetTitle("US MINERvA Tracks");
startXMC->Draw("HIST"); startXData->Draw("E0 P SAME");      tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}"); lVertex->Draw("SAME");
c1.Print("startXData.png"); c1.Print("startXData.pdf");


startZMC->GetYaxis()->SetRangeUser(0,startZMC->GetMaximum()*1.8);
startZMC->GetXaxis()->CenterTitle();
startZMC->GetYaxis()->CenterTitle();
startZMC->SetLineColor(kRed);
startZMC->SetMarkerColor(kRed);
startZMC->GetXaxis()->SetTitle("Starting Position Z [cm]");
startZMC->GetYaxis()->SetTitle("Number of Tracks");
startZMC->SetTitle("US MINERvA Tracks");
startZMC->Draw("HIST"); startZData->Draw("E0 P SAME");      tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}"); lVertex->Draw("SAME");
c1.Print("startZData.png"); c1.Print("startZData.pdf");

startYMC->GetYaxis()->SetRangeUser(0,startYMC->GetMaximum()*1.8);
startYMC->GetXaxis()->CenterTitle();
startYMC->GetYaxis()->CenterTitle();
startYMC->SetLineColor(kRed);
startYMC->SetMarkerColor(kRed);
startYMC->GetXaxis()->SetTitle("Starting Position Y [cm]");
startYMC->GetYaxis()->SetTitle("Number of Tracks");
startYMC->SetTitle("US MINERvA Tracks");
startYMC->Draw("HIST"); startYData->Draw("E0 P SAME");      tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}"); lVertex->Draw("SAME");
c1.Print("startYData.png"); c1.Print("startYData.pdf");

lengthMC->GetYaxis()->SetRangeUser(0,lengthMC->GetMaximum()*1.8);
lengthMC->GetXaxis()->CenterTitle();
lengthMC->GetYaxis()->CenterTitle();
lengthMC->SetLineColor(kRed);
lengthMC->SetMarkerColor(kRed);
lengthMC->GetXaxis()->SetTitle("Track Length [cm]");
lengthMC->GetYaxis()->SetTitle("Number of Tracks");
lengthMC->SetTitle("US MINERvA Tracks");
lengthMC->Draw("HIST"); lengthData->Draw("E0 P SAME");      tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}"); lVertex->Draw("SAME");
c1.Print("lengthData.png"); c1.Print("lengthData.pdf");



}

