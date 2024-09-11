
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

void plotSimpleData()
{



gROOT->LoadMacro("protoDUNEStyle.C");
gROOT->SetStyle("protoDUNEStyle");
gROOT->ForceStyle();
gStyle->SetTitleX(0.35);
gStyle->SetOptFit(111);


TCanvas c1=TCanvas();

TLatex tL;
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");

TFile fPurity("testPurity.root");
TFile fData("testData.root");

TFile fMx2("testMinerva.root");
TFile fMx2Data("testMinervaData.root");

TH1D* dotProduct=(TH1D*)fMx2.Get("dotProductUS");
dotProduct->SetName("dotProductUSMC");
TH1D* dotProductData=(TH1D*)fMx2Data.Get("dotProductUS");
TH1D* deltaXUSFront=(TH1D*)fMx2.Get("deltaXUSFront");
deltaXUSFront->SetName("deltaXUSFrontMC");
TH1D* deltaXUSFrontData=(TH1D*)fMx2Data.Get("deltaXUSFront");


TH1D* deltaYUSFront=(TH1D*)fMx2.Get("deltaYUSFront");
deltaYUSFront->SetName("deltaYUSFrontMC");
TH1D* deltaYUSFrontData=(TH1D*)fMx2Data.Get("deltaYUSFront");


TH1D* deltaXUS=(TH1D*)fMx2.Get("deltaXUS");
deltaXUS->SetName("deltaXUSMC");
TH1D* deltaXUSData=(TH1D*)fMx2Data.Get("deltaXUS");


TH1D* deltaYUS=(TH1D*)fMx2.Get("deltaYUS");
deltaYUS->SetName("deltaYUSMC");
TH1D* deltaYUSData=(TH1D*)fMx2Data.Get("deltaYUS");




TH1D* track_mult=(TH1D*)fPurity.Get("track_multShape");
track_mult->SetName("track_multSimShape");
TH1D* trackData=(TH1D*)fData.Get("track_multShape");
track_mult->GetXaxis()->CenterTitle();
track_mult->GetXaxis()->SetTitle("Number of Reconstructed Tracks");
track_mult->GetYaxis()->SetTitle("Fraction of Interactions");
track_mult->GetYaxis()->SetRangeUser(0.0,0.8);
track_mult->GetYaxis()->CenterTitle();
track_mult->SetTitle("Shape Only Comparisons");


track_mult->SetLineColor(kRed);
track_mult->Draw("HIST");
trackData->Draw("EO P SAME");

    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");


c1.Print("shapeData.png"); c1.Print("shapeData.pdf");



dotProduct->GetXaxis()->CenterTitle();
dotProduct->GetXaxis()->SetTitle("cos(#Delta #Theta)");
dotProduct->GetYaxis()->SetTitle("Fraction of Through-Going Mx2 Tracks");
dotProduct->GetYaxis()->SetRangeUser(0.0,0.8);
dotProduct->GetYaxis()->CenterTitle();
dotProduct->SetTitle("Matched Tracks: Shape Only");
dotProduct->Scale(1.f/dotProduct->GetEntries());
dotProductData->Scale(1.f/dotProductData->GetEntries());

dotProduct->SetLineColor(kRed);
dotProduct->Draw("HIST");
dotProductData->Draw("EO P SAME");

    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");


c1.Print("dotProductDataMx2Through.png"); c1.Print("dotProductDataMx2Through.pdf");



deltaXUS->GetXaxis()->CenterTitle();
deltaXUS->GetXaxis()->SetTitle("#Delta X_{DS}");
deltaXUS->GetYaxis()->SetTitle("Fraction of Through-Going Mx2 Tracks");
deltaXUS->GetYaxis()->SetRangeUser(0.0,0.8);
deltaXUS->GetYaxis()->CenterTitle();
deltaXUS->SetTitle("Matched Tracks: Shape Only");
deltaXUS->Scale(1.f/deltaXUS->GetEntries());
deltaXUSData->Scale(1.f/deltaXUSData->GetEntries());
deltaXUS->GetYaxis()->SetRangeUser(0.0,0.35);
deltaXUS->SetLineColor(kRed);
deltaXUS->Draw("HIST");
deltaXUSData->Draw("EO P SAME");

    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");


c1.Print("deltaXUSDataMx2Through.png"); c1.Print("deltaXUSDataMx2Through.pdf");


deltaYUS->GetXaxis()->CenterTitle();
deltaYUS->GetXaxis()->SetTitle("#Delta Y_{DS}");
deltaYUS->GetYaxis()->SetTitle("Fraction of Through-Going Mx2 Tracks");
deltaYUS->GetYaxis()->SetRangeUser(0.0,0.8);
deltaYUS->GetYaxis()->CenterTitle();
deltaYUS->SetTitle("Matched Tracks: Shape Only");
deltaYUS->Scale(1.f/deltaYUS->GetEntries());
deltaYUSData->Scale(1.f/deltaYUSData->GetEntries());

deltaYUS->SetLineColor(kRed);
deltaYUS->GetYaxis()->SetRangeUser(0.0,0.35);
deltaYUS->Draw("HIST");
deltaYUSData->Draw("EO P SAME");

    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");


c1.Print("deltaYUSDataMx2Through.png"); c1.Print("deltaYUSDataMx2Through.pdf");



deltaXUSFront->GetXaxis()->CenterTitle();
deltaXUSFront->GetXaxis()->SetTitle("#Delta X_{US}");
deltaXUSFront->GetYaxis()->SetTitle("Fraction of Through-Going Mx2 Tracks");
deltaXUSFront->GetYaxis()->SetRangeUser(0.0,0.8);
deltaXUSFront->GetYaxis()->CenterTitle();
deltaXUSFront->SetTitle("Matched Tracks: Shape Only");
deltaXUSFront->Scale(1.f/deltaXUSFront->GetEntries());
deltaXUSFrontData->Scale(1.f/deltaXUSFrontData->GetEntries());
deltaXUSFront->SetLineColor(kRed);
deltaXUSFront->GetYaxis()->SetRangeUser(0.0,0.35);
deltaXUSFront->Draw("HIST");
deltaXUSFrontData->Draw("EO P SAME");

    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");


c1.Print("deltaXUSFrontDataMx2Through.png"); c1.Print("deltaXUSFrontDataMx2Through.pdf");


deltaYUSFront->GetXaxis()->CenterTitle();
deltaYUSFront->GetXaxis()->SetTitle("#Delta Y_{US}");
deltaYUSFront->GetYaxis()->SetTitle("Fraction of Through-Going Mx2 Tracks");
deltaYUSFront->GetYaxis()->SetRangeUser(0.0,0.8);
deltaYUSFront->GetYaxis()->CenterTitle();
deltaYUSFront->SetTitle("Matched Tracks: Shape Only");
deltaYUSFront->Scale(1.f/deltaYUSFront->GetEntries());
deltaYUSFrontData->Scale(1.f/deltaYUSFrontData->GetEntries());
deltaYUSFront->GetYaxis()->SetRangeUser(0.0,0.35);
deltaYUSFront->SetLineColor(kRed);

deltaYUSFront->Draw("HIST");
deltaYUSFrontData->Draw("EO P SAME");

    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");


c1.Print("deltaYUSFrontDataMx2Through.png"); c1.Print("deltaYUSFrontDataMx2Through.pdf");



std::cout<<deltaXUSFront->GetName()<<","<<deltaXUSFront->GetMean()<<","<<deltaXUSFront->GetStdDev()<<","<<deltaXUSFrontData->GetMean()<<","<<deltaXUSFrontData->GetStdDev()<<std::endl;
std::cout<<deltaXUS->GetName()<<","<<deltaXUS->GetMean()<<","<<deltaXUS->GetStdDev()<<","<<deltaXUSData->GetMean()<<","<<deltaXUSData->GetStdDev()<<std::endl;
std::cout<<deltaYUSFront->GetName()<<","<<deltaYUSFront->GetMean()<<","<<deltaYUSFront->GetStdDev()<<","<<deltaYUSFrontData->GetMean()<<","<<deltaYUSFrontData->GetStdDev()<<std::endl;
std::cout<<deltaYUS->GetName()<<","<<deltaYUS->GetMean()<<","<<deltaYUS->GetStdDev()<<","<<deltaYUSData->GetMean()<<","<<deltaYUSData->GetStdDev()<<std::endl;





}

