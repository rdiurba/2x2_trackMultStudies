
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

TFile fPurity("testSPINE.root");
TFile fData("testSPINEData.root");

TH1D* totalPOTMC=(TH1D*)fPurity.Get("totalPOT");
totalPOTMC->SetName("totalPOTMC");


TH1D* totalPOTData=(TH1D*)fData.Get("totalPOT");

TFile fMx2("testMinervaSPINE.root");
TFile fMx2Data("testMinervaSPINEData.root");
TH2D* vertex2DData=(TH2D*)fData.Get("recoVertex2D");
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




TH1D* track_mult=(TH1D*)fPurity.Get("track_mult");
track_mult->SetName("track_multSimShape");
TH1D* trackData=(TH1D*)fData.Get("track_mult");


      TLegend *lVertex = new TLegend(0.3,0.6,0.5,0.85);
  lVertex->SetTextFont(133);
  lVertex->SetTextSize(25);
  lVertex->SetHeader(Form("%1.2fE17 POT",totalPOTData->GetBinContent(1)/1E4));
  lVertex->AddEntry(trackData,"Data","lp");
  lVertex->AddEntry(track_mult,"Simulation","l");

track_mult->Scale(totalPOTData->Integral()/totalPOTMC->Integral());
track_mult->GetXaxis()->CenterTitle();
track_mult->GetXaxis()->SetTitle("Number of Reconstructed Tracks");
track_mult->GetYaxis()->SetTitle("Number of Selected Interactions");
track_mult->GetYaxis()->SetRangeUser(0.0,track_mult->GetMaximum()*2.0);
track_mult->GetYaxis()->CenterTitle();

track_mult->SetTitle("SPINE Tracks in #nu_{#mu} CCINC");

track_mult->SetLineColor(kRed);
track_mult->Draw("HIST");
trackData->Draw("EO P SAME");

    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}"); lVertex->Draw("SAME");


c1.Print("shapeDataSPINE.png"); c1.Print("shapeDataSPINE.pdf");



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


c1.Print("dotProductDataMx2ThroughSPINE.png"); c1.Print("dotProductDataMx2ThroughSPINE.pdf");



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


c1.Print("deltaXUSDataMx2ThroughSPINE.png"); c1.Print("deltaXUSDataMx2ThroughSPINE.pdf");


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


c1.Print("deltaYUSDataMx2ThroughSPINE.png"); c1.Print("deltaYUSDataMx2ThroughSPINE.pdf");



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


c1.Print("deltaXUSFrontDataMx2ThroughSPINE.png"); c1.Print("deltaXUSFrontDataMx2ThroughSPINE.pdf");


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


c1.Print("deltaYUSFrontDataMx2ThroughSPINE.png"); c1.Print("deltaYUSFrontDataMx2ThroughSPINE.pdf");


    
vertex2DData->GetXaxis()->CenterTitle();
vertex2DData->GetXaxis()->SetTitle("Vertex Position in X [cm]");
vertex2DData->GetYaxis()->SetTitle("Vertex Position in Z [cm]");
vertex2DData->GetYaxis()->CenterTitle();
vertex2DData->SetTitle("Vertex Distribution");

vertex2DData->Draw("COLZ");


    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");


c1.Print("recoVertex2DDataSPINE.png"); c1.Print("recoVertex2DDataSPINE.pdf");



std::cout<<deltaXUSFront->GetName()<<","<<deltaXUSFront->GetMean()<<","<<deltaXUSFront->GetStdDev()<<","<<deltaXUSFrontData->GetMean()<<","<<deltaXUSFrontData->GetStdDev()<<std::endl;
std::cout<<deltaXUS->GetName()<<","<<deltaXUS->GetMean()<<","<<deltaXUS->GetStdDev()<<","<<deltaXUSData->GetMean()<<","<<deltaXUSData->GetStdDev()<<std::endl;
std::cout<<deltaYUSFront->GetName()<<","<<deltaYUSFront->GetMean()<<","<<deltaYUSFront->GetStdDev()<<","<<deltaYUSFrontData->GetMean()<<","<<deltaYUSFrontData->GetStdDev()<<std::endl;
std::cout<<deltaYUS->GetName()<<","<<deltaYUS->GetMean()<<","<<deltaYUS->GetStdDev()<<","<<deltaYUSData->GetMean()<<","<<deltaYUSData->GetStdDev()<<std::endl;

TH1D* dotProductMC=(TH1D*)fMx2.Get("dotProductUS");
dotProductMC->SetName("dotProductUSMC");



TH1D* startXMC=(TH1D*)fMx2.Get("trackStartXUS");
startXMC->SetName("startXMC");
TH1D* startYMC=(TH1D*)fMx2.Get("trackStartYUS");
startYMC->SetName("startYMC");
TH1D* startZMC=(TH1D*)fMx2.Get("trackStartZUS");
startZMC->SetName("startZMC");


TH1D* endXMC=(TH1D*)fMx2.Get("trackEndXUS");
endXMC->SetName("endXMC");
TH1D* endYMC=(TH1D*)fMx2.Get("trackEndYUS");
endYMC->SetName("endYMC");
TH1D* endZMC=(TH1D*)fMx2.Get("trackEndZUS");
endZMC->SetName("endZMC");

TH1D* endXData=(TH1D*)fMx2Data.Get("trackEndXUS");
TH1D* endYData=(TH1D*)fMx2Data.Get("trackEndYUS");
TH1D* endZData=(TH1D*)fMx2Data.Get("trackEndZUS");


TH1D* startXData=(TH1D*)fMx2Data.Get("trackStartXUS");
TH1D* startYData=(TH1D*)fMx2Data.Get("trackStartYUS");
TH1D* startZData=(TH1D*)fMx2Data.Get("trackStartZUS");
    
TH2D* startData=(TH2D*)fMx2Data.Get("trackStartUS");
TH2D* endData=(TH2D*)fMx2Data.Get("trackEndUS");


dotProductMC->GetXaxis()->CenterTitle();
dotProductMC->GetYaxis()->CenterTitle();
dotProductMC->SetLineColor(kRed);
dotProductMC->SetMarkerColor(kRed);



dotProductMC->Scale(totalPOTData->Integral()/totalPOTMC->Integral());
dotProductMC->GetYaxis()->SetRangeUser(0,dotProductMC->GetMaximum()*2);
dotProductMC->GetXaxis()->SetTitle("cos(#theta)");
dotProductMC->GetYaxis()->SetTitle("Number of Tracks");
dotProductMC->SetTitle("SPINE-Mx2 Through-Going");
dotProductMC->Draw("HIST"); dotProductData->Draw("E0 P SAME");      tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}"); lVertex->Draw("SAME");
c1.Print("dotProductDataSPINE.png"); c1.Print("dotProductDataSPINE.pdf");


startXMC->Scale(totalPOTData->Integral()/totalPOTMC->Integral());
startXMC->GetYaxis()->SetRangeUser(0,startXMC->GetMaximum()*2);
startXMC->GetXaxis()->CenterTitle();
startXMC->GetYaxis()->CenterTitle();
startXMC->SetLineColor(kRed);
startXMC->SetMarkerColor(kRed);
startXMC->GetXaxis()->SetTitle("Starting Poisiton X [cm]");
startXMC->GetYaxis()->SetTitle("Number of Tracks");
startXMC->SetTitle("SPINE-Mx2 Through-Going");
startXMC->Draw("HIST"); startXData->Draw("E0 P SAME");      tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");  lVertex->Draw("SAME");
c1.Print("startXDataSPINE.png"); c1.Print("startXDataSPINE.pdf");


startZMC->Scale(totalPOTData->Integral()/totalPOTMC->Integral());
startZMC->GetYaxis()->SetRangeUser(0,startZMC->GetMaximum()*2);

startZMC->GetXaxis()->CenterTitle();
startZMC->GetYaxis()->CenterTitle();
startZMC->SetLineColor(kRed);
startZMC->SetMarkerColor(kRed);
startZMC->GetXaxis()->SetTitle("Starting Position Z [cm]");
startZMC->GetYaxis()->SetTitle("Number of Tracks");
startZMC->SetTitle("SPINE-Mx2 Through-Going");
startZMC->Draw("HIST"); startZData->Draw("E0 P SAME");      tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");  lVertex->Draw("SAME");
c1.Print("startZDataSPINE.png"); c1.Print("startZDataSPINE.pdf");

startYMC->Scale(totalPOTData->Integral()/totalPOTMC->Integral());
startYMC->GetYaxis()->SetRangeUser(0,startYMC->GetMaximum()*2);
startYMC->GetXaxis()->CenterTitle();
startYMC->GetYaxis()->CenterTitle();
startYMC->SetLineColor(kRed);
startYMC->SetMarkerColor(kRed);
startYMC->GetXaxis()->SetTitle("Starting Position Y [cm]");
startYMC->GetYaxis()->SetTitle("Number of Tracks");
startYMC->SetTitle("SPINE-Mx2 Through-Going");
startYMC->Draw("HIST"); startYData->Draw("E0 P SAME");      tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}"); lVertex->Draw("SAME");
c1.Print("startYDataSPINE.png"); c1.Print("startYDataSPINE.pdf");


endXMC->Scale(totalPOTData->Integral()/totalPOTMC->Integral());
endXMC->GetYaxis()->SetRangeUser(0,endXMC->GetMaximum()*2);
endXMC->GetXaxis()->CenterTitle();
endXMC->GetYaxis()->CenterTitle();
endXMC->SetLineColor(kRed);
endXMC->SetMarkerColor(kRed);
endXMC->GetXaxis()->SetTitle("Ending Poisiton X [cm]");
endXMC->GetYaxis()->SetTitle("Number of Tracks");
endXMC->SetTitle("SPINE-Mx2 Through-Going");
endXMC->Draw("HIST"); endXData->Draw("E0 P SAME");      tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");  lVertex->Draw("SAME");
c1.Print("endXDataSPINE.png"); c1.Print("endXDataSPINE.pdf");


endZMC->Scale(totalPOTData->Integral()/totalPOTMC->Integral());
endZMC->GetYaxis()->SetRangeUser(0,endZMC->GetMaximum()*2);

endZMC->GetXaxis()->CenterTitle();
endZMC->GetYaxis()->CenterTitle();
endZMC->SetLineColor(kRed);
endZMC->SetMarkerColor(kRed);
endZMC->GetXaxis()->SetTitle("Ending Position Z [cm]");
endZMC->GetYaxis()->SetTitle("Number of Tracks");
endZMC->SetTitle("SPINE-Mx2 Through-Going");
endZMC->Draw("HIST"); endZData->Draw("E0 P SAME");      tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");  lVertex->Draw("SAME");
c1.Print("endZDataSPINE.png"); c1.Print("endZDataSPINE.pdf");

endYMC->Scale(totalPOTData->Integral()/totalPOTMC->Integral());
endYMC->GetYaxis()->SetRangeUser(0,endYMC->GetMaximum()*2);
endYMC->GetXaxis()->CenterTitle();
endYMC->GetYaxis()->CenterTitle();
endYMC->SetLineColor(kRed);
endYMC->SetMarkerColor(kRed);
endYMC->GetXaxis()->SetTitle("Ending Position Y [cm]");
endYMC->GetYaxis()->SetTitle("Number of Tracks");
endYMC->SetTitle("SPINE-Mx2 Through-Going");
endYMC->Draw("HIST"); endYData->Draw("E0 P SAME");      tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}"); lVertex->Draw("SAME");
c1.Print("endYDataSPINE.png"); c1.Print("endYDataSPINE.pdf");


startData->GetXaxis()->SetRangeUser(-60,60);
startData->GetYaxis()->SetRangeUser(-60,60);
startData->GetXaxis()->CenterTitle();
startData->GetYaxis()->CenterTitle();
startData->GetYaxis()->SetTitle("Entering Position Y [cm]");
startData->GetXaxis()->SetTitle("Entering Position X [cm]");
startData->GetZaxis()->SetTitle("Number of Tracks");
startData->SetTitle("SPINE-Mx2 Through-Going");
startData->Draw("COLZ");    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}"); 
c1.Print("startPositionDataUSData.png"); c1.Print("startPositionDataUSDataSPINE.pdf");

endData->GetXaxis()->SetRangeUser(-60,60);
endData->GetYaxis()->SetRangeUser(-60,60);
endData->GetXaxis()->CenterTitle();
endData->GetYaxis()->CenterTitle();
endData->GetYaxis()->SetTitle("Entering Position Y [cm]");
endData->GetXaxis()->SetTitle("Entering Position X [cm]");
endData->GetZaxis()->SetTitle("Number of Tracks");
endData->SetTitle("SPINE-Mx2 Through-Going");
endData->Draw("COLZ");    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}"); 
c1.Print("exitPositionDataUSData.png"); c1.Print("exitPositionDataUSDataSPINE.pdf");




}