	
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

void plotMinervaData()
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
gPad->SetLeftMargin(0.18);
gPad->SetRightMargin(0.15);
TLatex tL;
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
TFile fMx2("testMinervaPandora.root");
TFile fMx2Data("testMinervaPandoraData.root");

TH1D* totalPOTMC=(TH1D*)fMx2.Get("totalPOT");
totalPOTMC->SetName("totalPOTMC");


TH1D* totalPOTData=(TH1D*)fMx2Data.Get("totalPOT");


TH1D* dotProductMC=(TH1D*)fMx2.Get("dotProductUS");
dotProductMC->SetName("dotProductUSMC");
TH1D* dotProductData=(TH1D*)fMx2Data.Get("dotProductUS");

      TLegend *lVertex = new TLegend(0.3,0.6,0.5,0.85);
  lVertex->SetTextFont(133);
  lVertex->SetTextSize(25);
  lVertex->SetHeader(Form("%1.2fE17 POT",totalPOTData->GetBinContent(1)/1E4));
  lVertex->AddEntry(dotProductData,"Data","lp");
  lVertex->AddEntry(dotProductMC,"Simulation","l");

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

    TH2D* startMC=(TH2D*)fMx2.Get("trackStartUS");
TH2D* endMC=(TH2D*)fMx2.Get("trackEndUS");
endMC->SetName("endMC");
startMC->SetName("startMC");

        TH2D* mx2StartMC=(TH2D*)fMx2.Get("mx2StartUS");
TH2D* mx2EndMC=(TH2D*)fMx2.Get("mx2EndUS");
mx2EndMC->SetName("mx2endMC");
mx2StartMC->SetName("mx2startMC");

    
TH1D* endXData=(TH1D*)fMx2Data.Get("trackEndXUS");
TH1D* endYData=(TH1D*)fMx2Data.Get("trackEndYUS");
TH1D* endZData=(TH1D*)fMx2Data.Get("trackEndZUS");


TH1D* startXData=(TH1D*)fMx2Data.Get("trackStartXUS");
TH1D* startYData=(TH1D*)fMx2Data.Get("trackStartYUS");
TH1D* startZData=(TH1D*)fMx2Data.Get("trackStartZUS");
    
TH2D* startData=(TH2D*)fMx2Data.Get("trackStartUS");
TH2D* endData=(TH2D*)fMx2Data.Get("trackEndUS");

    TH2D* mx2StartData=(TH2D*)fMx2Data.Get("mx2StartUS");
TH2D* mx2EndData=(TH2D*)fMx2Data.Get("mx2EndUS");


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

dotProductMC->GetXaxis()->CenterTitle();
dotProductMC->GetYaxis()->CenterTitle();
dotProductMC->SetLineColor(kRed);
dotProductMC->SetMarkerColor(kRed);



dotProductMC->Scale(totalPOTData->Integral()/totalPOTMC->Integral());
dotProductMC->GetYaxis()->SetRangeUser(0,dotProductMC->GetMaximum()*2);
dotProductMC->GetXaxis()->SetTitle("cos(#theta)");
dotProductMC->GetYaxis()->SetTitle("Number of Tracks");
dotProductMC->SetTitle("Pandora-Mx2 Through-Going");
dotProductMC->Draw("HIST"); dotProductData->Draw("E0 P SAME");      tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}"); lVertex->Draw("SAME");
c1.Print("pandoraPlots/dotProductData.png"); c1.Print("pandoraPlots/dotProductData.pdf");


startXMC->Scale(totalPOTData->Integral()/totalPOTMC->Integral());
startXMC->GetYaxis()->SetRangeUser(0,startXMC->GetMaximum()*2);
startXMC->GetXaxis()->CenterTitle();
startXMC->GetYaxis()->CenterTitle();
startXMC->SetLineColor(kRed);
startXMC->SetMarkerColor(kRed);
startXMC->GetXaxis()->SetTitle("Starting Poisiton X [cm]");
startXMC->GetYaxis()->SetTitle("Number of Tracks");
startXMC->SetTitle("Pandora-Mx2 Through-Going");
startXMC->Draw("HIST"); startXData->Draw("E0 P SAME");      tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");  lVertex->Draw("SAME");
c1.Print("pandoraPlots/startXData.png"); c1.Print("pandoraPlots/startXData.pdf");


startZMC->Scale(totalPOTData->Integral()/totalPOTMC->Integral());
startZMC->GetYaxis()->SetRangeUser(0,startZMC->GetMaximum()*2);

startZMC->GetXaxis()->CenterTitle();
startZMC->GetYaxis()->CenterTitle();
startZMC->SetLineColor(kRed);
startZMC->SetMarkerColor(kRed);
startZMC->GetXaxis()->SetTitle("Starting Position Z [cm]");
startZMC->GetYaxis()->SetTitle("Number of Tracks");
startZMC->SetTitle("Pandora-Mx2 Through-Going");
startZMC->Draw("HIST"); startZData->Draw("E0 P SAME");      tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");  lVertex->Draw("SAME");
c1.Print("pandoraPlots/startZData.png"); c1.Print("pandoraPlots/startZData.pdf");

startYMC->Scale(totalPOTData->Integral()/totalPOTMC->Integral());
startYMC->GetYaxis()->SetRangeUser(0,startYMC->GetMaximum()*2);
startYMC->GetXaxis()->CenterTitle();
startYMC->GetYaxis()->CenterTitle();
startYMC->SetLineColor(kRed);
startYMC->SetMarkerColor(kRed);
startYMC->GetXaxis()->SetTitle("Starting Position Y [cm]");
startYMC->GetYaxis()->SetTitle("Number of Tracks");
startYMC->SetTitle("Pandora-Mx2 Through-Going");
startYMC->Draw("HIST"); startYData->Draw("E0 P SAME");      tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}"); lVertex->Draw("SAME");
c1.Print("pandoraPlots/startYData.png"); c1.Print("pandoraPlots/startYData.pdf");


endXMC->Scale(totalPOTData->Integral()/totalPOTMC->Integral());
endXMC->GetYaxis()->SetRangeUser(0,endXMC->GetMaximum()*2);
endXMC->GetXaxis()->CenterTitle();
endXMC->GetYaxis()->CenterTitle();
endXMC->SetLineColor(kRed);
endXMC->SetMarkerColor(kRed);
endXMC->GetXaxis()->SetTitle("Ending Poisiton X [cm]");
endXMC->GetYaxis()->SetTitle("Number of Tracks");
endXMC->SetTitle("Pandora-Mx2 Through-Going");
endXMC->Draw("HIST"); endXData->Draw("E0 P SAME");      tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");  lVertex->Draw("SAME");
c1.Print("pandoraPlots/endXData.png"); c1.Print("pandoraPlots/endXData.pdf");


endZMC->Scale(totalPOTData->Integral()/totalPOTMC->Integral());
endZMC->GetYaxis()->SetRangeUser(0,endZMC->GetMaximum()*2);

endZMC->GetXaxis()->CenterTitle();
endZMC->GetYaxis()->CenterTitle();
endZMC->SetLineColor(kRed);
endZMC->SetMarkerColor(kRed);
endZMC->GetXaxis()->SetTitle("Ending Position Z [cm]");
endZMC->GetYaxis()->SetTitle("Number of Tracks");
endZMC->SetTitle("Pandora-Mx2 Through-Going");
endZMC->Draw("HIST"); endZData->Draw("E0 P SAME");      tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");  lVertex->Draw("SAME");
c1.Print("pandoraPlots/endZData.png"); c1.Print("pandoraPlots/endZData.pdf");

endYMC->Scale(totalPOTData->Integral()/totalPOTMC->Integral());
endYMC->GetYaxis()->SetRangeUser(0,endYMC->GetMaximum()*2);
endYMC->GetXaxis()->CenterTitle();
endYMC->GetYaxis()->CenterTitle();
endYMC->SetLineColor(kRed);
endYMC->SetMarkerColor(kRed);
endYMC->GetXaxis()->SetTitle("Ending Position Y [cm]");
endYMC->GetYaxis()->SetTitle("Number of Tracks");
endYMC->SetTitle("Pandora-Mx2 Through-Going");
endYMC->Draw("HIST"); endYData->Draw("E0 P SAME");      tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}"); lVertex->Draw("SAME");
c1.Print("pandoraPlots/endYData.png"); c1.Print("pandoraPlots/endYData.pdf");


startData->GetXaxis()->SetRangeUser(-60,60);
startData->GetYaxis()->SetRangeUser(-60,60);
startData->GetXaxis()->CenterTitle();
startData->GetYaxis()->CenterTitle();
startData->GetYaxis()->SetTitle("Entering Position Y [cm]");
startData->GetXaxis()->SetTitle("Entering Position X [cm]");
startData->GetZaxis()->SetTitle("Number of Tracks");
startData->SetTitle("Pandora-Mx2 Through-Going");
startData->Draw("COLZ");    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}"); 
c1.Print("pandoraPlots/startPositionDataUS.png"); c1.Print("pandoraPlots/startPositionDataUS.pdf");

endData->GetXaxis()->SetRangeUser(-60,60);
endData->GetYaxis()->SetRangeUser(-60,60);
endData->GetXaxis()->CenterTitle();
endData->GetYaxis()->CenterTitle();
endData->GetYaxis()->SetTitle("Exiting Position Y [cm]");
endData->GetXaxis()->SetTitle("Exiting Position X [cm]");
endData->GetZaxis()->SetTitle("Number of Tracks");
endData->SetTitle("Pandora-Mx2 Through-Going");
endData->Draw("COLZ");    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}"); 
c1.Print("pandoraPlots/exitPositionDataUS.png"); c1.Print("pandoraPlots/exitPositionDataUS.pdf");

    std::cout<<startData->GetNbinsX()<<std::endl;
    std::cout<<"Total Efficiency"<<std::endl;
std::cout<<startData->Integral()/mx2StartData->Integral()<<","<<endData->Integral()/mx2EndData->Integral()<<std::endl;
    std::cout<<"Efficiency x<0"<<std::endl;

std::cout<<startData->Integral(1,7,1,15)/mx2StartData->Integral(1,7,1,15)<<","<<endData->Integral(1,7,1,15)/mx2EndData->Integral(1,7,1,15)<<std::endl;
        std::cout<<"Efficiency x>0"<<std::endl;

std::cout<<startData->Integral(9,15,1,15)/mx2StartData->Integral(9,15,1,15)<<","<<endData->Integral(9,15,1,15)/mx2EndData->Integral(9,15,1,15)<<std::endl;

    
startData->Divide(mx2StartData);
endData->Divide(mx2EndData);


startData->GetXaxis()->SetRangeUser(-60,60);
startData->GetYaxis()->SetRangeUser(-60,60);
startData->GetXaxis()->CenterTitle();
startData->GetYaxis()->CenterTitle();
startData->GetYaxis()->SetTitle("Entering Position Y [cm]");
startData->GetXaxis()->SetTitle("Entering Position X [cm]");
startData->GetZaxis()->SetTitle("Matched TPC Tracks/Mx2 Tracks");
    startData->SetTitle("Pandora-Mx2 Through-Going");
    startData->GetZaxis()->SetRangeUser(0,1.0);
startData->Draw("COLZ");    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}"); 
c1.Print("pandoraPlots/startPositionDataUSEff.png"); c1.Print("pandoraPlots/startPositionDataUSEff.pdf");

endData->GetXaxis()->SetRangeUser(-60,60);
endData->GetYaxis()->SetRangeUser(-60,60);
endData->GetXaxis()->CenterTitle();
endData->GetYaxis()->CenterTitle();

endData->GetZaxis()->SetTitle("Matched TPC Tracks/Mx2 Tracks");
endData->SetTitle("Pandora-Mx2 Through-Going");
        endData->GetZaxis()->SetRangeUser(0,1.0);

endData->Draw("COLZ");    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}"); 
c1.Print("pandoraPlots/exitPositionDataUSEff.png"); c1.Print("pandoraPlots/exitPositionDataUSEff.pdf");

std::cout<<startMC->GetNbinsX()<<std::endl;
    std::cout<<startMC->Integral()<<","<<mx2StartMC->Integral()<<std::endl;

std::cout<<startMC->Integral()/mx2StartMC->Integral()<<","<<endMC->Integral()/mx2EndMC->Integral()<<std::endl;
    startMC->Divide(mx2StartMC);
endMC->Divide(mx2EndMC);


startMC->GetXaxis()->SetRangeUser(-60,60);
startMC->GetYaxis()->SetRangeUser(-60,60);
startMC->GetXaxis()->CenterTitle();
startMC->GetYaxis()->CenterTitle();
startMC->GetYaxis()->SetTitle("Entering Position Y [cm]");
startMC->GetXaxis()->SetTitle("Entering Position X [cm]");
startMC->GetZaxis()->SetTitle("Matched TPC Tracks/Mx2 Tracks");
    startMC->SetTitle("Pandora Through-Going MC");
    startMC->GetZaxis()->SetRangeUser(0,1.0);
startMC->Draw("COLZ");    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}"); 
c1.Print("pandoraPlots/startPositionMCUSEff.png"); c1.Print("pandoraPlots/startPositionMCUSEff.pdf");

endMC->GetXaxis()->SetRangeUser(-60,60);
endMC->GetYaxis()->SetRangeUser(-60,60);
endMC->GetXaxis()->CenterTitle();
endMC->GetYaxis()->CenterTitle();
endMC->GetYaxis()->SetTitle("Exiting Position Y [cm]");
endMC->GetXaxis()->SetTitle("Exiting Position X [cm]");
endMC->GetZaxis()->SetTitle("Matched TPC Tracks/Mx2 Tracks");
endMC->SetTitle("Pandora Through-Going MC");
        endMC->GetZaxis()->SetRangeUser(0,1.0);

endMC->Draw("COLZ");    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}"); 
c1.Print("pandoraPlots/exitPositionMCUSEff.png"); c1.Print("pandoraPlots/exitPositionMCUSEff.pdf");


}
