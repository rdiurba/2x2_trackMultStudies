	
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

void plotGenericSelection()
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
TFile fPurity("testPurity.root");
TFile fBad("testNoMinerva.root");
TFile fMinerva("testMinerva.root");

TH1D* backtrackEl=(TH1D*)fPurity.Get("recoBacktrackElAr");
backtrackEl->GetXaxis()->SetTitle("E_{#mu} [GeV]");
backtrackEl->GetYaxis()->SetTitle("Number of Interactions");
backtrackEl->SetTitle("Selected Interactions");
backtrackEl->Draw("E0 HIST");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("backtrackEl.png"); c1.Print("backtrackEl.pdf");

TH1D* backtrackCosl=(TH1D*)fPurity.Get("recoBacktrackCoslAr");
backtrackCosl->GetXaxis()->SetTitle("cos(#theta_{#mu})");
backtrackCosl->GetYaxis()->SetTitle("Number of Interactions");
backtrackCosl->SetTitle("Selected Interactions");
backtrackCosl->GetXaxis()->SetRangeUser(0.85,1);
backtrackCosl->Draw("E0 HIST");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("backtrackCosl.png"); c1.Print("backtrackCosl.pdf");



TH2D* confusionMatrix=(TH2D*)fPurity.Get("confusionMatrix");
confusionMatrix->SetTitle("Confusion Matrix for Particle ID");
confusionMatrix->GetXaxis()->SetTitle("Reco. Particle ID");
confusionMatrix->GetYaxis()->SetTitle("True Particle ID");
std::vector<string> labels={"Muon", "Proton","Pion","Other"};
for (int i=0; i<4; i++){
std::cout<<labels.at(i)<<std::endl;
confusionMatrix->GetXaxis()->SetBinLabel(i+1,labels.at(i).c_str());
confusionMatrix->GetYaxis()->SetBinLabel(i+1,labels.at(i).c_str());
double integral=confusionMatrix->Integral(i+1, i+1, 0,7);
for (int j=0; j<4; j++){
confusionMatrix->SetBinContent(i+1,j+1,confusionMatrix->GetBinContent(i+1,j+1)/integral);
}

}

confusionMatrix->SetMarkerSize(2.0);
confusionMatrix->GetXaxis()->CenterTitle();
confusionMatrix->GetYaxis()->CenterTitle();
confusionMatrix->GetYaxis()->LabelsOption("v");
gPad->SetLeftMargin(0.18);
gPad->SetRightMargin(0.15);
confusionMatrix->GetYaxis()->SetTitleOffset(1.5);
gStyle->SetPaintTextFormat("1.3f");
confusionMatrix->Draw("COLZ TEXT");
c1.Print("confusionMatrix.pdg");
c1.Print("confusionMatrix.png");

TH2D* responseMatrix=(TH2D*)fPurity.Get("responseMult");
responseMatrix->SetTitle("Response Matrix for Track Mult.");
responseMatrix->GetXaxis()->SetTitle("Reco. Charged Part.");
responseMatrix->GetYaxis()->SetTitle("True Charged Part.");

for (int i=0; i<20; i++){
double integral=responseMatrix->Integral(0,20, i+1,i+1);
for (int j=0; j<20; j++){
responseMatrix->SetBinContent(j+1,i+1,responseMatrix->GetBinContent(j+1,i+1)/integral);
}

}

responseMatrix->SetMarkerSize(2.0);
responseMatrix->GetXaxis()->CenterTitle();
responseMatrix->GetYaxis()->CenterTitle();
responseMatrix->GetXaxis()->SetRangeUser(1,8);
responseMatrix->GetYaxis()->SetRangeUser(1,8);
//responseMatrix->GetYaxis()->LabelsOption("v");
gPad->SetLeftMargin(0.18);
gPad->SetRightMargin(0.15);
responseMatrix->GetYaxis()->SetTitleOffset(1.2);
gStyle->SetPaintTextFormat("1.3f");
responseMatrix->Draw("COLZ TEXT");
c1.Print("responseMatrix.pdg");
c1.Print("responseMatrix.png");


TH1D* true_multTrkOnly=(TH1D*)fPurity.Get("true_multTrkOnly");
true_multTrkOnly->GetXaxis()->CenterTitle(); true_multTrkOnly->GetYaxis()->CenterTitle();
true_multTrkOnly->GetXaxis()->SetTitle("Number of True Tracks");
true_multTrkOnly->GetYaxis()->SetTitle("Number of Interactions"); true_multTrkOnly->SetTitle("Truth");
true_multTrkOnly->Draw("E0 HIST");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("true_multTrk.png"); c1.Print("true_multTrk.pdf");


TEfficiency* protonEff=(TEfficiency*)fPurity.Get("trueProtonEWithRecoInt_clone");
protonEff->SetTitle("Proton Eff. (GeV)");
protonEff->Draw();
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("protonEff.png");


TEfficiency* pionEff=(TEfficiency*)fPurity.Get("truePionEWithRecoInt_clone");
pionEff->SetTitle("Pion Eff. (GeV)");
pionEff->Draw();
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("pionEff.png");

TEfficiency* pionDirX=(TEfficiency*)fPurity.Get("truePionWithRecoIntDirX_clone");
pionDirX->SetTitle("Pion Eff. for cos(#theta_{X})");
pionDirX->Draw();
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("pionDirX.png");


TEfficiency* pionDirY=(TEfficiency*)fPurity.Get("truePionWithRecoIntDirY_clone");
pionDirY->SetTitle("Pion Eff. for cos(#theta_{Y})");
pionDirY->Draw();
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("pionDirY.png");


TEfficiency* pionDirZ=(TEfficiency*)fPurity.Get("truePionWithRecoIntDirZ_clone");
pionDirZ->SetTitle("Pion Eff. for cos(#theta_{Z})");
pionDirZ->Draw();
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("pionDirZ.png");


TEfficiency* protonDirX=(TEfficiency*)fPurity.Get("trueProtonWithRecoIntDirX_clone");
protonDirX->SetTitle("Proton Eff. for cos(#theta_{X})");
protonDirX->Draw();
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("protonDirX.png");


TEfficiency* protonDirY=(TEfficiency*)fPurity.Get("trueProtonWithRecoIntDirY_clone");
protonDirY->SetTitle("Proton Eff. for cos(#theta_{Y})");
protonDirY->Draw();
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("protonDirY.png");


TEfficiency* protonDirZ=(TEfficiency*)fPurity.Get("trueProtonWithRecoIntDirZ_clone");
protonDirZ->SetTitle("Proton Eff. for cos(#theta_{Z})");
protonDirZ->Draw();
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("protonDirZ.png");





TH1D* track_mult=(TH1D*)fPurity.Get("track_mult");
TH1D* track_multQE=(TH1D*)fPurity.Get("track_multQE");
track_multQE->SetFillColor(kRed); track_multQE->SetLineColor(kRed); track_multQE->SetLineWidth(0);
TH1D* track_multMEC=(TH1D*)fPurity.Get("track_multMEC");
track_multMEC->SetFillColor(kBlue); track_multMEC->SetLineColor(kBlue); track_multMEC->SetLineWidth(0);
TH1D* track_multDIS=(TH1D*)fPurity.Get("track_multDIS");
track_multDIS->SetFillColor(kYellow); track_multDIS->SetLineColor(kYellow); 
TH1D* track_multRES=(TH1D*)fPurity.Get("track_multRES");
track_multRES->SetFillColor(kGreen+2); track_multRES->SetLineColor(kGreen+2); 
TH1D* track_multNC=(TH1D*)fPurity.Get("track_multNC");
track_multNC->SetFillColor(kGray); track_multNC->SetLineColor(kGray);
TH1D* track_multCOH=(TH1D*)fPurity.Get("track_multCOH");
track_multCOH->SetFillColor(kCyan); track_multCOH->SetLineColor(kCyan);
TH1D* track_multRock=(TH1D*)fPurity.Get("track_multRock");
track_multRock->SetFillColor(kViolet); track_multRock->SetLineColor(kViolet);
TH1D* track_multSec=(TH1D*)fPurity.Get("track_multSec");
track_multSec->SetFillColor(kBlack); track_multSec->SetLineColor(kBlack);
track_multRES->SetLineWidth(0);
track_multDIS->SetLineWidth(0);
track_multNC->SetLineWidth(0);
track_multCOH->SetLineWidth(0);
track_multRock->SetLineWidth(0);
track_multSec->SetLineWidth(0);
  TLegend *lTrack = new TLegend(0.5,0.3,0.8,0.9);
  lTrack->SetTextFont(133);
  lTrack->SetTextSize(25);
  lTrack->SetHeader("RHC #nu_{#mu}/#bar{#nu_{#mu}}");
  lTrack->AddEntry(track_multQE,"CC-QE","f");
  lTrack->AddEntry(track_multMEC,"CC-MEC","f");
  lTrack->AddEntry(track_multDIS,"CC-DIS","f");
  lTrack->AddEntry(track_multRES,"CC-RES","f");
  lTrack->AddEntry(track_multCOH,"CC-COH","f");
  lTrack->AddEntry(track_multNC,"NC","f");
  lTrack->AddEntry(track_multRock,"Rock","f");
  lTrack->AddEntry(track_multSec,"Secondary Int.","f");
THStack* h1d_track=new THStack("hs","");
h1d_track->Add(track_multQE);
h1d_track->Add(track_multMEC);
h1d_track->Add(track_multDIS);
h1d_track->Add(track_multRES);
h1d_track->Add(track_multCOH);
h1d_track->Add(track_multNC);
h1d_track->Add(track_multRock);
h1d_track->Add(track_multSec);







TH1D* track_multGood=(TH1D*)fPurity.Get("track_multGood");
track_mult->GetXaxis()->CenterTitle();
track_mult->GetXaxis()->SetTitle("Number of Reconstructed Tracks");
track_mult->GetYaxis()->SetTitle("Number of Interactions");
track_mult->GetYaxis()->CenterTitle();
track_multGood->GetXaxis()->CenterTitle();
track_multGood->GetYaxis()->CenterTitle();
track_multGood->GetXaxis()->SetTitle("Number of Reconstructed Tracks");
track_multGood->GetYaxis()->SetTitle("Number of Interactions");
track_multGood->SetTitle("Selected Interactions: Cheated");
track_multGood->SetLineColor(kRed);
track_multGood->SetMarkerColor(kRed);
track_mult->SetLineColor(kRed);
track_mult->SetMarkerColor(kRed);
track_mult->SetTitle("Selected Interactions");
track_mult->Draw("E0 HIST");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("track_mult.png"); c1.Print("track_mult.pdf");
track_mult->GetYaxis()->SetRangeUser(0,1000);
track_mult->SetLineWidth(0); track_mult->GetXaxis()->SetRangeUser(0,10);
track_mult->Draw("HIST");
h1d_track->Draw("SAME");
lTrack->Draw("SAME");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("track_multStack.png"); c1.Print("track_multStack.pdf");
track_mult->SetLineWidth(track_multGood->GetLineWidth());


track_multGood->Draw("E0 HIST");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("track_multGood.png"); c1.Print("track_multGood.pdf");


track_mult->SetLineWidth(0); track_mult->GetXaxis()->SetRangeUser(0,10); track_mult->SetTitle("Scaled to 1 Day");
track_mult->Scale(24.0/82.0);
track_multQE->Scale(24.0/82.0);
track_multMEC->Scale(24.0/82.0);
track_multDIS->Scale(24.0/82.0);
track_multRES->Scale(24.0/82.0);
track_multCOH->Scale(24.0/82.0);
track_multNC->Scale(24.0/82.0);
track_multRock->Scale(24.0/82.0);
track_multSec->Scale(24.0/82.0);

track_multQE->SetFillColor(kRed); track_multQE->SetLineColor(kRed); track_multQE->SetLineWidth(0);

track_multMEC->SetFillColor(kBlue); track_multMEC->SetLineColor(kBlue); track_multMEC->SetLineWidth(0);

track_multDIS->SetFillColor(kYellow); track_multDIS->SetLineColor(kYellow); 

track_multRES->SetFillColor(kGreen+2); track_multRES->SetLineColor(kGreen+2); 

track_multNC->SetFillColor(kGray); track_multNC->SetLineColor(kGray);

track_multCOH->SetFillColor(kCyan); track_multCOH->SetLineColor(kCyan);

track_multRock->SetFillColor(kViolet); track_multRock->SetLineColor(kViolet);

track_multSec->SetFillColor(kBlack); track_multSec->SetLineColor(kBlack);
THStack* h1d_track2=new THStack("hs2","hs2");
h1d_track2->Add(track_multQE);
h1d_track2->Add(track_multMEC);
h1d_track2->Add(track_multDIS);
h1d_track2->Add(track_multRES);
h1d_track2->Add(track_multCOH);
h1d_track2->Add(track_multNC);
h1d_track2->Add(track_multRock);
h1d_track2->Add(track_multSec);
std::cout<<track_multQE->GetMaximum()<<std::endl;

track_mult->Draw("HIST");
h1d_track2->Draw("SAME HIST");
lTrack->Draw("SAME");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("track_multStack1Day.png"); c1.Print("track_multStack1Day.pdf");
track_mult->SetLineWidth(track_multGood->GetLineWidth());

track_mult->SetLineWidth(0); track_mult->GetXaxis()->SetRangeUser(0,10); track_mult->SetTitle("Scaled to 10 Days");
track_mult->Scale(2.927*82.0/24.0);
track_multQE->Scale(2.927*82.0/24.0);
track_multMEC->Scale(2.927*82.0/24.0);
track_multDIS->Scale(2.9277*82.0/24.0);
track_multRES->Scale(2.927*82.0/24.0);
track_multCOH->Scale(2.927*82.0/24.0);
track_multNC->Scale(2.927*82.0/24.0);
track_multRock->Scale(2.927*82.0/24.0);
track_multSec->Scale(2.927*82.0/24.0);


track_multQE->SetFillColor(kRed); track_multQE->SetLineColor(kRed); track_multQE->SetLineWidth(0);

track_multMEC->SetFillColor(kBlue); track_multMEC->SetLineColor(kBlue); track_multMEC->SetLineWidth(0);

track_multDIS->SetFillColor(kYellow); track_multDIS->SetLineColor(kYellow); 

track_multRES->SetFillColor(kGreen+2); track_multRES->SetLineColor(kGreen+2); 

track_multNC->SetFillColor(kGray); track_multNC->SetLineColor(kGray);

track_multCOH->SetFillColor(kCyan); track_multCOH->SetLineColor(kCyan);

track_multRock->SetFillColor(kViolet); track_multRock->SetLineColor(kViolet);

track_multSec->SetFillColor(kBlack); track_multSec->SetLineColor(kBlack);
std::cout<<track_multQE->GetMaximum()<<std::endl;
THStack* h1d_track3=new THStack("hs3","hs3");
h1d_track3->Add(track_multQE);
h1d_track3->Add(track_multMEC);
h1d_track3->Add(track_multDIS);
h1d_track3->Add(track_multRES);
h1d_track3->Add(track_multCOH);
h1d_track3->Add(track_multNC);
h1d_track3->Add(track_multRock);
h1d_track3->Add(track_multSec);
track_mult->GetYaxis()->SetTitleOffset(1.4);
//track_mult->GetYaxis()->SetRangeUser(0,600*240.0/82.0);
track_mult->Draw("HIST");
h1d_track3->Draw("SAME HIST");
lTrack->Draw("SAME");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("track_multStack10Days.png"); c1.Print("track_multStack10Days.pdf");
track_mult->SetLineWidth(track_multGood->GetLineWidth());



track_mult->Scale(82.0/240.0);
track_multQE->Scale(82.0/240.0);
track_multMEC->Scale(82.0/240.0);
track_multDIS->Scale(82.0/240.0);
track_multRES->Scale(82.0/240.0);
track_multCOH->Scale(82.0/240.0);
track_multNC->Scale(82.0/240.0);
track_multRock->Scale(82.0/240.0);
track_multSec->Scale(82.0/240.0);
track_multRock->Scale(2);
track_multQE->SetFillColor(kRed); track_multQE->SetLineColor(kRed); track_multQE->SetLineWidth(0);

track_multMEC->SetFillColor(kBlue); track_multMEC->SetLineColor(kBlue); track_multMEC->SetLineWidth(0);

track_multDIS->SetFillColor(kYellow); track_multDIS->SetLineColor(kYellow);

track_multRES->SetFillColor(kGreen+2); track_multRES->SetLineColor(kGreen+2);

track_multNC->SetFillColor(kGray); track_multNC->SetLineColor(kGray);

track_multCOH->SetFillColor(kCyan); track_multCOH->SetLineColor(kCyan);

track_multRock->SetFillColor(kViolet); track_multRock->SetLineColor(kViolet);

track_multSec->SetFillColor(kBlack); track_multSec->SetLineColor(kBlack);
track_mult->SetLineWidth(0);
track_mult->GetYaxis()->SetRangeUser(0,1000);
THStack* h1d_track4=new THStack("hs4","hs4");
h1d_track4->Add(track_multQE);
h1d_track4->Add(track_multMEC);
h1d_track4->Add(track_multDIS);
h1d_track4->Add(track_multRES);
h1d_track4->Add(track_multCOH);
h1d_track4->Add(track_multNC);
h1d_track4->Add(track_multRock);
h1d_track4->Add(track_multSec);

track_mult->SetTitle("Scaled Rock by 2");

track_mult->Draw("HIST");
h1d_track4->Draw("SAME HIST");
lTrack->Draw("SAME");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("track_multStackRockTimes2.png"); c1.Print("track_multStackRockTimes2.pdf");
track_mult->SetLineWidth(track_multGood->GetLineWidth());




TH1D* rock=(TH1D*)fPurity.Get("rock");
rock->GetXaxis()->CenterTitle();
rock->SetLineColor(kRed);
rock->SetMarkerColor(kRed);
rock->GetXaxis()->SetBinLabel(1,"Non-Rock");
rock->GetXaxis()->SetBinLabel(2,"Rock");
rock->GetYaxis()->SetTitle("Number of Interactions");
rock->SetTitle("Classification of Background Int.");
rock->GetYaxis()->CenterTitle();
rock->Draw("E0 P");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("rock.png"); c1.Print("rock.pdf");


TH1D* rockWithoutMNV=(TH1D*)fBad.Get("rock");
rockWithoutMNV->GetXaxis()->CenterTitle();
rockWithoutMNV->GetXaxis()->SetTitle("Background Event Comes from Rock");
rockWithoutMNV->GetXaxis()->SetBinLabel(1,"Non-Rock");
rockWithoutMNV->GetXaxis()->SetBinLabel(2,"Rock");
rockWithoutMNV->GetYaxis()->SetTitle("Number of Interactions");
rockWithoutMNV->SetTitle("Classification Without Using MINERvA");
rockWithoutMNV->GetYaxis()->CenterTitle();
rockWithoutMNV->SetLineColor(kRed);
rockWithoutMNV->SetMarkerColor(kRed);
rockWithoutMNV->Draw("E0 P");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("rockWithoutMNV.png"); c1.Print("rockWithoutMNV.pdf");


TH1D* track_multWithoutMNV=(TH1D*)fBad.Get("track_mult");
track_multWithoutMNV->SetLineColor(kRed);
track_multWithoutMNV->SetMarkerColor(kRed);
track_multWithoutMNV->GetXaxis()->CenterTitle();
track_multWithoutMNV->GetXaxis()->SetTitle("Number of Reconstructed Tracks");
track_multWithoutMNV->GetYaxis()->SetTitle("Number of Interactions");
track_multWithoutMNV->GetYaxis()->CenterTitle();
track_multWithoutMNV->SetTitle("Without Using MINERvA");

track_multWithoutMNV->Draw("E0 HIST");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("track_multWithoutMNV.png"); c1.Print("track_multWithoutMNV.pdf");


TH1D* deltaX=(TH1D*)fMinerva.Get("deltaX");
TH1D* deltaY=(TH1D*)fMinerva.Get("deltaY");
deltaX->GetXaxis()->CenterTitle();
deltaX->GetXaxis()->SetTitle("#Delta X");
deltaX->GetYaxis()->SetTitle("Number of Tracks");
deltaX->SetTitle("DS MINERvA-TPC Matching");
deltaY->SetTitle("DS MINERvA-TPC Matching");
deltaX->SetLineColor(kRed);
deltaX->SetMarkerColor(kRed);
deltaX->GetYaxis()->CenterTitle();
deltaY->GetXaxis()->CenterTitle();
deltaY->GetYaxis()->CenterTitle();
deltaY->GetXaxis()->SetTitle("#Delta Y");
deltaY->GetYaxis()->SetTitle("Number of Tracks");
deltaY->SetLineColor(kRed);
deltaY->SetMarkerColor(kRed);


TH1D* dotProduct=(TH1D*)fMinerva.Get("dotProduct");
dotProduct->GetXaxis()->CenterTitle();
dotProduct->GetYaxis()->CenterTitle();
dotProduct->SetLineColor(kRed);
dotProduct->SetMarkerColor(kRed);
dotProduct->GetXaxis()->SetTitle("cos(#theta)");
dotProduct->GetYaxis()->SetTitle("Number of Interactions");
dotProduct->SetTitle("DS MINERvA-TPC Matching");


TH1D* deltaXGood=(TH1D*)fMinerva.Get("deltaXGood");
TH1D* deltaYGood=(TH1D*)fMinerva.Get("deltaYGood");

deltaXGood->GetXaxis()->CenterTitle();
deltaXGood->GetXaxis()->SetTitle("#Delta X");
deltaXGood->SetTitle("Matching: Cheated");
deltaYGood->SetTitle("Matching: Cheated");
deltaXGood->GetYaxis()->SetTitle("Number of Tracks");
deltaXGood->GetYaxis()->CenterTitle();
deltaYGood->GetXaxis()->CenterTitle();
deltaYGood->GetYaxis()->CenterTitle();
deltaYGood->GetXaxis()->SetTitle("#Delta Y");
deltaYGood->GetYaxis()->SetTitle("Number of Tracks");
deltaXGood->SetTitle("Cheated");
deltaYGood->SetTitle("Cheated");
deltaXGood->SetLineColor(kRed);
deltaXGood->SetMarkerColor(kRed);
deltaYGood->SetLineColor(kRed);
deltaYGood->SetMarkerColor(kRed);



TH1D* dotProductGood=(TH1D*)fMinerva.Get("dotProductGood");
dotProductGood->GetXaxis()->CenterTitle();
dotProductGood->GetYaxis()->CenterTitle();
dotProductGood->GetXaxis()->SetTitle("cos(#theta)");
dotProductGood->GetYaxis()->SetTitle("Number of Tracks");
dotProductGood->SetTitle("Matching: Cheated");
dotProductGood->SetLineColor(kRed);
dotProductGood->SetMarkerColor(kRed);

TH1D* deltaXBad=(TH1D*)fMinerva.Get("deltaXBad");
TH1D* deltaYBad=(TH1D*)fMinerva.Get("deltaYBad");
deltaXBad->SetTitle("DS MINERvA-TPC Matching: Background");
deltaYBad->SetTitle("DS MINERvA-TPC Matching: Background");

deltaXBad->GetXaxis()->CenterTitle();
deltaXBad->GetXaxis()->SetTitle("#Delta X");
deltaXBad->GetYaxis()->SetTitle("Number of Tracks");
deltaXBad->GetYaxis()->CenterTitle();
deltaYBad->GetXaxis()->CenterTitle();
deltaYBad->GetYaxis()->CenterTitle();
deltaYBad->GetXaxis()->SetTitle("#Delta Y");
deltaYBad->GetYaxis()->SetTitle("Number of Tracks");


deltaX->Draw("E0 P");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("deltaX.png"); c1.Print("deltaX.pdf");

deltaY->Draw("E0 P");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("deltaY.png"); c1.Print("deltaY.pdf");


deltaXGood->Draw("E0 P");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("deltaXGood.png"); c1.Print("deltaXGood.pdf");


deltaYGood->Draw("E0 P");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("deltaYGood.png"); c1.Print("deltaYGood.pdf");



deltaXBad->Draw("E0 P");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("deltaXBad.png"); c1.Print("deltaXBad.pdf");


deltaYBad->Draw("E0 P");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("deltaYBad.png"); c1.Print("deltaYBad.pdf");


dotProduct->Draw("E0 P");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("dotProduct.png"); c1.Print("dotProduct.pdf");


dotProductGood->Draw("E0 P");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("dotProductGood.png"); c1.Print("dotProductGood.pdf");



}

