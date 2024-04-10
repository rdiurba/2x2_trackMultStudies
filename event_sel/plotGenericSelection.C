	
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





}

