	
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



gROOT->LoadMacro("protoDUNEStyle.C");
gROOT->SetStyle("protoDUNEStyle");
gROOT->ForceStyle();
gStyle->SetTitleX(0.35);
gStyle->SetOptFit(111);


TCanvas c1=TCanvas();

TLatex tL;
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
TFile fPurity("testPurityMiniRun5.root");
TFile fWion22("testOutputWion_22.7.root");
TFile fWion25("testOutputWion_25.1.root");

//TFile fPurity("trkEff1Sigma.root");
TFile fBad("testNoMinerva.root");
TFile fMinerva("testMinerva.root");


TH1D* diffAngle=(TH1D*)fPurity.Get("trueDiffPosvsPDirZ");
diffAngle->GetYaxis()->SetTitle("Number of Interactions");
diffAngle->SetTitle("All True CC #nu_{#mu}");
diffAngle->GetXaxis()->SetTitle("Diff. True Angle (Pos.-Momentum)");
diffAngle->GetYaxis()->SetTitleOffset(1.4);
diffAngle->Draw("HIST");
c1.Print("diffAngleMR5.png");



TH1D* backtrackEl=(TH1D*)fPurity.Get("recoBacktrackElAr");
backtrackEl->GetXaxis()->SetTitle("True E_{#mu} [GeV]");
backtrackEl->GetYaxis()->SetTitle("Number of Interactions");
backtrackEl->SetTitle("Selected Interactions");
backtrackEl->Draw("E0 HIST");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("backtrackElMR5.png"); c1.Print("backtrackElMR5.pdf");

TH1D* backtrackCosl=(TH1D*)fPurity.Get("recoBacktrackCoslAr");
backtrackCosl->GetXaxis()->SetTitle("True cos(#theta_{#mu})");
backtrackCosl->GetYaxis()->SetTitle("Number of Interactions");
backtrackCosl->SetTitle("Selected Interactions");
backtrackCosl->GetXaxis()->SetRangeUser(0.85,1);
backtrackCosl->Draw("E0 HIST");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("backtrackCoslMR5.png"); c1.Print("backtrackCoslMR5.pdf");



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
c1.Print("confusionMatrixMR5.png");

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
c1.Print("responseMatrixMR5.png");


TH2D* responseMatrix22=(TH2D*)fWion22.Get("responseMult");
responseMatrix22->SetTitle("Response W_{ion}=22.7 eV");
responseMatrix22->GetXaxis()->SetTitle("Reco. Charged Part.");
responseMatrix22->GetYaxis()->SetTitle("True Charged Part.");

for (int i=0; i<20; i++){
double integral=responseMatrix22->Integral(0,20, i+1,i+1);
for (int j=0; j<20; j++){
responseMatrix22->SetBinContent(j+1,i+1,responseMatrix22->GetBinContent(j+1,i+1)/integral);
}

}

responseMatrix22->SetMarkerSize(2.0);
responseMatrix22->GetXaxis()->CenterTitle();
responseMatrix22->GetYaxis()->CenterTitle();
responseMatrix22->GetXaxis()->SetRangeUser(1,8);
responseMatrix22->GetYaxis()->SetRangeUser(1,8);
//responseMatrix->GetYaxis()->LabelsOption("v");
gPad->SetLeftMargin(0.18);
gPad->SetRightMargin(0.15);
responseMatrix22->GetYaxis()->SetTitleOffset(1.2);
gStyle->SetPaintTextFormat("1.3f");
responseMatrix22->Draw("COLZ TEXT");
c1.Print("responseMatrix22MR5.pdf");
c1.Print("responseMatrix22MR5.png");



TH2D* responseMatrixWion25=(TH2D*)fWion25.Get("responseMult");
responseMatrixWion25->SetTitle("Response W_{ion}=25.1 eV");
responseMatrixWion25->GetXaxis()->SetTitle("Reco. Charged Part.");
responseMatrixWion25->GetYaxis()->SetTitle("True Charged Part.");

for (int i=0; i<20; i++){
double integral=responseMatrixWion25->Integral(0,20, i+1,i+1);
for (int j=0; j<20; j++){
responseMatrixWion25->SetBinContent(j+1,i+1,responseMatrixWion25->GetBinContent(j+1,i+1)/integral);
}

}

responseMatrixWion25->SetMarkerSize(2.0);
responseMatrixWion25->GetXaxis()->CenterTitle();
responseMatrixWion25->GetYaxis()->CenterTitle();
responseMatrixWion25->GetXaxis()->SetRangeUser(1,8);
responseMatrixWion25->GetYaxis()->SetRangeUser(1,8);
//responseMatrix->GetYaxis()->LabelsOption("v");
gPad->SetLeftMargin(0.18);
gPad->SetRightMargin(0.15);
responseMatrixWion25->GetYaxis()->SetTitleOffset(1.2);
gStyle->SetPaintTextFormat("1.3f");
responseMatrixWion25->Draw("COLZ TEXT");
c1.Print("responseMatrixWion25MR5.pdf");
c1.Print("responseMatrixWion25MR5.png");



TH1D* true_multTrkOnly=(TH1D*)fPurity.Get("true_multTrkOnly");
true_multTrkOnly->GetXaxis()->CenterTitle(); true_multTrkOnly->GetYaxis()->CenterTitle();
true_multTrkOnly->GetXaxis()->SetTitle("Number of True Tracks");
true_multTrkOnly->GetYaxis()->SetTitle("Number of Interactions"); true_multTrkOnly->SetTitle("Truth");
true_multTrkOnly->Draw("E0 HIST");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("true_multTrkMR5.png"); c1.Print("true_multTrkMR5.pdf");


TEfficiency* trkLenEff=(TEfficiency*)fPurity.Get("trueTrkLen_clone");
trkLenEff->SetTitle("Trk. Len (cm)");
trkLenEff->Draw();
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("trueTrkLenEffMR5.png");

TEfficiency* protonEff=(TEfficiency*)fPurity.Get("trueProtonEWithRecoInt_clone");
TH1D* denomWion25=(TH1D*)fWion25.Get("trueProtonEWithRecoInt");
TH1D* denomWion22=(TH1D*)fWion22.Get("trueProtonEWithRecoInt");
TH1D* protonEffWion25=(TH1D*)fWion25.Get("selProtonE");
TH1D* protonEffWion22=(TH1D*)fWion22.Get("selProtonE");
protonEffWion22->SetName("protonEffWion22");
protonEffWion25->SetName("protonEffWion25");
protonEffWion22->Divide(denomWion22);
protonEffWion25->Divide(denomWion25);
protonEff->SetLineColor(kRed);
protonEff->SetMarkerColor(kRed);
protonEffWion22->SetMarkerColor(kBlue);
protonEffWion22->SetLineColor(kBlue);
protonEffWion25->SetLineColor(kBlack);
protonEffWion25->SetLineColor(kBlack);
protonEffWion22->GetXaxis()->SetRangeUser(0,0.30);
  TLegend *lSyst2 = new TLegend(0.2,0.7,0.3,0.9);
  lSyst2->SetTextFont(133);
  lSyst2->SetTextSize(25);

  lSyst2->AddEntry(protonEffWion25,"W_{ion}=25.1 eV","l");
  lSyst2->AddEntry(protonEff,"W_{ion}=23.6 eV","l");
  lSyst2->AddEntry(protonEffWion22,"W_{ion}=22.7 eV","l");

protonEffWion22->SetLineColor(0);


protonEffWion22->GetYaxis()->SetRangeUser(0,1);
protonEffWion22->SetTitle("Proton Eff. in KE (GeV)");
protonEffWion22->Draw("HIST");
protonEff->Draw("SAME");

//protonEffWion25->Draw("SAME HIST");
//lSyst2->Draw("SAME");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("protonEffMR5.png");


TEfficiency* pionEff=(TEfficiency*)fPurity.Get("truePionEWithRecoInt_clone");
TH1D* denomWion252=(TH1D*)fWion25.Get("truePionEWithRecoInt");
TH1D* denomWion222=(TH1D*)fWion22.Get("truePionEWithRecoInt");
TH1D* pionEffWion25=(TH1D*)fWion25.Get("selPionE");
TH1D* pionEffWion22=(TH1D*)fWion22.Get("selPionE");
pionEffWion22->SetName("pionEffWion22");
pionEffWion25->SetName("pionEffWion25");
pionEffWion22->Divide(denomWion222);
pionEffWion25->Divide(denomWion252);
pionEff->SetMarkerColor(kRed);
pionEff->SetLineColor(kRed);
pionEffWion22->SetMarkerColor(kBlue);
pionEffWion22->SetLineColor(kBlue);
pionEffWion25->SetLineColor(kBlack);
pionEffWion25->SetLineColor(kBlack);
pionEffWion22->SetTitle("Pion Eff. (GeV)");
pionEffWion22->GetYaxis()->SetRangeUser(0,1);
pionEffWion22->SetLineColor(0);
pionEffWion22->Draw("HIST");
pionEff->Draw("SAME");

    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
//lSyst2->Draw("SAME");
c1.Print("pionEffMR5.png");

TEfficiency* pionDirX=(TEfficiency*)fPurity.Get("truePionWithRecoIntDirX_clone");
pionDirX->SetTitle("Pion Eff. for cos(#theta_{X})");
pionDirX->Draw();
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("pionDirXMR5.png");


TEfficiency* pionDirY=(TEfficiency*)fPurity.Get("truePionWithRecoIntDirY_clone");
pionDirY->SetTitle("Pion Eff. for cos(#theta_{Y})");
pionDirY->Draw();
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("pionDirYMR5.png");


TEfficiency* pionDirZ=(TEfficiency*)fPurity.Get("truePionWithRecoIntDirZ_clone");
pionDirZ->SetTitle("Pion Eff. for cos(#theta_{Z})");
pionDirZ->Draw();
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("pionDirZMR5.png");


TEfficiency* protonDirX=(TEfficiency*)fPurity.Get("trueProtonWithRecoIntDirX_clone");
protonDirX->SetTitle("Proton Eff. for cos(#theta_{X})");
protonDirX->Draw();
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("protonDirXMR5.png");


TEfficiency* protonDirY=(TEfficiency*)fPurity.Get("trueProtonWithRecoIntDirY_clone");
protonDirY->SetTitle("Proton Eff. for cos(#theta_{Y})");
protonDirY->Draw();
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("protonDirYMR5.png");


TEfficiency* protonDirZ=(TEfficiency*)fPurity.Get("trueProtonWithRecoIntDirZ_clone");
protonDirZ->SetTitle("Proton Eff. for cos(#theta_{Z})");
protonDirZ->Draw();
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("protonDirZMR5.png");



TH1D* track_mult=(TH1D*)fPurity.Get("track_mult");
TH1D* track_multWion22=(TH1D*)fWion22.Get("track_mult");
track_multWion22->SetName("track_multWion22");
TH1D* track_multWion25=(TH1D*)fWion25.Get("track_mult");
track_multWion25->SetName("track_multWion25");
track_multWion22->SetLineStyle(2);
track_multWion25->SetLineStyle(3);
track_multWion22->SetLineColor(kBlue);
track_multWion25->SetLineColor(kBlack);
track_mult->SetLineColor(kRed);
track_mult->SetMarkerColor(kRed);
track_mult->GetYaxis()->SetRangeUser(0,1.2*track_multWion25->GetMaximum());
track_mult->SetTitle("Selected Int: Stat. Error");
track_mult->GetXaxis()->CenterTitle();
track_mult->GetXaxis()->SetTitle("Number of Reconstructed Tracks");
track_mult->GetYaxis()->SetTitle("Number of Interactions");
track_mult->GetYaxis()->CenterTitle();

track_mult->Draw("HIST E0");
track_multWion25->Draw("HIST SAME");
track_multWion22->Draw("HIST SAME");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
  TLegend *lSyst = new TLegend(0.3,0.5,0.8,0.9);
  lSyst->SetTextFont(133);
  lSyst->SetTextSize(25);
  lSyst->SetHeader("RHC #nu_{#mu}/#bar{#nu_{#mu}}");
  lSyst->AddEntry(track_multWion25,"W_{ion}=25.1 eV","l");
  lSyst->AddEntry(track_mult,"W_{ion}=23.6 eV (Nominal)","l");
  lSyst->AddEntry(track_multWion22,"W_{ion}=22.7 eV","l");
  lSyst->Draw("SAME");
    
c1.Print("track_multWionSystMR5.png"); c1.Print("track_multWionSystMR5.pdf");




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
  TLegend *lTrack = new TLegend(0.5,0.5,0.8,0.9);
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
track_mult->SetTitle("Selected Int: Stat. Error");
track_mult->GetYaxis()->SetRangeUser(0,800);
track_mult->Draw("E0 HIST");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");

c1.Print("track_multMR5.png"); c1.Print("track_multMR5.pdf");




track_mult->GetYaxis()->SetRangeUser(0,2000);
track_mult->SetLineWidth(0); track_mult->GetXaxis()->SetRangeUser(0,10);
track_mult->Draw("HIST");
track_mult->GetXaxis()->SetTitle("Number of Reconstructed Tracks");
track_mult->GetYaxis()->SetTitle("Number of Interactions");
track_mult->SetTitle("Selected Interactions");
h1d_track->Draw("SAME");
lTrack->Draw("SAME");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("track_multStackMR5.png"); c1.Print("track_multStackMR5.pdf");
track_mult->SetLineWidth(track_multGood->GetLineWidth());



std::vector<double> diffCCDIS;
std::vector<double> diffCCQE;
std::vector<double> diffCCMEC;
std::vector<double> diffCCRES;
std::vector<double> diffRock;
std::vector<double> diffTrkEff;
TFile fRW("trkEff1Sigma.root");
TH1D* rwTracks=(TH1D*)fRW.Get("track_multRW");
for (int i=0; i<track_mult->GetNbinsX(); i++){
diffRock.push_back(track_multRock->GetBinContent(i)*2-track_multRock->GetBinContent(i));
diffTrkEff.push_back(abs(rwTracks->GetBinContent(i)-track_mult->GetBinContent(i)));
diffCCDIS.push_back(track_multDIS->GetBinContent(i)-track_multDIS->GetBinContent(i)*0.5);
diffCCMEC.push_back(track_multMEC->GetBinContent(i)-track_multMEC->GetBinContent(i)*0.5);
diffCCQE.push_back(track_multQE->GetBinContent(i)-track_multQE->GetBinContent(i)*0.8);
diffCCRES.push_back(track_multRES->GetBinContent(i)-track_multRES->GetBinContent(i)*0.8);
}
std::vector<double> percentSyst;
std::vector<double> percentStat;
for (int i=0; i<track_mult->GetNbinsX(); i++){
double errRock=diffRock.at(i)*diffRock.at(i);
double errTrkEff=diffTrkEff.at(i)*diffTrkEff.at(i);
double errCCQE=diffCCQE.at(i)*diffCCQE.at(i);
double errCCMEC=diffCCMEC.at(i)*diffCCMEC.at(i);
double errCCRES=diffCCRES.at(i)*diffCCRES.at(i);
double errCCDIS=diffCCDIS.at(i)*diffCCDIS.at(i);
double errStat=track_mult->GetBinError(i)*track_mult->GetBinError(i);
percentStat.push_back(TMath::Sqrt(errStat)/track_mult->GetBinContent(i));
//std::cout<<errRock<<","<<errTrkEff<<","<<errCCQE<<","<<errCCMEC<<","<<errCCRES<<","<<errCCDIS<<","<<errStat<<std::endl;
double systOnly=TMath::Sqrt(errRock+errTrkEff+errCCQE+errCCMEC+errCCRES+errCCDIS);
track_mult->SetBinError(i,TMath::Sqrt(errStat+errRock+errTrkEff+errCCQE+errCCMEC+errCCRES+errCCDIS));
std::cout<<i-1<<" & "<<track_mult->GetBinContent(i)<<" & "<<track_mult->GetBinError(i)<<" & "<<TMath::Sqrt(errStat)<<" & "<<systOnly<<std::endl;
percentSyst.push_back(systOnly/track_mult->GetBinContent(i));
}

TH1D* relativeError=new TH1D("relativeError","relativeError",7,0,7);
TH1D* relativeErrorStat=new TH1D("relativeErrorStat","relativeErrorStat",7,0,7);
for (int i=2; i<track_mult->GetNbinsX()-9; i++){

relativeError->SetBinContent(i,percentSyst.at(i));
relativeErrorStat->SetBinContent(i,percentStat.at(i));
}
relativeErrorStat->GetYaxis()->SetTitleOffset(1.3);
relativeErrorStat->GetXaxis()->SetRangeUser(1,7);
relativeError->SetLineColor(kBlue); relativeError->SetMarkerColor(kBlue);
relativeError->SetTitle("Relative Error");
relativeError->GetYaxis()->SetTitle("Fractional Error");
relativeError->GetXaxis()->SetTitle("Track Multiplicity");
relativeError->Draw("HIST");
relativeErrorStat->Draw("HIST SAME");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");

c1.Print("relativeErrorMR5.png"); c1.Print("relativeErrorMR5.pdf");





track_mult->SetLineColor(kBlack);
track_mult->SetMarkerColor(kBlack);
track_mult->GetYaxis()->SetRangeUser(0,800);
track_mult->SetTitle("Selected Int.: Prelim. Errors");
track_mult->Draw("E0 HIST");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("track_multPrelimErrMR5.png"); c1.Print("track_multPrelimErrMR5.pdf");
track_mult->Scale(5.71);
rwTracks->Scale(5.71);
for (int i=2; i<11; i++){
double syst=TMath::Power(track_mult->GetBinContent(i)*percentSyst.at(i),2);
double stat=track_mult->GetBinContent(i);
track_mult->SetBinError(i,TMath::Sqrt(stat));
std::cout<<i-1<<" & "<<track_mult->GetBinContent(i)<<" & "<<track_mult->GetBinError(i)<<" & "<<TMath::Sqrt(stat)<<" & "<<TMath::Sqrt(syst)<<std::endl;

relativeError->SetBinContent(i,percentSyst.at(i));
relativeErrorStat->SetBinContent(i,1.0/TMath::Sqrt(stat));
}
relativeError->SetLineColor(kBlue); relativeError->SetMarkerColor(kBlue);
relativeErrorStat->SetTitle("Relative Error: Stat. Err.");
relativeErrorStat->GetYaxis()->SetTitle("Fractional Error");
relativeErrorStat->GetXaxis()->SetTitle("Track Multiplicity");
relativeErrorStat->Draw("HIST");
//relativeErrorStat->Draw("HIST SAME");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("relativeError20daysMR5.png"); c1.Print("relativeError20daysMR5.pdf");


for (int i=2; i<11; i++){
double syst=TMath::Power(track_mult->GetBinContent(i)*percentSyst.at(i),2);
double stat=rwTracks->GetBinContent(i);
relativeErrorStat->SetBinContent(i,1.0/TMath::Sqrt(stat));
}
relativeErrorStat->GetXaxis()->SetRangeUser(1,6);
relativeError->SetLineColor(kBlue); relativeError->SetMarkerColor(kBlue);
relativeErrorStat->SetTitle("Data Est. Eff.: Stat. Err.");
relativeErrorStat->GetYaxis()->SetTitle("Fractional Error");
relativeErrorStat->GetXaxis()->SetTitle("Track Multiplicity");
relativeErrorStat->Draw("HIST");
//relativeErrorStat->Draw("HIST SAME");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
   c1.Print("relativeError20daysBadTrkEffMR5.png"); c1.Print("relativeError20daysBadTrkEffMR5.pdf");






std::cout<<track_multQE->GetEntries()<<","<<track_multMEC->GetEntries()<<","<<track_multRES->GetEntries()<<","<<track_multDIS->GetEntries()<<std::endl;
track_mult->GetYaxis()->SetTitleOffset(1.4);
track_mult->GetYaxis()->SetRangeUser(0,5000);
track_mult->SetTitle("Sel. Int.: 20 Days (Stat. Err.)");
track_mult->Draw("E0 HIST");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");

c1.Print("track_mult20daysStatMR5.png"); c1.Print("track_mult20daysStatMR5.pdf");





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
c1.Print("rockWithoutMNVMR5.png"); c1.Print("rockWithoutMNVMR5.pdf");


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
c1.Print("track_multWithoutMNVMR5.png"); c1.Print("track_multWithoutMNVMR5.pdf");


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
dotProduct->GetXaxis()->SetTitle("cos(#Delta #theta)");
dotProduct->GetYaxis()->SetTitle("Number of Tracks");
dotProduct->SetTitle("DS MINERvA-TPC Matching");



TH1D* deltaXUS=(TH1D*)fMinerva.Get("deltaXUS");
TH1D* deltaYUS=(TH1D*)fMinerva.Get("deltaYUS");
deltaXUS->GetXaxis()->CenterTitle();
deltaXUS->GetXaxis()->SetTitle("DS #Delta X");
deltaXUS->GetYaxis()->SetTitle("Number of Tracks");
deltaXUS->SetTitle("Through-Going MINERvA-TPC Matching");
deltaYUS->SetTitle("Through-Going MINERvA-TPC Matching");
deltaXUS->SetLineColor(kRed);
deltaXUS->SetMarkerColor(kRed);
deltaXUS->GetYaxis()->CenterTitle();
deltaYUS->GetXaxis()->CenterTitle();
deltaYUS->GetYaxis()->CenterTitle();
deltaYUS->GetXaxis()->SetTitle("DS #Delta Y");
deltaYUS->GetYaxis()->SetTitle("Number of Tracks");
deltaYUS->SetLineColor(kRed);
deltaYUS->SetMarkerColor(kRed);


TH1D* dotProductUS=(TH1D*)fMinerva.Get("dotProductUS");
dotProductUS->GetXaxis()->CenterTitle();
dotProductUS->GetYaxis()->CenterTitle();
dotProductUS->SetLineColor(kRed);
dotProductUS->SetMarkerColor(kRed);
dotProductUS->GetXaxis()->SetTitle("cos(#Delta #theta)");
dotProductUS->GetYaxis()->SetTitle("Number of Tracks");
dotProductUS->SetTitle("Through-Going MINERvA-TPC Matching");



deltaXUS->GetYaxis()->SetTitleOffset(1.4);
dotProductUS->GetYaxis()->SetTitleOffset(1.4);
deltaYUS->GetYaxis()->SetTitleOffset(1.4);

deltaXUS->Draw("E0 P");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("deltaXUSMR5.png"); c1.Print("deltaXUSMR5.pdf");

deltaYUS->Draw("E0 P");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("deltaYUSMR5.png"); c1.Print("deltaYUSMR5.pdf");




dotProductUS->Draw("E0 P");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("dotProductUSMR5.png"); c1.Print("dotProductUSMR5.pdf");







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
c1.Print("deltaXMR5.png"); c1.Print("deltaXMR5.pdf");

deltaY->Draw("E0 P");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("deltaYMR5.png"); c1.Print("deltaYMR5.pdf");


deltaXGood->Draw("E0 P");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("deltaXGoodMR5.png"); c1.Print("deltaXGoodMR5.pdf");


deltaYGood->Draw("E0 P");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("deltaYGoodMR5.png"); c1.Print("deltaYGoodMR5.pdf");



deltaXBad->Draw("E0 P");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("deltaXBadMR5.png"); c1.Print("deltaXBadMR5.pdf");


deltaYBad->Draw("E0 P");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("deltaYBadMR5.png"); c1.Print("deltaYBadMR5.pdf");

dotProduct->GetYaxis()->SetTitleOffset(1.4);
dotProduct->Draw("E0 P");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("dotProductMR5.png"); c1.Print("dotProductMR5.pdf");


dotProductGood->Draw("E0 P");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("dotProductGoodMR5.png"); c1.Print("dotProductGoodMR5.pdf");

track_mult->GetYaxis()->SetRangeUser(0,1000);
track_mult->SetLineWidth(0); track_mult->GetXaxis()->SetRangeUser(0,10);
track_mult->Draw("HIST");
h1d_track->Draw("SAME");
lTrack->Draw("SAME");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("track_multStack20MR5.png"); c1.Print("track_multStack20MR5.pdf");
track_mult->SetLineWidth(track_multGood->GetLineWidth());


track_multGood->Draw("E0 HIST");
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("track_multGoodMR5.png"); c1.Print("track_multGoodMR5.pdf");


track_mult->SetLineWidth(0); track_mult->GetXaxis()->SetRangeUser(0,10); track_mult->SetTitle("Scaled to 1 Day");
track_mult->Scale(24.0/(82.0*5.71));
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
c1.Print("track_multStack1DayMR5.png"); c1.Print("track_multStack1DayMR5.pdf");
track_mult->SetLineWidth(track_multGood->GetLineWidth());

track_mult->SetLineWidth(0); track_mult->GetXaxis()->SetRangeUser(0,10); track_mult->SetTitle("Sel. Int.: 20 Days");
track_mult->Scale(5.71*82.0/24.0);
track_multQE->Scale(5.71*82.0/24.0);
track_multMEC->Scale(5.71*82.0/24.0);
track_multDIS->Scale(5.71*82.0/24.0);
track_multRES->Scale(5.71*82.0/24.0);
track_multCOH->Scale(5.71*82.0/24.0);
track_multNC->Scale(5.71*82.0/24.0);
track_multRock->Scale(5.71*82.0/24.0);
track_multSec->Scale(5.71*82.0/24.0);


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
track_mult->GetXaxis()->SetRangeUser(1,6);
track_mult->GetYaxis()->SetRangeUser(1,1E8);
//track_mult->GetYaxis()->SetRangeUser(0,600*240.0/82.0);
//track_mult->SetTitle("Data Est. Eff.");
//c1.SetLogy();
c1.SetLogy();
track_mult->Draw("HIST");
h1d_track3->Draw("SAME HIST");
lTrack->Draw("SAME");
c1.SetLogy();
    tL.SetNDC();
    tL.DrawLatex(0.20,0.94,"#bf{DUNE:ND-LAr 2x2}");
c1.Print("track_multStack20DaysMR5.png"); c1.Print("track_multStack20DaysMR5.pdf");
track_mult->SetLineWidth(track_multGood->GetLineWidth());



track_mult->Scale(1.0/5.71);
track_multQE->Scale(1.0/5.71);
track_multMEC->Scale(1.0/5.71);
track_multDIS->Scale(1.0/5.71);
track_multRES->Scale(1.0/5.71);
track_multCOH->Scale(1.0/5.71);
track_multNC->Scale(1.0/5.71);
track_multRock->Scale(1.0/5.71);
track_multSec->Scale(1.0/5.71);
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
c1.Print("track_multStackRockTimes2MR5.png"); c1.Print("track_multStackRockTimes2MR5.pdf");
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
c1.Print("rockMR5.png"); c1.Print("rockMR5.pdf");


}

