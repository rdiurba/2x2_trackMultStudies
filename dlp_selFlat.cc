// Macro to read the CAF files and select muon neutrino interactions and save histograms for unfolding and plots.
// Guide
// Loop over truth information at line 290
// Loop over reco at around line 410
// Find the muon that passes through MINERvA at line 560
// Select neutrinos witha muon selected at line 690
// Make effificiency plots for hadrons at around line 800
// Fill track mult. in line 1146
// Fill cosLReco in line 706
#include "TChain.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "duneanaobj/StandardRecord/Proxy/SRProxy.h"
#include "TMath.h"
#include "TVector.h"
#include "TVector3.h"
#include <fstream>
#include <iostream>
#include <string>

double minTrkLength = 3;


int caf_plotter(std::string input_file_list, std::string output_rootfile,
                bool mcOnly) {
  int muons=0; int pions=0;
  double minLength = 0;
  int goodIntNum = 0;
  int badIntNum = 0;
  int goodIntInFidVol = 0;
  int badIntInFidVol = 0;
  int goodMINERvAMatch = 0;
  int totalMINERvAMatch = 0;
  int totalNuAr=0;
  int totalNuArFidVol=0;
  double offset = 0;
  std::vector<int> goodEvents;
  // Set a bunch of histograms and variables

  // Set the beam direction  
  const auto& beam_dir = TVector3(0.0, -0.05836, 1.0); //(-3.34 degrees in Y)
  double beam_x = 0;
  double beam_y = -0.05836 / TMath::Sqrt(0.05836 * 0.05836 + 1);
  double beam_z = 1.0 / TMath::Sqrt(0.05836 * 0.05836 + 1);

  // Define histograms
  TH1D *part_energy_hist =
      new TH1D("recpart_energy", "Reco particle energy in GeV", 100, 0, 1);
  TH1D *trueEnu =
      new TH1D("trueEnu", "trueEnu", 100, 0, 20);

  int binsMult = 7;
  Double_t edgesMult[8] = {1, 2, 3, 4, 5, 6, 8, 10};
  TH1D *track_length = new TH1D("track_length", "track_length", 100, 0, 100);
  TH1D *part_multTrkOnly =
      new TH1D("part_multTrkOnly", "part_multTrkOnly", binsMult, edgesMult);
  TH1D *part_mult = new TH1D("part_mult", "part_mult", binsMult, edgesMult);
  TH1D *bad_origin = new TH1D("rock", "rock", 2, 0, 2);
  TH1D *track_mult = new TH1D("track_mult", "track_mult", binsMult, edgesMult);
  TH1D *track_multMode =
      new TH1D("track_multMode", "track_multMode", binsMult, edgesMult);
  TH1D *track_multQE =
      new TH1D("track_multQE", "track_multQE", binsMult, edgesMult);
  TH1D *track_multMEC =
      new TH1D("track_multMEC", "track_multMEC", binsMult, edgesMult);
  TH1D *track_multRES =
      new TH1D("track_multRES", "track_multRES", binsMult, edgesMult);
  TH1D *track_multCOH =
      new TH1D("track_multCOH", "track_multCOH", binsMult, edgesMult);
  TH1D *track_multDIS =
      new TH1D("track_multDIS", "track_multDIS", binsMult, edgesMult);
  TH1D *track_multNC =
      new TH1D("track_multNC", "track_multNC", binsMult, edgesMult);
  TH1D *track_multRock =
      new TH1D("track_multRock", "track_multRock", binsMult, edgesMult);
  TH1D *track_multSec =
      new TH1D("track_multSec", "track_multSec", binsMult, edgesMult);

  TH1D *track_multBad =
      new TH1D("track_multBad", "track_multBad", binsMult, edgesMult);
  TH1D *track_multGood =
      new TH1D("track_multGood", "track_multGood", binsMult, edgesMult);
  TH1D *trackCorrectness =
      new TH1D("track_correctness", "track_correctness", 4, 0, 4);
  TH1D *showerCorrectness =
      new TH1D("shower_corrrectness", "shower_correctness", 4, 0, 4);
  TH1D *recoHistVertexY =
      new TH1D("recoHistVertexY", "recoHistVertexY", 140, -70, 70);
  TH1D *recoHistVertexX =
      new TH1D("recoHistVertexX", "recoHistVertexX", 140, -70, 70);
  TH1D *recoHistVertexZ =
      new TH1D("recoHistVertexZ", "recoHistVertexZ", 140, -70, 70);
  TH1D *diffVertexHist =
      new TH1D("diffVertexHist", "diffVertexHist", 100, 0, 10);

  TH1D *recoBacktrackCCAr =
      new TH1D("recoBacktrackCCAr", "recoBacktrackCCAr", 2, 0, 2);
  TH1D *recoBacktrackCCRock =
      new TH1D("recoBacktrackCCRock", "recoBacktrackCCRock", 2, 0, 2);
  TH1D *recoBacktrackCCSec =
      new TH1D("recoBacktrackCCSec", "recoBacktrackCCSec", 2, 0, 2);

  TH2D *responseMult = new TH2D("responseMult", "responseMult", binsMult,
                                edgesMult, binsMult, edgesMult);
  TH1D *backtrackTrueMult =
      new TH1D("backtrackTrueMult", "backtrackTrueMult", binsMult, edgesMult);
  TH1D *recoBacktrackPDGAr =
      new TH1D("recoBacktrackPDGAr", "recoBacktrackPDGAr", 40, -20, 20);
  TH1D *recoBacktrackPDGRock =
      new TH1D("recoBacktrackPDGRock", "recoBacktrackPDGRock", 40, -20, 20);
  TH1D *recoBacktrackPDGSec =
      new TH1D("recoBacktrackPDGSec", "recoBacktrackPDGSec", 40, -20, 20);

  TH1D *recoBacktrackElAr =
      new TH1D("recoBacktrackElAr", "recoBacktrackElAr", 40, 0, 20);

  Double_t edges[7] = {0.91, 0.96, 0.98, 0.9887, 0.994, 0.9974, 1};
    TH1D *trueBacktrackedCosL_zoomOut = new TH1D(
      "trueBacktrackedCosL_zoomOut", "trueBacktrackedCosL_zoomOut", 50, 0.8, 1);
  TH1D *trueCosL_zoomOut = new TH1D(
      "trueCosL_zoomOut", "trueCosL_zoomOut", 50, 0.8, 1);
  TH1D *trueCosLNumubar_zoomOut = new TH1D(
      "trueCosLNumubar_zoomOut", "trueCosLNumubar_zoomOut", 50, 0.8, 1);
      TH1D *trueCosLNumu_zoomOut = new TH1D(
      "trueCosLNumu_zoomOut", "trueCosLNumu_zoomOut", 50, 0.8, 1);
  TH1D *trueBacktrackedCosL =
      new TH1D("trueBacktrackedCosL", "trueBacktrackedCosL", 6, edges);
  TH1D *minervaMatchPDG = new TH1D("minervaMatchPDG", "minervaMatchPDG", 6000, -3000, 3000);
  TH1D *minervaMatchE = new TH1D("minervaMatchE", "minervaMatchE", 50, 0, 20);
  TH1D *minervaMatchCos = new TH1D("minervaMatchCos", "minervaMatchCos", 100, -1, 1);
  TH1D *selPionE = new TH1D("selPionE", "selPionE", 15, 0, 0.2);
  TH1D *selProtonE = new TH1D("selProtonE", "selProtonE", 15, 0, 0.2);

  TH1D *selPionDirX = new TH1D("selPionDirX", "selPionDirX", 20, -1, 1);
  TH1D *selPionDirY = new TH1D("selPionDirY", "selPionDirY", 20, -1, 1);
  TH1D *selPionDirZ = new TH1D("selPionDirZ", "selPionDirZ", 20, -1, 1);

  TH1D *selProtonDirX = new TH1D("selProtonDirX", "selProtonDirX", 20, -1, 1);
  TH1D *selProtonDirY = new TH1D("selProtonDirY", "selProtonDirY", 20, -1, 1);
  TH1D *selProtonDirZ = new TH1D("selProtonDirZ", "selProtonDirZ", 20, -1, 1);

  TH1D *selTrkLen = new TH1D("selTrkLen", "selTrkLen", 20, 0, 10);
  TH1D *selTrkLenNonQE = new TH1D("selTrkLenNonQE", "selTrkLenNonQE", 20, 0, 10);
  TH1D *selTrkLenQE = new TH1D("selTrkLenQE", "selTrkLenQE", 20, 0, 10);
  TH1D *selProtonXZ = new TH1D("selProtonXZ", "selProtonXZ", 20, -1, 1);
  TH1D *selProtonYZ = new TH1D("selProtonYZ", "selProtonYZ", 20, -1, 1);
  TH1D *selProtonXY = new TH1D("selProtonXY", "selProtonXY", 20, -1, 1);

  TH1D *selPionXZ = new TH1D("selPionXZ", "selPionXZ", 20, -1, 1);
  TH1D *selPionYZ = new TH1D("selPionYZ", "selPionYZ", 20, -1, 1);
  TH1D *selPionXY = new TH1D("selPionXY", "selPionXY", 20, -1, 1);


  TH2D *confusionMatrix =  new TH2D("confusionMatrix", "confusionMatrix", 4, 0, 4, 4, 0, 4);

  TH1D *nPi = new TH1D("nPi", "nPi", binsMult, edges);
  TH1D *nP = new TH1D("nP", "nP", 20, 0, 20);
  TH1D *escapePi = new TH1D("escapePi", "escapePi", 20, 0, 20);
  TH1D *containPiLen = new TH1D("containPiLen", "containPiLen", 100, 0, 20);
  TH1D *containPLen = new TH1D("containP", "containP", 100, 0, 20);

  TH1D *truePionEWithRecoInt =  new TH1D("truePionEWithRecoInt", "truePionEWithRecoInt", 15, 0, 0.2);
  TH1D *trueProtonEWithRecoInt =  new TH1D("trueProtonEWithRecoInt", "trueProtonEWithRecoInt", 15, 0, 0.2);

  TH1D *trueTrkLenQE = new TH1D("trueTrkLenQE", "trueTrkLenQE", 20, 0, 10);
  TH1D *trueTrkLenNonQE =    new TH1D("trueTrkLenNonQE", "trueTrkLenNonQE", 20, 0, 10);
  TH1D *trueTrkLen = new TH1D("trueTrkLen", "trueTrkLen", 20, 0, 10);
  TH1D *trueProtonWithRecoIntDirX = new TH1D(
      "trueProtonWithRecoIntDirX", "trueProtonWithRecoIntDirX", 20, -1, 1);
  TH1D *trueProtonWithRecoIntDirY = new TH1D(
      "trueProtonWithRecoIntDirY", "trueProtonWithRecoIntDirY", 20, -1, 1);
  TH1D *trueProtonWithRecoIntDirZ = new TH1D(
      "trueProtonWithRecoIntDirZ", "trueProtonWithRecoIntDirZ", 20, -1, 1);

  TH1D *truePionWithRecoIntDirX =
      new TH1D("truePionWithRecoIntDirX", "truePionWithRecoIntDirX", 20, -1, 1);
  TH1D *truePionWithRecoIntDirY =
      new TH1D("truePionWithRecoIntDirY", "truePionWithRecoIntDirY", 20, -1, 1);
  TH1D *truePionWithRecoIntDirZ =
      new TH1D("truePionWithRecoIntDirZ", "truePionWithRecoIntDirZ", 20, -1, 1);



  TH1D *trueProtonWithRecoIntXY =
      new TH1D("trueProtonWithRecoIntXY", "trueProtonWithRecoIntXY", 20, -1, 1);
  TH1D *trueProtonWithRecoIntXZ =
      new TH1D("trueProtonWithRecoIntXZ", "trueProtonWithRecoIntXZ", 20, -1, 1);
  TH1D *trueProtonWithRecoIntYZ =
      new TH1D("trueProtonWithRecoIntYZ", "trueProtonWithRecoIntYZ", 20, -1, 1);

  TH1D *truePionWithRecoIntXY =
      new TH1D("truePionWithRecoIntXY", "truePionWithRecoIntXY", 20, -1, 1);
  TH1D *truePionWithRecoIntXZ =
      new TH1D("truePionWithRecoIntXZ", "truePionWithRecoIntXZ", 20, -1, 1);
  TH1D *truePionWithRecoIntYZ =
      new TH1D("truePionWithRecoIntYZ", "truePionWithRecoIntYZ", 20, -1, 1);



  TH1D *true_mult = new TH1D("true_mult", "true_mult", 15, 0, 15);
  TH1D *true_multTrkOnly =
      new TH1D("true_multTrkOnly", "true_multTrkOnly", binsMult, edgesMult);
  TH1D *trueDiffPosvsPDirZ =
      new TH1D("trueDiffPosvsPDirZ", "trueDiffPosvsPDirZ", 20, -1, 1);

  TH1D *true_multGENIE =      new TH1D("true_multGENIE", "true_multGENIE", binsMult, edgesMult);
  TH2D *responseGenieToG4 = new TH2D("responseGenieToG4", "responseGenieToG4",
                                     binsMult, edges, binsMult, edgesMult);
  TH1D *matchTrue_mult =
      new TH1D("matchTrue_mult", "matchTrue_mult", binsMult, edgesMult);
  TH1D *matchTrue_multTrkOnly = new TH1D(
      "matchTrue_multTrkOnly", "matchTrue_multTrkOnly", binsMult, edgesMult);
  TH1D *matchTrue_multGENIE = new TH1D(
      "matchTrue_multGENIE", "matchTrue_multGENIE", binsMult, edgesMult);
  TH1D *matchHistEl = new TH1D("matchHistEl", "matchHistEl", 100, 0, 20);
  TH1D *matchHistCosl = new TH1D("matchHistCosl", "matchHistCosl", 20, -1, 1);

  TH1D *histEl = new TH1D("histEl", "histEl", 50, 0, 20);

  TH2D *recoVertex2DNoCuts = new TH2D(
      "recoVertex2DNoCuts", "recoVertex2DNoCuts", 70, -70, 70, 70, -70, 70);
  TH2D *recoVertex2D =
      new TH2D("recoVertex2D", "recoVertex2D", 60, -60, 60, 60, -60, 60);
  TH2D *recoVertex2DBadYZ = new TH2D("recoVertex2DBadYZ", "recoVertex2DBadYZ",
                                     70, -70, 70, 70, -70, 70);

  TH1D *recoCosL = new TH1D("recoCosL", "recoCosL", 6, edges);
  TH1D *badRecoCosL = new TH1D("badRecoCosL", "badRecoCosL", 6, edges);
  TH1D *trueBacktrackedCosL_unbinned =
      new TH1D("trueBacktrackedCosL_unbinned", "trueBacktrackedCosL_unbinned",
               50, 0.8, 1);
  TH1D *trueCosL_unbinned =
      new TH1D("trueCosL_unbinned", "trueCosL_unbinned", 50, 0.8, 1);
  TH1D *trueBacktrackedMult_unbinned =
      new TH1D("trueBacktrackedMult_unbinned", "trueBacktrackedMult_unbinned",
               20, 0.0, 20);
  TH1D *trueMult_unbinned =
      new TH1D("trueMult_unbinned", "trueMult_unbinned", 20, 0, 20);
  TH1D *trueCosL = new TH1D("trueCosL", "trueCosL", 6, edges);
  TH2D *responseCosL =
      new TH2D("responseCosL", "responseCosL", 6, edges, 6, edges);

  TH2D *recoVertex2DBad =
      new TH2D("recoVertex2DBad", "recoVertex2DBad", 70, -70, 70, 70, -70, 70);
  TH2D *trueVertex2DBad =
      new TH2D("trueVertex2DBad", "trueVertex2DBad", 60, -60, 60, 60, -60, 60);

  TH1D *longest_track = new TH1D("longest_track", "longest_track", 100, 0, 100);
  TH1D *totalPOT = new TH1D("totalPOT", "totalPOT", 1, 0, 1);
  TH1D *totalSpills = new TH1D("totalSpilles", "totalSpills", 1, 0, 1);
  TH1D *trueIntPerSpill = new TH1D("trueIntPerSpill", "trueIntPerSpill", 10,0,10);
  double sumPOT = 0;


  // Give an input list
  std::ifstream caf_list(input_file_list.c_str());
  int totalEntries=0;
  // Check if input list is present
  if (!caf_list.is_open()) {
    std::cerr << Form("File %s not found", input_file_list.c_str())
              << std::endl;
    return 1;
  }

  // Add files to CAF chain from input list
  std::string tmp;


  while(caf_list >> tmp){

    TFile * f = new TFile(tmp.c_str());
    TTree * caf_chain = (TTree*) f->Get("cafTree");  

  //Check if CAF tree is present
  long Nentries = caf_chain->GetEntries();

  //Define Standard Record and link it to the CAF tree branch "rec"
  const auto& sr = new caf::SRProxy(caf_chain, "rec");
  totalEntries=totalEntries+Nentries;
  totalSpills->Fill(0.5, Nentries);
  bool skipEvent = false;
  // Run over the events

  for (long n = 0; n < Nentries; ++n) {

    double trueVtxX = -9999;
    double trueVtxY = -999;
    double trueVtxZ = -9999;
    skipEvent = false;

    caf_chain->GetEntry(n); // Get spill from tree
    sumPOT = sr->beam.pulsepot / 1e13 + sumPOT;
    int intAboveThresh=0;

    bool hasANeutrino = false;
    double mnvOffsetX = -10;
    double mnvOffsetY = 5;
    if (mcOnly) {
      mnvOffsetX = 0;
      mnvOffsetY = 0;
      // Save truth distributions
      for (long unsigned ntrue = 0; ntrue < sr->mc.nu.size(); ntrue++) {
        const auto& vertex = sr->mc.nu[ntrue].vtx;
        trueVtxX = vertex.x;
        trueVtxY = vertex.y;
        trueVtxZ = vertex.z;

        const auto& truePrimary = sr->mc.nu[ntrue].prim;
        const auto& trueSecondary = sr->mc.nu[ntrue].sec;
        int aboveThresh=0;
        for (long unsigned int k = 0; k < sr->mc.nu[ntrue].prim.size(); k++) {
        int pdg=abs(sr->mc.nu[ntrue].prim[k].pdg);
        const auto& start_pos = sr->mc.nu[ntrue].prim[k].start_pos;
        if (abs(start_pos.x) > 59 || abs(start_pos.x) < 5)
          continue;
        if (abs(start_pos.z) > 59.5 ||
            abs(start_pos.z) < 5)
          continue;
        if (abs(start_pos.y) > 57)
          continue;
        if (pdg!=211 && pdg!=2212 && pdg!=111 && pdg!=13 /*&& pdg!=11*/) continue; 
        const auto& p = sr->mc.nu[ntrue].prim[k].p;
            
        double momentum=TMath::Sqrt(p.px * p.px + p.py * p.py + p.pz * p.pz);
        double energy=sr->mc.nu[ntrue].prim[k].p.E;
        double mass=TMath::Sqrt(energy*energy-momentum*momentum);
        double KE=energy-mass;
        if (KE>0.01) aboveThresh=1;
        }
        for (long unsigned int k = 0; k < sr->mc.nu[ntrue].sec.size(); k++) {
        int pdg=abs(sr->mc.nu[ntrue].sec[k].pdg);
        const auto& start_pos = sr->mc.nu[ntrue].sec[k].start_pos;
        if (abs(start_pos.x) > 59 || abs(start_pos.x) < 5)
          continue;
        if (abs(start_pos.z) > 59.5 ||
            abs(start_pos.z) < 5)
          continue;
        if (abs(start_pos.y) > 57)
          continue;
        if (pdg!=211 && pdg!=2212 && pdg!=111 && pdg!=13 /*&& pdg!=11*/) continue; 
        const auto& p = sr->mc.nu[ntrue].sec[k].p;
            
        double momentum=TMath::Sqrt(p.px * p.px + p.py * p.py + p.pz * p.pz);
        double energy=sr->mc.nu[ntrue].sec[k].p.E;
        double mass=TMath::Sqrt(energy*energy-momentum*momentum);
        double KE=energy-mass;
        if (KE>0.01) aboveThresh=1;
        }
        if (aboveThresh>0) intAboveThresh++;

          
        int truePart = 0;
        int truePartTrkOnly = 0;
        int truePartNoG4 = 0;
        // We only want CC numu interactions on argon in the fiducial volume
        if (sr->mc.nu[ntrue].targetPDG != 1000180400)
          continue;
        if (abs(sr->mc.nu[ntrue].vtx.x) > 63 || abs(sr->mc.nu[ntrue].vtx.x) < 2)
          continue;
        if (abs(sr->mc.nu[ntrue].vtx.z) > 64 ||
            abs(sr->mc.nu[ntrue].vtx.z) < 2)
          continue;
        if (abs(sr->mc.nu[ntrue].vtx.y) > 63)
          continue;
        totalNuAr++;



        if (abs(sr->mc.nu[ntrue].pdg) != 14)
          continue;
        if (abs(sr->mc.nu[ntrue].vtx.x) > 59 || abs(sr->mc.nu[ntrue].vtx.x) < 5)
          continue;
        if (abs(sr->mc.nu[ntrue].vtx.z) > 59.5 ||
            abs(sr->mc.nu[ntrue].vtx.z) < 5)
          continue;
        if (abs(sr->mc.nu[ntrue].vtx.y) > 57)
          continue;
        totalNuArFidVol++;
        if (sr->mc.nu[ntrue].iscc == false)
          continue;
        // Sum the number of hadrons in case you need it
        int npipPrimaries = 0;
        int npimPrimaries = 0;
        int npPrimaries = 0;
        int nneutronsPrimaries = 0;
        int npi0Primaries = 0; 
        for (long unsigned int k = 0; k < sr->mc.nu[ntrue].prim.size(); k++) {
          if (sr->mc.nu[ntrue].prim[k].pdg == 211)
            npipPrimaries++;
          if (sr->mc.nu[ntrue].prim[k].pdg == 2212)
            npPrimaries++;
          if (sr->mc.nu[ntrue].prim[k].pdg == -211)
            npimPrimaries++;
          if (sr->mc.nu[ntrue].prim[k].pdg == 111)
              npi0Primaries++;
          if (sr->mc.nu[ntrue].prim[k].pdg == 2112)
            nneutronsPrimaries++;
        }

        int nProton = 0;
        int nPion = 0;
        int escapingPions = 0;
        double Elep = -999;
        double cosL = -999;
        for (long unsigned primaries = 0; primaries < truePrimary.size();
             primaries++) {
          int pdg = truePrimary[primaries].pdg;
          double px = truePrimary[primaries].p.px;
          double py = truePrimary[primaries].p.py;
          double pz = truePrimary[primaries].p.pz;

          if ((abs(pdg) == 13 || abs(pdg) == 2212 || abs(pdg) == 211 ||
               abs(pdg) == 321)) {
            truePartNoG4++;
            if (abs(pdg) == 211)
              nPion++;
            if (abs(pdg) == 2212)
              nProton++;
            const auto& start_pos = sr->mc.nu[ntrue].prim[primaries].start_pos;
            const auto& end_pos = sr->mc.nu[ntrue].prim[primaries].end_pos;
            // if (std::isnan(start_pos.z)){  continue;}
            const auto& p = sr->mc.nu[ntrue].prim[primaries].p;
            double dX = (end_pos.x - start_pos.x);
            double dY = (end_pos.y - start_pos.y);
            double dZ = (end_pos.z - start_pos.z);
            double dPX = p.px;
            double dPY = p.py;
            double dPZ = p.pz;
            double length = TMath::Sqrt(dX * dX + dY * dY + dZ * dZ);
            double lengthP =
                TMath::Sqrt(p.px * p.px + p.py * p.py + p.pz * p.pz);
            // Save the muon information for the CCINC analysis
            if (sr->mc.nu[ntrue].iscc && abs(sr->mc.nu[ntrue].pdg) == 14 &&
                abs(pdg) == 13) {

              Elep = sr->mc.nu[ntrue].prim[primaries].p.E;
              double dir_x = dPX / lengthP;
              double dir_y = dPY / lengthP;
              double dir_z = dPZ / lengthP;
              cosL = dir_x * beam_x + dir_y * beam_y + dir_z * beam_z;
              histEl->Fill(Elep);
              trueCosL_unbinned->Fill(cosL); trueCosL_zoomOut->Fill(cosL);
                if ((sr->mc.nu[ntrue].pdg) == 14){ trueCosLNumu_zoomOut->Fill(cosL);}
                else{ trueCosLNumubar_zoomOut->Fill(cosL);}
              if (cosL > 0.91 && Elep>1) {
                trueCosL->Fill(cosL);
                trueEnu->Fill(sr->mc.nu[ntrue].E);
              }
              // std::cout<<n<<","<<ntrue<<","<<primaries<<","<<dX<<","<<dY<<","<<dZ<<","<<p.pz/lengthP<<std::endl;

              trueDiffPosvsPDirZ->Fill(cosL - p.pz / lengthP);
            }
            // Some information on if the pion stays or leaves
            if ((abs(end_pos.z) > 62 || abs(end_pos.x) > 62) && abs(pdg) == 211)
              escapingPions++;
            else if (abs(pdg) == 211)
              containPiLen->Fill(length);
            if (abs(end_pos.z) < 62 && abs(end_pos.x) < 62 && abs(pdg) == 2212)
              containPLen->Fill(length);
              
            // Save the track if it is above the minimum threshold
            if (length > minTrkLength) {
              truePartTrkOnly++;
            }
          }
          truePart++;
        }
        // Save hadron information if it is above the cosL signal definition
        if (cosL < 0.9 || Elep<1)
          continue;
        escapePi->Fill(escapingPions);
        nPi->Fill(nPion);
        nP->Fill(nProton);
        true_mult->Fill(truePart);
        true_multTrkOnly->Fill(truePartTrkOnly);
        true_multGENIE->Fill(truePartNoG4);
        responseGenieToG4->Fill(truePart, truePartTrkOnly);
        trueMult_unbinned->Fill(truePartTrkOnly);
        hasANeutrino = true;
      }
    }
    int rock = 0;
    bool goodInteraction = false;
    std::vector<int> trueInteractionIndex;
    std::vector<int> primaryTrkIndex;
    trueIntPerSpill->Fill(intAboveThresh);
    for (long unsigned nixn = 0; nixn < sr->common.ixn.dlp.size(); nixn++) {
      double cosL = -999;
      double endZTemp=-999;
      bool oneContained = false;
      bool oneNotContained = false;
      goodInteraction = false;
      double biggestMatch = -999;
      int biggestMatchIndex = -999;
      double maxDotProductDS = -999;
      double maxDotProductUS = -999;
      int maxEventPar = -999;
      int maxEventTyp = -9999;
      int maxEventIxn = -999;
      double dirZExiting = -999;
      double endZMuonCand=-999;
      double startZMuonCand = -999;
      double dirXExiting = -9999;
      double dirYExiting = -9999;
      std::vector<int> mx2IntCandidateVector, mx2IdxCandidateVector, ndlarTrkCandidateVector;

      double recoVertexX = sr->common.ixn.dlp[nixn].vtx.x;
      double recoVertexY = sr->common.ixn.dlp[nixn].vtx.y;
      double recoVertexZ = sr->common.ixn.dlp[nixn].vtx.z;

      recoVertex2DNoCuts->Fill(sr->common.ixn.dlp[nixn].vtx.x,
                               sr->common.ixn.dlp[nixn].vtx.z);
      if (/*abs(abs(sr->common.ixn.dlp[nixn].vtx.x)-33)<1 || */ abs(
              sr->common.ixn.dlp[nixn].vtx.x) > 59 ||
          abs(sr->common.ixn.dlp[nixn].vtx.x) < 5 ||
          abs(sr->common.ixn.dlp[nixn].vtx.y) > 57 ||
          abs(sr->common.ixn.dlp[nixn].vtx.z) < 5 ||
          abs(sr->common.ixn.dlp[nixn].vtx.z) > 59.5)
        continue;

      if (mcOnly) {

        for (long unsigned int ntruth = 0; ntruth < sr->common.ixn.dlp[nixn].truth.size();
             ntruth++) {

          if (biggestMatch < sr->common.ixn.dlp[nixn].truthOverlap.at(ntruth)) {
            biggestMatch = sr->common.ixn.dlp[nixn].truthOverlap.at(ntruth);
            biggestMatchIndex = sr->common.ixn.dlp[nixn].truth.at(ntruth);
          }
        }
        if (biggestMatchIndex>-1 && sr->mc.nu.size()>0){
        if (sr->mc.nu[biggestMatchIndex].id > 1E9) {
          rock = 1;
        }
        trueVtxX = sr->mc.nu[biggestMatchIndex].vtx.x;
        trueVtxY = sr->mc.nu[biggestMatchIndex].vtx.y;
        trueVtxZ = sr->mc.nu[biggestMatchIndex].vtx.z;
        double diffVtx =
            TMath::Sqrt((recoVertexY - trueVtxY) * (recoVertexY - trueVtxY) +
                        (recoVertexX - trueVtxX) * (recoVertexX - trueVtxX) +
                        (recoVertexZ - trueVtxZ) * (recoVertexZ - trueVtxZ));
        if (abs(sr->mc.nu[biggestMatchIndex].pdg) == 14 &&
            abs(sr->mc.nu[biggestMatchIndex].iscc) == true &&
            sr->mc.nu[biggestMatchIndex].targetPDG == 1000180400 &&
            diffVtx < 5 && abs(trueVtxX) < 59 && abs(trueVtxZ) > 5 &&
            abs(trueVtxZ) < 59.5 && abs(trueVtxY) < 57){

                        for (int primaries = 0;
                 primaries < sr->mc.nu[biggestMatchIndex].prim.size();
                 primaries++) {
            if (abs(sr->mc.nu[biggestMatchIndex].prim[primaries].pdg) != 13)
                continue;
              double Elep = sr->mc.nu[biggestMatchIndex].prim[primaries].p.E;
              const auto& start_pos =
                  sr->mc.nu[biggestMatchIndex].prim[primaries].start_pos;
              const auto& end_pos =
                  sr->mc.nu[biggestMatchIndex].prim[primaries].end_pos;
              const auto& p = sr->mc.nu[biggestMatchIndex].prim[primaries].p;
              // if (std::isnan(start_pos.z)) continue;
              double dX = (p.px);
              double dY = (p.py);
              double dZ = (p.pz);
              double length =
                  TMath::Sqrt(p.px * p.px + p.py * p.py + p.pz * p.pz);
              double dir_x = dX / length;
              double dir_y = dY / length;
              double dir_z = dZ / length;
              double cosLPrelim = dir_x * beam_x + dir_y * beam_y + dir_z * beam_z;

             if (cosLPrelim<0.91 || Elep<1) goodInteraction = false; else goodInteraction=true;}
        }}
        recoHistVertexY->Fill(sr->common.ixn.dlp[nixn].vtx.y);
        recoHistVertexX->Fill(sr->common.ixn.dlp[nixn].vtx.x);
        recoHistVertexZ->Fill(sr->common.ixn.dlp[nixn].vtx.z);
        if (goodInteraction == true)
          goodIntNum++;
        else
          badIntNum++;



        trueInteractionIndex.push_back(biggestMatchIndex);

        if (goodInteraction == true)
          goodIntInFidVol++;
        else
          badIntInFidVol++;
      }
      int partMult = 0;
      int partMultTrkOnly = 0;
      bool hasAtLeastOneProton = false;
      bool hasAtLeastOneMuon = false;
      bool hasAtLeastOne = false;
      // std::cout<<"Interaction: "<<nixn<<std::endl;
      int correctTrack = 2;
      int correctShower = 2;
      double longestTrk = -9999;
      int trackMult = 0;
      int trackMultExit = 0;
      int minervaTracks = 0;
      int minervaThrough = 0;
      int mx2IntCandidate=-9;
      int mx2IdxCandidate=-9;
      for (long unsigned npart = 0;
           npart < sr->common.ixn.dlp[nixn].part.dlp.size();
           npart++) { // loop over particles
      //  if (!sr->common.ixn.dlp[nixn].part.dlp[npart].primary)
      //    continue;
        int pdg = sr->common.ixn.dlp[nixn].part.dlp[npart].pdg;
        // Loop over primary tracks
        if ((abs(pdg) == 2212 || abs(pdg) == 13 || abs(pdg) == 211 ||
             abs(pdg) == 321)) {

          const auto& start_pos = sr->common.ixn.dlp[nixn].part.dlp[npart].start;
          const auto& end_pos = sr->common.ixn.dlp[nixn].part.dlp[npart].end;
          double diffVertexdZ =
              abs(start_pos.z - sr->common.ixn.dlp[nixn].vtx.z);
          double diffVertexdX =
              abs(start_pos.x - sr->common.ixn.dlp[nixn].vtx.x);
          double diffVertexdY =
              abs(start_pos.y - sr->common.ixn.dlp[nixn].vtx.y);
          double diffVertex = TMath::Sqrt(diffVertexdZ * diffVertexdZ +
                                          diffVertexdY * diffVertexdY +
                                          diffVertexdX * diffVertexdX);
          if (diffVertex > 5)
            continue;
          // Make sure it is near the vertex
          double dX = (end_pos.x - start_pos.x);
          double dY = (end_pos.y - start_pos.y);
          double dZ = (end_pos.z - start_pos.z);
          double length = TMath::Sqrt(dX * dX + dY * dY + dZ * dZ);
          double dirX = dX / length;
          double dirY = dY / length;
          double dirZ = dZ / length;


          if (std::isnan(start_pos.z))
            length = -999;
          if (length > longestTrk)
            longestTrk = length;
          // Make sure it is above the track threshold
          if (sr->common.ixn.dlp[nixn].part.dlp[npart].primary == true &&
              length > minTrkLength) {
            partMult++;
            trackMult++;
            // see if it punches out and match it to MINERvA
            if ((start_pos.z) > 62 || (end_pos.z) > 62)
              trackMultExit++;
            int maxPartMinerva = -999;
            int maxTypeMinerva = -999;
            int maxIxnMinerva = -999;
            int maxPartMinervaUS = -999;
            int maxTypeMinervaUS = -999;
            if ((abs(start_pos.z) > 62 || abs(end_pos.z) > 62)) {
              int minervaPass = 0;
              double dotProductDS = -999;
              double deltaExtrapYUS = -999;
              double deltaExtrapY = -999;
              double dotProductUS = -999;
              double deltaExtrapX = -999;
              double deltaExtrapXUS = -999;
              double dotProductFromCAF=-999;

        	 	for(int k=0; k<sr->nd.trkmatch.extrap.size(); k++){
                if (sr->nd.trkmatch.extrap[k].larid.ixn!=nixn || sr->nd.trkmatch.extrap[k].larid.reco!=1) continue;
                    int index=sr->nd.trkmatch.extrap[k].larid.idx;
                    if (sr->nd.lar.dlp[nixn].tracks[index].start.z==sr->common.ixn.dlp[nixn].part.dlp[npart].start.z && sr->nd.lar.dlp[nixn].tracks[index].end.z==sr->common.ixn.dlp[nixn].part.dlp[npart].end.z){
                   double dotProductTemp=abs(sr->nd.trkmatch.extrap[k].angdispl);
                    if (dotProductTemp<0.99) continue;
                 int i=sr->nd.trkmatch.extrap[k].minervaid.ixn;
                 int j=sr->nd.trkmatch.extrap[k].minervaid.idx;


                  double dir_z = sr->nd.minerva.ixn[i].tracks[j].dir.z;
                  double end_z = sr->nd.minerva.ixn[i].tracks[j].end.z;
                  double start_z = sr->nd.minerva.ixn[i].tracks[j].start.z;
                  double end_x = sr->nd.minerva.ixn[i].tracks[j].end.x;
                  double start_x = sr->nd.minerva.ixn[i].tracks[j].start.x;

                  double end_y = sr->nd.minerva.ixn[i].tracks[j].end.y;
                  double start_y = sr->nd.minerva.ixn[i].tracks[j].start.y;
                  // Try to match a track starting in 2x2 to MINERvA
                  if (start_z > 0 && start_z < 170 && end_z > 170 &&
                      ((end_pos.z) > 62)) {
                    int truthPart =
                        sr->nd.minerva.ixn[i].tracks[j].truth[0].part;
                    double dXMnv = (sr->nd.minerva.ixn[i].tracks[j].end.x -
                                    sr->nd.minerva.ixn[i].tracks[j].start.x);
                    double dYMnv = (sr->nd.minerva.ixn[i].tracks[j].end.y -
                                    sr->nd.minerva.ixn[i].tracks[j].start.y);
                    double dZMnv = (sr->nd.minerva.ixn[i].tracks[j].end.z -
                                    sr->nd.minerva.ixn[i].tracks[j].start.z);
                    double lengthMinerva = TMath::Sqrt(
                        dXMnv * dXMnv + dYMnv * dYMnv + dZMnv * dZMnv);
                    if (lengthMinerva < 10)
                      continue;
                    double dirXMinerva = dXMnv / lengthMinerva;
                    double dirYMinerva = dYMnv / lengthMinerva;
                    double dirZMinerva = dZMnv / lengthMinerva;
                    double dotProduct = dirXMinerva * dirX +
                                        dirYMinerva * dirY + dirZ * dirZMinerva;
                    double extrapdZ = start_z - end_pos.z;
                    double extrapY =
                        dirY / dirZ * (extrapdZ) + end_pos.y - start_y;
                    double extrapX =
                        dirX / dirZ * (extrapdZ) + end_pos.x - start_x;
                    double diffExtrap =
                        TMath::Sqrt(TMath::Power(extrapY - start_y, 2));

                    if (dotProductDS < dotProductTemp) {
                      dotProductDS = dotProductTemp;
                      deltaExtrapY = extrapY;
                      deltaExtrapX = extrapX;
                      endZTemp=end_z;
                      mx2IntCandidate=i;
                      mx2IdxCandidate=j;
                      dirZExiting = dirZ;
                      dirXExiting = dirX;
                      dirYExiting = dirY;
                      if (mcOnly) {
                        maxPartMinerva =
                            sr->nd.minerva.ixn[i].tracks[j].truth[0].part;
                        maxTypeMinerva =
                            sr->nd.minerva.ixn[i].tracks[j].truth[0].type;
                        maxIxnMinerva =
                            sr->nd.minerva.ixn[i].tracks[j].truth[0].ixn;
                      }
                      // if (end_z>300){ minervaPass=1;} if(dirZExiting<dirZ){
                      // dirZExiting=dirZ;}
                    }
                  }
                  // Try matching a through-going particle of both to each other
                  if (start_z < 0 && end_z > 0 &&
                      ((start_pos.z < -62 && end_pos.z > 62) ||
                       (start_pos.z > 62 && end_pos.z < -62))) {
                    int truthPart =
                        sr->nd.minerva.ixn[i].tracks[j].truth[0].part;
                    double dXMnv = (sr->nd.minerva.ixn[i].tracks[j].end.x -
                                    sr->nd.minerva.ixn[i].tracks[j].start.x);
                    double dYMnv = (sr->nd.minerva.ixn[i].tracks[j].end.y -
                                    sr->nd.minerva.ixn[i].tracks[j].start.y);
                    double dZMnv = (sr->nd.minerva.ixn[i].tracks[j].end.z -
                                    sr->nd.minerva.ixn[i].tracks[j].start.z);
                    double lengthMinerva = TMath::Sqrt(
                        dXMnv * dXMnv + dYMnv * dYMnv + dZMnv * dZMnv);
                    double dirXMinerva = dXMnv / lengthMinerva;
                    double dirYMinerva = dYMnv / lengthMinerva;
                    double dirZMinerva = dZMnv / lengthMinerva;
                    double dotProduct = dirXMinerva * dirX +
                                        dirYMinerva * dirY + dirZ * dirZMinerva;

                    double extrapdZUS = end_z - end_pos.z;
                    double extrapYUS =
                        dirY / dirZ * (extrapdZUS) + end_pos.y - end_y;
                    double extrapXUS =
                        dirX / dirZ * (extrapdZUS) + end_pos.x - end_x;

                    // double
                    // diffExtrap=TMath::Sqrt(TMath::Power(extrapY-end_y,2));
                    if (dotProductUS < dotProductTemp)
                      dotProductUS = dotProductTemp;
                    deltaExtrapYUS = extrapYUS;
                    deltaExtrapXUS = extrapXUS;
                    if (mcOnly) {
                      maxPartMinervaUS =
                          sr->nd.minerva.ixn[i].tracks[j].truth[0].part;
                      maxTypeMinervaUS =
                          sr->nd.minerva.ixn[i].tracks[j].truth[0].type;
                    }
                  }
                }
              } // Minerva
              // See how good the best one was
              if (dotProductDS > maxDotProductDS) {
                maxDotProductDS = dotProductDS;
                mx2IntCandidateVector.push_back(mx2IntCandidate);
                mx2IdxCandidateVector.push_back(mx2IdxCandidate);
                ndlarTrkCandidateVector.push_back(npart);
                maxEventPar = maxPartMinerva;
                maxEventTyp = maxTypeMinerva;
                maxEventIxn = maxIxnMinerva;
                
              }
              if (dotProductDS > 0.99) {
                minervaTracks++;
                if (minervaPass == 1) {
                  minervaThrough++;
                  startZMuonCand = start_pos.z;
                  if (start_pos.z > end_pos.z)
                    startZMuonCand = end_pos.z;
                }
              }
              if (dotProductUS > maxDotProductUS)
                maxDotProductUS = dotProductUS;
              if (dotProductUS > 0.99) { //could add a cout here to get rid of any overlapping through-going muons
        	 }
            } // exiting particles
              // else oneContained=true;
          }   // primary of particle greater than 5
        }
      } // particles
     if (/*minervaThrough<1  ||*/ maxDotProductDS < 0.99 || mx2IntCandidateVector.size()==0)
        continue;
      // If you got one matched to downstream Mx2, then it likely was the muon
      int muonIndex=-1; double muonCandidateZ=-170;
      for (int k=0; k<mx2IntCandidateVector.size(); k++){
            int i=mx2IntCandidateVector.at(k); int j=mx2IdxCandidateVector.at(k); 
            double end_z = sr->nd.minerva.ixn[i].tracks[j].end.z;
            if (end_z>muonCandidateZ){ muonIndex=k; muonCandidateZ=end_z;}

            }
     // Calculate muon direction based on distance
      int npart=ndlarTrkCandidateVector.at(muonIndex);
      const auto& start_pos = sr->common.ixn.dlp[nixn].part.dlp[npart].start;
      const auto& end_pos = sr->common.ixn.dlp[nixn].part.dlp[npart].end;
      if (mcOnly==1){
                      int i=mx2IntCandidateVector.at(muonIndex); int j=mx2IdxCandidateVector.at(muonIndex); 
                      int  ixnMinerva=sr->nd.minerva.ixn[i].tracks[j].truth[0].ixn;
                         int       partMinerva =
                          sr->nd.minerva.ixn[i].tracks[j].truth[0].part;
                      int typeMinerva =
                          sr->nd.minerva.ixn[i].tracks[j].truth[0].type;
          if (typeMinerva==1){
        //std::cout<<sr->mc.nu[maxEventIxn].prim[maxEventPar].pdg<<std::endl;
        if (abs(sr->mc.nu[ixnMinerva].prim[partMinerva].pdg)==13) muons++;
                if (abs(sr->mc.nu[ixnMinerva].prim[partMinerva].pdg)==211) pions++;}
         else{
        if (abs(sr->mc.nu[ixnMinerva].sec[partMinerva].pdg)==13) muons++;
                if (abs(sr->mc.nu[ixnMinerva].sec[partMinerva].pdg)==211) pions++;

         }
          
      }

      double dX = (end_pos.x - start_pos.x);
      double dY = (end_pos.y - start_pos.y);
      double dZ = (end_pos.z - start_pos.z);
      double length = TMath::Sqrt(dX * dX + dY * dY + dZ * dZ);
      double dirX = dX / length;
      double dirY = dY / length;
      double dirZ = dZ / length;
      
      double cosLReco =
          dirX * beam_x + dirY * beam_y + dirZ * beam_z;
      // If you got one matched to downstream Mx2, then it likely was the muon

      // Change the muon to make sure it is in the frame of the downward going beam

      recoCosL->Fill(cosLReco);


      goodEvents.push_back(n);
      double trueCosL = -999;
      // Now that we got a good reco interaction, just check the hadron info and muon info
      if (mcOnly) {
        if (maxEventTyp == 1) {
          minervaMatchPDG->Fill(sr->mc.nu[maxEventIxn].prim[maxEventPar].pdg);
          int minervaPDG = sr->mc.nu[maxEventIxn].prim[maxEventPar].pdg;
          if (abs(minervaPDG) == 13) {
            minervaMatchE->Fill(sr->mc.nu[maxEventIxn].prim[maxEventPar].p.E);
            const auto& start_posMnv =
                sr->mc.nu[maxEventIxn].prim[maxEventPar].start_pos;
            const auto& end_posMnv = sr->mc.nu[maxEventIxn].prim[maxEventPar].end_pos;

            double dX = (end_posMnv.x - start_posMnv.x);
            double dY = (end_posMnv.y - start_posMnv.y);
            double dZ = (end_posMnv.z - start_posMnv.z);

            const auto& p = sr->mc.nu[maxEventIxn].prim[maxEventPar].p;
            double length =
                TMath::Sqrt(p.px * p.px + p.py * p.py + p.pz * p.pz);
            cosL = p.pz / length;

            if (abs(startZMuonCand - start_posMnv.z) < 5)
              minervaMatchCos->Fill(cosL);
          }
        }

        // std::cout<<"counting hadrons (kaons are negligible)"<<std::endl;
        if (sr->mc.nu[biggestMatchIndex].prim.size() < 500 &&
            sr->mc.nu[biggestMatchIndex].targetPDG == 1000180400) {
          for (int primaries = 0;
               primaries < sr->mc.nu[biggestMatchIndex].prim.size();
               primaries++) {
            const auto& start_pos =
                sr->mc.nu[biggestMatchIndex].prim[primaries].start_pos;
            const auto& end_pos = sr->mc.nu[biggestMatchIndex].prim[primaries].end_pos;
            // if (std::isnan(start_pos.z)) continue;
            // std::cout<<"Got a position"<<std::endl;
            const auto& p = sr->mc.nu[biggestMatchIndex].prim[primaries].p;
            double dX = p.px; //(end_pos.x-start_pos.x);
            double dY = p.py; //(end_pos.y-start_pos.y);
            double dZ = p.pz; //(end_pos.z-start_pos.z);
            double length = TMath::Sqrt(dX * dX + dY * dY + dZ * dZ);
            double dirX = dX / length;
            double dirY = dY / length;
            double dirZ = dZ / length;
            double dXLen = (end_pos.x - start_pos.x);
            double dYLen = (end_pos.y - start_pos.y);
            double dZLen = (end_pos.z - start_pos.z);
            double lengthPos =
                TMath::Sqrt(dXLen * dXLen + dYLen * dYLen + dZLen * dZLen);

            if (!std::isnan(start_pos.z) &&
                (abs(sr->mc.nu[biggestMatchIndex].prim[primaries].pdg) ==
                     211 /*||
                            abs(sr->mc.nu[biggestMatchIndex].prim[primaries].pdg)==13*/
                 || abs(sr->mc.nu[biggestMatchIndex].prim[primaries].pdg) ==
                        2212)) {
              trueTrkLen->Fill(lengthPos);
              int cc = sr->mc.nu[biggestMatchIndex].iscc;
              int mode = sr->mc.nu[biggestMatchIndex].mode;
              if (cc == 1 && (mode == 1 || mode == 1001 || mode == 10))
                trueTrkLenQE->Fill(lengthPos);
              if (cc == 1 && (mode == 3 || mode == 4))
                trueTrkLenNonQE->Fill(lengthPos);
            }
            if (abs(sr->mc.nu[biggestMatchIndex].prim[primaries].pdg) == 2212) {
              trueProtonEWithRecoInt->Fill(
                  sr->mc.nu[biggestMatchIndex].prim[primaries].p.E - 0.938272);
              trueProtonWithRecoIntDirZ->Fill(dirZ);
              trueProtonWithRecoIntDirX->Fill(dirX);
              trueProtonWithRecoIntDirY->Fill(dirY);
            }
            if (abs(sr->mc.nu[biggestMatchIndex].prim[primaries].pdg) == 211) {
              truePionEWithRecoInt->Fill(
                  sr->mc.nu[biggestMatchIndex].prim[primaries].p.E - 0.13957);
              truePionWithRecoIntDirZ->Fill(dirZ);
              truePionWithRecoIntDirX->Fill(dirX);
              truePionWithRecoIntDirY->Fill(dirY);
            }

          }
        }

        for (long unsigned npart = 0;
             npart < sr->common.ixn.dlp[nixn].part.dlp.size();
             npart++) { // loop over particles
          // std::cout<<"Containment variable:
          // "<<sr->common.ixn.dlp[nixn].part.dlp[npart].contained<<std::endl;
          if (!sr->common.ixn.dlp[nixn].part.dlp[npart].primary)
            continue;
          int pdg = sr->common.ixn.dlp[nixn].part.dlp[npart].pdg;

          const auto& start_pos = sr->common.ixn.dlp[nixn].part.dlp[npart].start;
         const auto& end_pos = sr->common.ixn.dlp[nixn].part.dlp[npart].end;
          double dX = (end_pos.x - start_pos.x);
          double dY = (end_pos.y - start_pos.y);
          double dZ = (end_pos.z - start_pos.z);
          double length = TMath::Sqrt(dX * dX + dY * dY + dZ * dZ);
          double dirX = dX / length;
          double dirY = dY / length;
          double dirZ = dZ / length;
          auto lengthDist = length;
          if (abs(pdg) == 22 || abs(pdg) == 11 || abs(pdg) == 111) {
            auto truthSize =
                sr->common.ixn.dlp[nixn].part.dlp[npart].truth.size();
            double maxPartTruthOverlap = 0;

            for (int backTrack = 0; backTrack < truthSize; backTrack++) {
              int parType = sr->common.ixn.dlp[nixn]
                                .part.dlp[npart]
                                .truth[backTrack]
                                .type;
              int partNumber = sr->common.ixn.dlp[nixn]
                                   .part.dlp[npart]
                                   .truth[backTrack]
                                   .part;
              int interactionNumber =
                  sr->common.ixn.dlp[nixn].part.dlp[npart].truth[backTrack].ixn;

              double partTruthOverlap = sr->common.ixn.dlp[nixn]
                                            .part.dlp[npart]
                                            .truthOverlap[backTrack];
              if (maxPartTruthOverlap < partTruthOverlap) {
                maxPartTruthOverlap = partTruthOverlap;
              }
            }

            for (int backTrack = 0; backTrack < truthSize; backTrack++) {
              int parType = sr->common.ixn.dlp[nixn]
                                .part.dlp[npart]
                                .truth[backTrack]
                                .type;
              int partNumber = sr->common.ixn.dlp[nixn]
                                   .part.dlp[npart]
                                   .truth[backTrack]
                                   .part;
              int interactionNumber =
                  sr->common.ixn.dlp[nixn].part.dlp[npart].truth[backTrack].ixn;
              double partTruthOverlap = sr->common.ixn.dlp[nixn]
                                            .part.dlp[npart]
                                            .truthOverlap[backTrack];

              int backtracked = -9999;
              if (partTruthOverlap > maxPartTruthOverlap) {
                if (parType < 3) {
                  backtracked =
                      sr->mc.nu[interactionNumber].prim[partNumber].pdg;
                  primaryTrkIndex.push_back(partNumber);
                } else
                  backtracked =
                      sr->mc.nu[interactionNumber].sec[partNumber].pdg;
                correctShower = 0;
                if (parType > 2)
                  correctShower = 3;
                else if (backtracked == 11 || backtracked == 22 ||
                         backtracked == 111 || backtracked == 2112) {
                  correctShower = 1;
                }

              } else
                backtracked = -1;
            }
            if (sr->common.ixn.dlp[nixn].part.dlp[npart].primary == true)
              showerCorrectness->Fill(correctShower);
          }
          int maxPartNumber = -999;
          int maxTypeNumber = -999;
          if ((abs(pdg) == 2212 || abs(pdg) == 13 || abs(pdg) == 211 ||
               abs(pdg) == 321)) {
            auto truthSize =
                sr->common.ixn.dlp[nixn].part.dlp[npart].truth.size();
            double maxPartTruthOverlap = 0;
            double diffVertexdZ =
                abs(start_pos.z - sr->common.ixn.dlp[nixn].vtx.z);
            double diffVertexdX =
                abs(start_pos.x - sr->common.ixn.dlp[nixn].vtx.x);
            double diffVertexdY =
                abs(start_pos.y - sr->common.ixn.dlp[nixn].vtx.y);
            double diffVertex = TMath::Sqrt(diffVertexdZ * diffVertexdZ +
                                            diffVertexdY * diffVertexdY +
                                            diffVertexdX * diffVertexdX);
            diffVertexHist->Fill(diffVertex);

            for (int backTrack = 0; backTrack < truthSize; backTrack++) {
              int parType = sr->common.ixn.dlp[nixn]
                                .part.dlp[npart]
                                .truth[backTrack]
                                .type;
              int partNumber = sr->common.ixn.dlp[nixn]
                                   .part.dlp[npart]
                                   .truth[backTrack]
                                   .part;
              int interactionNumber =
                  sr->common.ixn.dlp[nixn].part.dlp[npart].truth[backTrack].ixn;
              double partTruthOverlap = sr->common.ixn.dlp[nixn]
                                            .part.dlp[npart]
                                            .truthOverlap[backTrack];
              // std::cout<<"Truth has overlap of:
              // "<<partTruthOverlap<<std::endl;
              if (maxPartTruthOverlap < partTruthOverlap) {
                maxPartTruthOverlap = partTruthOverlap;
              }
            }

            for (int backTrack = 0; backTrack < truthSize; backTrack++) {
              int parType = sr->common.ixn.dlp[nixn]
                                .part.dlp[npart]
                                .truth[backTrack]
                                .type;
              int partNumber = sr->common.ixn.dlp[nixn]
                                   .part.dlp[npart]
                                   .truth[backTrack]
                                   .part;
              int interactionNumber =
                  sr->common.ixn.dlp[nixn].part.dlp[npart].truth[backTrack].ixn;
              double partTruthOverlap = sr->common.ixn.dlp[nixn]
                                            .part.dlp[npart]
                                            .truthOverlap[backTrack];
              if (maxPartTruthOverlap > partTruthOverlap)
                continue;
              int backtracked = -9999;
              if (partTruthOverlap > 0.5) {
                maxPartNumber = partNumber;
                maxTypeNumber = parType;
                if (parType < 3)
                  backtracked =
                      sr->mc.nu[interactionNumber].prim[partNumber].pdg;
                else
                  backtracked =
                      sr->mc.nu[interactionNumber].sec[partNumber].pdg;
                if (parType < 3) {
                  int pdgNumber = 0;
                  if (abs(pdg) == 13) {
                    pdgNumber = 0;
                  } else if (abs(pdg) == 2212) {
                    pdgNumber = 1;
                  } else if (abs(pdg) == 211) {
                    pdgNumber = 2;
                  } else {
                    pdgNumber = 3;
                  }
                  int backtrackedPDG = 0;
                  if (abs(backtracked) == 13)
                    backtrackedPDG = 0;
                  else if (abs(backtracked) == 2212) {
                    backtrackedPDG = 1;
                  } else if (abs(backtracked) == 211) {
                    backtrackedPDG = 2;
                  } else {
                    backtrackedPDG = 3;
                  }

                  const auto& start_pos =
                      sr->mc.nu[interactionNumber].prim[partNumber].start_pos;
                  // if (std::isnan(start_pos.z)) continue;
                  const auto& p = sr->mc.nu[interactionNumber].prim[partNumber].p;
                  const auto& end_pos =
                      sr->mc.nu[interactionNumber].prim[partNumber].end_pos;

                  double dX = p.px; //(end_pos.x-start_pos.x);
                  double dY = p.py; //(end_pos.y-start_pos.y);
                  double dZ = p.pz; //(end_pos.z-start_pos.z);
                  double length = TMath::Sqrt(dX * dX + dY * dY + dZ * dZ);

                  double dXLen = (end_pos.x - start_pos.x);
                  double dYLen = (end_pos.y - start_pos.y);
                  double dZLen = (end_pos.z - start_pos.z);
                  double lengthPos = TMath::Sqrt(dXLen * dXLen + dYLen * dYLen +
                                                 dZLen * dZLen);
                  if (sr->common.ixn.dlp[nixn].part.dlp[npart].contained ==
                          true &&
                      lengthPos > minTrkLength)
                    confusionMatrix->Fill(pdgNumber, backtrackedPDG);

                  double dirX = dX / length;
                  double dirY = dY / length;
                  double dirZ = dZ / length;
                  if (!std::isnan(start_pos.z) &&
                      (backtracked == 2212 /*|| abs(backtracked)==13*/ ||
                       abs(backtracked) == 211) &&
                      parType < 3) {
                    selTrkLen->Fill(lengthPos);

                    int cc = sr->mc.nu[biggestMatchIndex].iscc;
                    int mode = sr->mc.nu[biggestMatchIndex].mode;
                    if (cc == 1 && (mode == 1 || mode == 1001 || mode == 10))
                      selTrkLenQE->Fill(lengthPos);
                    if (cc == 1 && (mode == 3 || mode == 4))
                      selTrkLenNonQE->Fill(lengthPos);
                  }
                  if (sr->mc.nu[biggestMatchIndex].targetPDG == 1000180400 &&
                      backtracked == 2212 && parType < 3) {
                    selProtonE->Fill(
                        sr->mc.nu[interactionNumber].prim[partNumber].p.E -
                        0.938272);
                    selProtonDirX->Fill(dirX);
                    selProtonDirZ->Fill(dirZ);
                    selProtonDirY->Fill(dirY);
                  }
                  if (sr->mc.nu[biggestMatchIndex].targetPDG == 1000180400 &&
                      abs(backtracked) == 211 && parType < 3) {
                    selPionE->Fill(
                        sr->mc.nu[interactionNumber].prim[partNumber].p.E -
                        0.13957039);
                    selPionDirX->Fill(dirX);
                    selPionDirZ->Fill(dirZ);
                    selPionDirY->Fill(dirY);
                  }
                }
                correctTrack = 0;
                // std::cout<<"Backtrackilng for "<<parType<<" :
                // "<<pdg<<","<<backtracked<<std::endl;
                if (parType > 2)
                  correctTrack = 3;
                else if (backtracked == 11 || backtracked == 22 ||
                         backtracked == 111 || backtracked == 2112) {
                  correctTrack = 1;
                }
              }
              if (partTruthOverlap < 0.5)
                correctTrack = 2;
              else
                backtracked = -1;

              if (lengthDist > minTrkLength)
                trackCorrectness->Fill(correctTrack);
            }
          }
        }
        int trueMatchMult = 0;

        if (goodInteraction) {
          recoBacktrackCCAr->Fill(sr->mc.nu[biggestMatchIndex].iscc);
          recoBacktrackPDGAr->Fill(sr->mc.nu[biggestMatchIndex].pdg);
          if (sr->mc.nu[biggestMatchIndex].iscc &&
              abs(sr->mc.nu[biggestMatchIndex].pdg) == 14) {

            for (int primaries = 0;
                 primaries < sr->mc.nu[biggestMatchIndex].prim.size();
                 primaries++) {
              double Elep = sr->mc.nu[biggestMatchIndex].prim[primaries].p.E;
              const auto& start_pos =
                  sr->mc.nu[biggestMatchIndex].prim[primaries].start_pos;
              const auto& end_pos =
                  sr->mc.nu[biggestMatchIndex].prim[primaries].end_pos;
              const auto& p = sr->mc.nu[biggestMatchIndex].prim[primaries].p;
              // if (std::isnan(start_pos.z)) continue;
              double dX = (p.px);
              double dY = (p.py);
              double dZ = (p.pz);
              double length =
                  TMath::Sqrt(p.px * p.px + p.py * p.py + p.pz * p.pz);
              double dir_x = dX / length;
              double dir_y = dY / length;
              double dir_z = dZ / length;
              cosL = dir_x * beam_x + dir_y * beam_y + dir_z * beam_z;
              dX = end_pos.x - start_pos.x;
              dY = end_pos.y - start_pos.y;
              dZ = end_pos.z - start_pos.z;
              double lengthTrk = TMath::Sqrt(dX * dX + dY * dY + dZ * dZ);

              int truePDG = sr->mc.nu[biggestMatchIndex].prim[primaries].pdg;
              if ((abs(truePDG) == 13 || abs(truePDG) == 211 ||
                   abs(truePDG) == 2212) &&
                  lengthTrk > minTrkLength)
                trueMatchMult++;
              // if (goodInteraction)
              // std::cout<<biggestMatch<<","<<sr->common.ixn.dlp[nixn].vtx.z<<","<<sr->mc.nu[biggestMatchIndex].vtx.z<<","<<start_pos.z<<","<<end_pos.z<<std::endl;
              if (abs(sr->mc.nu[biggestMatchIndex].prim[primaries].pdg) != 13)
                continue;
              recoBacktrackElAr->Fill(Elep);
              trueCosL = cosL;
              //  std::cout<<cosL<<std::endl;
            }
            // Fill response matrix for the lepton angle
            if (goodInteraction) {
              trueBacktrackedCosL_zoomOut->Fill(trueCosL);
              trueBacktrackedCosL_unbinned->Fill(trueCosL);
              if (trueCosL<0.91)       recoVertex2D->Fill(sr->common.ixn.dlp[nixn].vtx.x,
                         sr->common.ixn.dlp[nixn].vtx.z);
              if (trueCosL>0.91){trueBacktrackedCosL->Fill(trueCosL);
              responseCosL->Fill(cosLReco, trueCosL);
            }}
            // Fill the response matrix for the track multiplicity
            if (cosLReco > 0.9 && trueCosL > 0.9 &&
                goodInteraction) { /*std::cout<<trackMult<<","<<trueMatchMult<<std::endl;*/
              responseMult->Fill(trackMult, trueMatchMult);
              backtrackTrueMult->Fill(trueMatchMult);
              trueBacktrackedMult_unbinned->Fill(trueMatchMult);
            }
          }
        }
        // Fill the stack histograms of the track multiplicity by truth info
        if (cosLReco > 0.9 && (!goodInteraction || trueCosL < 0.9)) {
          track_multBad->Fill(trackMult);
        }
        if (!goodInteraction && sr->mc.nu[biggestMatchIndex].id > 1E9) {
          recoBacktrackCCRock->Fill(sr->mc.nu[biggestMatchIndex].iscc);
          recoBacktrackPDGRock->Fill(sr->mc.nu[biggestMatchIndex].pdg);
        }
        if (!goodInteraction && sr->mc.nu[biggestMatchIndex].id < 1E9) {
          recoBacktrackCCSec->Fill(sr->mc.nu[biggestMatchIndex].iscc);
          recoBacktrackPDGSec->Fill(sr->mc.nu[biggestMatchIndex].pdg);
        }

        if (!goodInteraction || trueCosL < 0.91)
          badRecoCosL->Fill(cosLReco);
        if (cosLReco > 0.9) {
          if (goodInteraction && trueCosL > 0.9 &&
              sr->mc.nu[biggestMatchIndex].iscc) {
            track_multGood->Fill(trackMult);
            int cc = sr->mc.nu[biggestMatchIndex].iscc;
            int mode = sr->mc.nu[biggestMatchIndex].mode;
            track_multMode->Fill(mode);
            if (cc == 1) {
              if (mode == 1 || mode == 1001)
                track_multQE->Fill(trackMult);
              if (mode == 10)
                track_multMEC->Fill(trackMult);
              if (mode == 3)
                track_multDIS->Fill(trackMult);
              if (mode == 4)
                track_multRES->Fill(trackMult);
              if (mode == 5)
                track_multCOH->Fill(trackMult);

            } else
              track_multNC->Fill(trackMult);
          } else {
            if (rock == 1)
              track_multRock->Fill(trackMult);
            else{
              track_multSec->Fill(trackMult);
            }
            bad_origin->Fill(rock);
            recoVertex2DBad->Fill(sr->common.ixn.dlp[nixn].vtx.x,
                                  sr->common.ixn.dlp[nixn].vtx.z);
            trueVertex2DBad->Fill(trueVtxX, trueVtxZ);
            recoVertex2DBadYZ->Fill(sr->common.ixn.dlp[nixn].vtx.y,
                                    sr->common.ixn.dlp[nixn].vtx.z);
          }
        }
      }
      // Fill the reco information for track multiplicity
      if (cosLReco > 0.9) {
        track_mult->Fill(trackMult);
        part_mult->Fill(partMult);
      }

    } // end for interaction
  }}   // end for spills
  std::cout << "Total: " << recoCosL->GetEntries() << std::endl;
  std::cout << "All Interactions No Additional Selection: "
            << double(goodIntNum) / trueCosL->GetEntries() << ","
            << float(goodIntNum) / float(goodIntNum + badIntNum) << ","<<  float(goodIntNum + badIntNum) << std::endl;
  std::cout << "In fiducial volume: "
            << double(goodIntInFidVol) / trueCosL->GetEntries() << ","
            << float(goodIntInFidVol) / float(goodIntInFidVol + badIntInFidVol) << "," << float(goodIntInFidVol + badIntInFidVol)
            << std::endl;
  std::cout << "All events with MINERvA match: " << (recoCosL->GetEntries()-badRecoCosL->GetEntries())
            << "/" << trueCosL->GetEntries() << ","
            << (recoCosL->GetEntries()-badRecoCosL->GetEntries()) / trueCosL->GetEntries() << ","
            << (recoCosL->GetEntries()-badRecoCosL->GetEntries())/recoCosL->GetEntries()<< "," <<recoCosL->GetEntries()<< std::endl;
  // Create output file and write your histograms
  TFile *caf_out_file = new TFile(output_rootfile.c_str(), "recreate");
  matchHistCosl->Write();
  matchHistEl->Write();
  histEl->Write();
  trueIntPerSpill->Write();
  part_mult->Write();
  part_multTrkOnly->Write();
  part_energy_hist->Write();
  track_mult->Write();
  track_multGood->Write();
  track_multBad->Write();
  track_length->Write();
  track_multNC->Write();
  track_multMode->Write();
  track_multQE->Write();
  track_multDIS->Write();
  track_multMEC->Write();
  track_multRES->Write();
  track_multCOH->Write();
  track_multRock->Write();
  track_multSec->Write();
  diffVertexHist->Write();
  minervaMatchPDG->Write();
  minervaMatchCos->Write();
  minervaMatchE->Write();
  responseCosL->Write();
  recoBacktrackCCAr->Write();
  recoBacktrackCCRock->Write();
  recoBacktrackCCSec->Write();
  recoBacktrackPDGSec->Write();
  recoBacktrackPDGAr->Write();
  recoBacktrackPDGRock->Write();

  recoBacktrackElAr->Write();
  trueBacktrackedCosL->Write();
  trueBacktrackedCosL_zoomOut->Write();
  trueCosL_unbinned->Write();
  trueCosL_zoomOut->Write();
  trueCosLNumubar_zoomOut->Write();
  trueCosLNumu_zoomOut->Write();
  trueBacktrackedCosL_unbinned->Write();
  true_mult->Write();
  true_multTrkOnly->Write();
  true_multGENIE->Write();
  recoHistVertexX->Write();
  recoHistVertexZ->Write();
  recoHistVertexY->Write();

  matchTrue_mult->Write();
  matchTrue_multTrkOnly->Write();
  matchTrue_multGENIE->Write();
  bad_origin->Write();
  recoVertex2D->Write();

  showerCorrectness->Write();
  trackCorrectness->Write();
  selPionE->Write();
  selProtonE->Write();

  truePionEWithRecoInt->Write();
  trueProtonEWithRecoInt->Write();
  trueEnu->Write();



  auto protonKEEff = new TEfficiency(*selProtonE, *trueProtonEWithRecoInt);
  auto pionKEEff = new TEfficiency(*selPionE, *truePionEWithRecoInt);

  auto protonDirX = new TEfficiency(*selProtonDirX, *trueProtonWithRecoIntDirX);
  auto protonDirY = new TEfficiency(*selProtonDirY, *trueProtonWithRecoIntDirY);
  auto protonDirZ = new TEfficiency(*selProtonDirZ, *trueProtonWithRecoIntDirZ);

  auto pionDirX = new TEfficiency(*selPionDirX, *truePionWithRecoIntDirX);
  auto pionDirY = new TEfficiency(*selPionDirY, *truePionWithRecoIntDirY);
  auto pionDirZ = new TEfficiency(*selPionDirZ, *truePionWithRecoIntDirZ);


  auto effTrkLen = new TEfficiency(*selTrkLen, *trueTrkLen);
  auto effTrkLenQE = new TEfficiency(*selTrkLenQE, *trueTrkLenQE);
  auto effTrkLenNonQE = new TEfficiency(*selTrkLenNonQE, *trueTrkLenNonQE);

  effTrkLen->Write();
  effTrkLenQE->Write();
  effTrkLenNonQE->Write();
  protonDirX->Write();
  protonDirY->Write();
  protonDirZ->Write();
  pionDirX->Write();
  pionDirY->Write();
  pionDirZ->Write();

  protonKEEff->Write();
  pionKEEff->Write();

  responseMult->Write();
  trueMult_unbinned->Write();
  trueBacktrackedMult_unbinned->Write();
  backtrackTrueMult->Write();
  responseGenieToG4->Write();

  confusionMatrix->Write();
  recoCosL->Write();
  badRecoCosL->Write();
  trueCosL->Write();
  trueDiffPosvsPDirZ->Write();
  totalSpills->Write();
  totalPOT->Fill(0.5, sumPOT);
  std::cout<<totalPOT->GetBinContent(1)<<std::endl;
  std::cout<<totalSpills->GetBinContent(1)<<std::endl;
  totalPOT->Write();
  caf_out_file->Close();
  std::cout<<"Total number of neutrinos in Ar (and Fid. Vol.): "<<totalNuAr<<","<<totalNuArFidVol<<std::endl;
  std::cout<<"Muon vs. Pions vs. other: "<<muons<<","<<pions<<","<<recoCosL->GetEntries()-pions-muons<<std::endl;
  return 1;
}

int main(int argc, char **argv) {

  if (argc != 4) {
    std::cout << "\n USAGE: " << argv[0]
              << "input_caf_file_list output_root_file\n"
              << std::endl;
    return 1;
  }

  std::string input_file_list = argv[1];
  std::string output_rootfile = argv[2];
  std::string mcOnlyString = argv[3];
  bool mcOnly = true;
  if (mcOnlyString == "0")
    mcOnly = false;
  std::cout << mcOnly << "," << argv[3] << std::endl;
  caf_plotter(input_file_list, output_rootfile, mcOnly);

  return 0;
}
