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
#include "TProfile.h"
#include "../fluxSyst/FluxCovarianceReweight.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "duneanaobj/StandardRecord/StandardRecord.h" //Ideally, this should be SRProxy.h, but there is an include error for that now. Alternatively, you can use SetBranchStatus function in TreeLoader, but it does not work for the common branch (to do)
#include <fstream>
#include <iostream>
#include <string>

double minTrkLength = 3;

  bool Passes_cut(caf::SRTrack track_minerva, double x1_lar, double x2_lar, double y1_lar, double y2_lar, double z1_lar, double z2_lar, double &costheta, double &residual)
  {

    double z_extr=-70;
    double d_x = 17;
    double d_y = 19;
    double d_thetax =0.08;
    double d_thetay=0.09;
      
    double x1_minerva = track_minerva.start.x;
    double x2_minerva = track_minerva.end.x;
    double y1_minerva = track_minerva.start.y;
    double y2_minerva = track_minerva.end.y;
    double z1_minerva = track_minerva.start.z;
    double z2_minerva = track_minerva.end.z;

    double dX=x2_lar-x1_lar;
    double dY=y2_lar-y1_lar;
    double dZ=z2_lar-z1_lar;
    double track_LarLen=TMath::Sqrt(dX*dX+dY*dY+dZ*dZ);

    /*
    The experimental setup: Liquid Argon Detector is placed betbeen two MINERvA planes.
    To define matching criteria it is needed to find angles between LAr and MINERvA tracks.
    For LAr detector resolution is diffeent in X and Y direction, therefore it is needed to find angles between tracks
    as finction of the angle in X direction and as the function of an angle in Y direction. Distances between tracks
    will be calculated as distancec between extrapolated points - points of intersection of LAr and MINERvA tracls with
    the plane (parallel to plane XY) of LAr detector.
    */

    double tg_theta_mn_x = (x2_minerva - x1_minerva) / (z2_minerva - z1_minerva); // tangent of an angle between minerva track and X-axis
    double tg_theta_mn_y = (y2_minerva - y1_minerva) / (z2_minerva - z1_minerva); // tangent of an angle between minerva track and Y-axis
    double theta_mn_x = atan(tg_theta_mn_x);                                      // angle between minerva track and X-axis
    double theta_mn_y = atan(tg_theta_mn_y);                                      // angle between minerva track and Y-axis

    double tg_theta_nd_x = (x2_lar - x1_lar) / (z2_lar - z1_lar); // tangent of the angle between LAr track and X-axis
    double tg_theta_nd_y = (y2_lar - y1_lar) / (z2_lar - z1_lar); // tangent of the angle between LAr track and Y-axis
    double theta_nd_x = atan(tg_theta_nd_x);                      // angle between LAr track and X-axis
    double theta_nd_y = atan(tg_theta_nd_y);                      // angle between LAr track and Y-axis

    double delta_theta_x = theta_mn_y - theta_nd_y;
    double delta_theta_y = theta_mn_x - theta_nd_x;

    // Extrapolating Both tracks to the same point z = zextr (here it's the front of Lar)
    double t_mn = (z_extr - z1_minerva) / (z2_minerva - z1_minerva);
    double x_mn = t_mn * (x2_minerva - x1_minerva) + x1_minerva; // X-coordinate of extrapolated point of LAr track
    double y_mn = t_mn * (y2_minerva - y1_minerva) + y1_minerva; // Y-coordinate of extrapolated point of LAr track

    double t_nd = (z_extr - z1_lar) / (z2_lar - z1_lar); // parametr of the equation of the line (LAr track)
    double x_nd = t_nd * (x2_lar - x1_lar) + x1_lar;     // X-coordinate of extrapolated point of LAr track
    double y_nd = t_nd * (y2_lar - y1_lar) + y1_lar;     // Y-coordinate of extrapolated point of LAr track

    double dist_x = (x_mn - x_nd); // distance between X-coordinates of extrapolated points of minerva and LAr tracks
    double dist_y = (y_mn - y_nd); // distance between Y-coordinates of extrapolated points of minerva and LAr tracks

    residual = sqrt(pow(dist_x, 2) + pow(dist_y, 2));
    costheta = ((x2_minerva - x1_minerva) * (x2_lar - x1_lar) +
                (y2_minerva - y1_minerva) * (y2_lar - y1_lar) +
                (z2_minerva - z1_minerva) * (z2_lar - z1_lar)) /
               (track_LarLen * track_minerva.len_cm); // angle between minerva and Lar tracks

    return (abs(delta_theta_x) < d_thetax && abs(delta_theta_y) < d_thetay && abs(dist_x) < d_y && abs(dist_y) < d_x);
  }

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
  TFile fFlux("../fluxSyst/numiFluxSyst.root");
  TH1D* hflux_rhc_nue=(TH1D*)fFlux.Get("flux_prediction/hflux_rhc_nue");
  TH1D* hflux_rhc_nuebar=(TH1D*)fFlux.Get("flux_prediction/hflux_rhc_nuebar");
      TH1D* hflux_rhc_numu=(TH1D*)fFlux.Get("flux_prediction/hflux_rhc_numu");
      TH1D* hflux_rhc_numubar=(TH1D*)fFlux.Get("flux_prediction/hflux_rhc_numubar");
    FluxCovarianceReweight rw;
    
    rw.LoadBinning("../fluxSyst/flux_covariance_binning_NuMI_GeV.txt");
    rw.LoadCovariance();
rw.LoadFluxHistograms(
  hflux_rhc_nue,
  hflux_rhc_nuebar,
  hflux_rhc_numu,
  hflux_rhc_numubar
);
rw.GenerateThrows(1000);
    double Phi_nom = rw.GetIntegratedFluxNominal();


  // Set the beam direction  
  const auto beam_dir = TVector3(0.0, -0.05836, 1.0); //(-3.34 degrees in Y)
  double beam_x = 0;
  double beam_y = -0.05836 / TMath::Sqrt(0.05836 * 0.05836 + 1);
  double beam_z = 1.0 / TMath::Sqrt(0.05836 * 0.05836 + 1);
  // Define histograms
  int binsMult = 7;
  Double_t edgesMult[8] = {1, 2, 3, 4, 5, 6, 8, 10};
  TH1D *track_mult = new TH1D("track_mult", "track_mult", binsMult, edgesMult);
  TH1D *true_mult = new TH1D("true_mult", "true_mult", binsMult, edgesMult);
  TH1D *true_multTrkOnly = new TH1D("true_multTrkOnly", "true_multTrkOnly", binsMult, edgesMult);

  TH1D *track_multBad =
      new TH1D("track_multBad", "track_multBad", binsMult, edgesMult);
  TH1D *track_multGood =
      new TH1D("track_multGood", "track_multGood", binsMult, edgesMult);

  TH2D *responseMult = new TH2D("responseMult", "responseMult", binsMult,
                                edgesMult, binsMult, edgesMult);
  TH1D *backtrackTrueMult =
      new TH1D("backtrackTrueMult", "backtrackTrueMult", binsMult, edgesMult);
  Double_t edges[7] = {0.91, 0.96, 0.98, 0.9887, 0.994, 0.9974, 1};
  

  TH1D *recoCosL = new TH1D("recoCosL", "recoCosL", 6, edges);
  TH1D *badRecoCosL = new TH1D("badRecoCosL", "badRecoCosL", 6, edges);
  TH1D *trueBacktrackedCosL =
      new TH1D("backtrackTrueCosL", "backtrackTrueCosL", 6, edges);
  TH1D *trueCosL = new TH1D("trueCosL", "trueCosL", 6, edges);
  TH2D *responseCosL =
      new TH2D("responseCosL", "responseCosL", 6, edges, 6, edges);

  TH1D *totalPOT = new TH1D("totalPOT", "totalPOT", 1, 0, 1);
  TH1D *totalSpills = new TH1D("totalSpilles", "totalSpills", 1, 0, 1);
  std::vector<TH1D> trueCosLVec, true_multVec, recoCosLVec, recoMultVec, backtrackCosLVec, backtrackMultVec, badRecoCosLVec, badRecoMultVec;
  std::vector<TH2D> responseCosLVec, responseMultVec;

  for (int index=0; index<1000; index++){
   TH1D* tempHist=(TH1D*)trueCosL->Clone(Form("trueCosL_%d",index));
   TH1D* multTemp=(TH1D*)true_mult->Clone(Form("true_mult_%d",index));

   TH2D* tempResp=(TH2D*)responseCosL->Clone(Form("responseCosL_%d",index));
   TH2D* tempRespMult=(TH2D*)responseMult->Clone(Form("responseMult_%d",index));
   TH1D* tempReco=(TH1D*)recoCosL->Clone(Form("recoCosL_%d",index));
   TH1D* tempMult=(TH1D*)track_mult->Clone(Form("recoMult_%d",index));

   TH1D* tempBadReco=(TH1D*)recoCosL->Clone(Form("badRecoCosL_%d",index));
   TH1D* tempBadMult=(TH1D*)track_mult->Clone(Form("badRecoMult_%d",index));

   TH1D* tempMatchCosL=(TH1D*)trueBacktrackedCosL->Clone(Form("backtrackTrueCosL_%d",index));
   TH1D* tempMatchMult=(TH1D*)backtrackTrueMult->Clone(Form("backtrackTrueMult_%d",index));
  trueCosLVec.push_back(*tempHist);
    true_multVec.push_back(*multTemp);
    responseCosLVec.push_back(*tempResp);
    responseMultVec.push_back(*tempRespMult);
  recoCosLVec.push_back(*tempReco);
  recoMultVec.push_back(*tempMult);
  backtrackCosLVec.push_back(*tempMatchCosL);
     backtrackMultVec.push_back(*tempMatchMult);
    badRecoCosLVec.push_back(*tempBadReco);
    badRecoMultVec.push_back(*tempBadMult);
  }
  recoCosL->SetTitle(Form("%0.08f",Phi_nom)); 
  double sumPOT = 0;

  // Give an input list
  std::ifstream caf_list(input_file_list.c_str());

  // Check if input list is present
  if (!caf_list.is_open()) {
    std::cerr << Form("File %s not found", input_file_list.c_str())
              << std::endl;
    return 1;
  }

  // Add files to CAF chain from input list
  std::string tmp;
  TChain *caf_chain = new TChain("cafTree");
  int index=0;
  while(caf_list >> tmp){

    index++;

    TFile * f = new TFile(tmp.c_str());
    TTree * caf_chain = (TTree*) f->Get("cafTree");  

  //Check if CAF tree is present
  long Nentries = caf_chain->GetEntries();
  // Check if CAF tree is present
  if (!caf_chain) {
    std::cerr << Form("There is no tree in %s", tmp.c_str()) << std::endl;
    return 1;
  }


  // Define Standard Record and link it to the CAF tree branch "rec"
  auto sr = new caf::StandardRecord;
  caf_chain->SetBranchAddress("rec", &sr);
  totalSpills->Fill(0.5, Nentries);
  bool skipEvent = false;
  // Run over the events


  for (long n = 0; n < Nentries; ++n) {

    double trueVtxX = -9999;
    double trueVtxY = -999;
    double trueVtxZ = -9999;
    skipEvent = false;
   // if (n % 10000 == 0)
     // std::cout << Form("Processing trigger %ld of %ld", n, Nentries)
       //         << std::endl;
    caf_chain->GetEntry(n); // Get spill from tree
    sumPOT = sr->beam.pulsepot / 1e13 + sumPOT;
    int intAboveThresh=0;

    bool hasANeutrino = false;
    double mnvOffsetX = -10;
    double mnvOffsetY = 5;

      // Save truth distributions
      for (long unsigned ntrue = 0; ntrue < sr->mc.nu.size(); ntrue++) {
        auto vertex = sr->mc.nu[ntrue].vtx;
        trueVtxX = vertex.x;
        trueVtxY = vertex.y;
        trueVtxZ = vertex.z;

        auto truePrimary = sr->mc.nu[ntrue].prim;
        auto trueSecondary = sr->mc.nu[ntrue].sec;
        int aboveThresh=0;


        int truePart = 0;
        int truePartTrkOnly = 0;
        int truePartNoG4 = 0;
        // We only want CC numu interactions on argon in the fiducial volume
        if (abs(sr->mc.nu[ntrue].pdg) != 14)
          continue;
        if (abs(sr->mc.nu[ntrue].targetPDG) != 1000180400)
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

            auto start_pos = sr->mc.nu[ntrue].prim[primaries].start_pos;
            auto end_pos = sr->mc.nu[ntrue].prim[primaries].end_pos;
            // if (std::isnan(start_pos.z)){  continue;}
            auto p = sr->mc.nu[ntrue].prim[primaries].p;
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
            }

              
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
        int gIdx=sr->mc.nu[ntrue].genieIdx;

        int trueMatchMult = 0;
        int covarIndex=-1;
        Double_t totalWeight[1000];

        covarIndex=rw.FindBinIndex(sr->mc.nu[ntrue].pdg,sr->mc.nu[ntrue].E);
        if (covarIndex>-1 && covarIndex<151)
        rw.FillAllThrowsForBin(covarIndex, totalWeight);

        if (cosL > 0.91 && Elep>1) {
            trueCosL->Fill(cosL);
            for (int index=0;  index<1000; index++){
            if (covarIndex>-1 && covarIndex<151){
           if (totalWeight[index]>0.01 && totalWeight[index]<20){  trueCosLVec.at(index).Fill(cosL,totalWeight[index]);
            true_multVec.at(index).Fill(truePartTrkOnly,totalWeight[index]);}
            else{
            trueCosLVec.at(index).Fill(cosL,1);
            true_multVec.at(index).Fill(truePartTrkOnly,1);



            }
                                                                
                                                                
                                                                }

            }
          }
        true_mult->Fill(truePartTrkOnly);
        true_multTrkOnly->Fill(truePartTrkOnly);

        hasANeutrino = true;
    }
    
    for (long unsigned nixn = 0; nixn < sr->common.ixn.dlp.size(); nixn++) {
      double cosL = -999;
      int rock=0;
      bool    goodInteraction = false;
      bool goodOutOfSigDef=false;
          double biggestMatch = -999;
      int biggestMatchIndex = -999;
          double maxDotProductDS = -999;
      std::vector<int> mx2IntCandidateVector, mx2IdxCandidateVector, ndlarTrkCandidateVector;
      double recoVertexX = sr->common.ixn.dlp[nixn].vtx.x;
      double recoVertexY = sr->common.ixn.dlp[nixn].vtx.y;
      double recoVertexZ = sr->common.ixn.dlp[nixn].vtx.z;


        for (long unsigned int ntruth = 0; ntruth < sr->common.ixn.dlp[nixn].truth.size();
             ntruth++) {

          if (biggestMatch < sr->common.ixn.dlp[nixn].truthOverlap.at(ntruth) && sr->common.ixn.dlp[nixn].truthOverlap.at(ntruth)>0.5) {
            biggestMatch = sr->common.ixn.dlp[nixn].truthOverlap.at(ntruth);
            biggestMatchIndex = sr->common.ixn.dlp[nixn].truth.at(ntruth);
            
          }
        }
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
            diffVtx < 20 && abs(trueVtxX) < 59 && abs(trueVtxZ) > 5 && abs(trueVtxX)>5 &&
            abs(trueVtxZ) < 59.5 && abs(trueVtxY) < 57){
            goodOutOfSigDef=true;
                        for (int primaries = 0;
                 primaries < sr->mc.nu[biggestMatchIndex].prim.size();
                 primaries++) {
            if (abs(sr->mc.nu[biggestMatchIndex].prim[primaries].pdg) != 13)
                continue;
              double Elep = sr->mc.nu[biggestMatchIndex].prim[primaries].p.E;
              auto start_pos =
                  sr->mc.nu[biggestMatchIndex].prim[primaries].start_pos;
              auto end_pos =
                  sr->mc.nu[biggestMatchIndex].prim[primaries].end_pos;
              auto p = sr->mc.nu[biggestMatchIndex].prim[primaries].p;
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
        }
        if (/*abs(abs(sr->common.ixn.dlp[nixn].vtx.x)-33)<1 || */ abs(
              sr->common.ixn.dlp[nixn].vtx.x) > 59 ||
          abs(sr->common.ixn.dlp[nixn].vtx.x) < 5 ||
          abs(sr->common.ixn.dlp[nixn].vtx.y) > 57 ||
          abs(sr->common.ixn.dlp[nixn].vtx.z) < 5 ||
          abs(sr->common.ixn.dlp[nixn].vtx.z) > 59.5)
        continue;


      int partMult = 0;
      int partMultTrkOnly = 0;
      int trackMult= 0;
      for (long unsigned npart = 0;
           npart < sr->common.ixn.dlp[nixn].part.dlp.size();
           npart++) { // loop over particles
      //  if (!sr->common.ixn.dlp[nixn].part.dlp[npart].primary)
      //    continue;
        int pdg = sr->common.ixn.dlp[nixn].part.dlp[npart].pdg;
        // Loop over primary tracks
        if ((abs(pdg) == 2212 || abs(pdg) == 13 || abs(pdg) == 211 ||
             abs(pdg) == 321)) {

          auto start_pos = sr->common.ixn.dlp[nixn].part.dlp[npart].start;
          auto end_pos = sr->common.ixn.dlp[nixn].part.dlp[npart].end;
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

          if (dirZ < 0) {
            dirZ = -dirZ;
            dirX = -dirX;
            dirY = -dirY;
            auto temp = start_pos;
            end_pos = start_pos;
            end_pos = temp;
          }
          if (std::isnan(start_pos.z))
            length = -999;

          // Make sure it is above the track threshold
          if (sr->common.ixn.dlp[nixn].part.dlp[npart].primary == true &&
              length > minTrkLength) {
            partMult++;
            trackMult++;
            // see if it punches out and match it to MINERvA
            int maxPartMinerva = -999;
            int maxTypeMinerva = -999;
            int maxIxnMinerva = -999;
            int maxPartMinervaUS = -999;
            int maxTypeMinervaUS = -999;
                  int mx2IntCandidate=-9;
      int mx2IdxCandidate=-9;
            if ((abs(start_pos.z) > 62 || abs(end_pos.z) > 62)){
              int minervaPass = 0;
              double dotProductDS = -999;
              double deltaExtrapYUS = -999;
              double deltaExtrapY = -999;
              double dotProductUS = -999;
              double deltaExtrapX = -999;
              double deltaExtrapXUS = -999;
              double dotProductFromCAF=-999;
        	 	for(int i=0; i<sr->nd.trkmatch.extrap.size(); i++){
                if (sr->nd.trkmatch.extrap[i].larid.ixn==nixn && sr->nd.trkmatch.extrap[i].larid.reco==1){
                    int index=sr->nd.trkmatch.extrap[i].larid.idx;
                    if (sr->nd.lar.dlp[nixn].tracks[index].start.z==sr->common.ixn.dlp[nixn].part.dlp[npart].start.z && sr->nd.lar.dlp[nixn].tracks[index].end.z==sr->common.ixn.dlp[nixn].part.dlp[npart].end.z){
                    double angldispl=abs(sr->nd.trkmatch.extrap[i].angdispl);
                    if (dotProductFromCAF<angldispl) dotProductFromCAF=angldispl;
                }}
                }
                            if (dotProductFromCAF<0) continue;
    
        
                    for(int i=0; i<sr->nd.minerva.ixn.size(); i++){
                
                    for (int j=0; j<sr->nd.minerva.ixn[i].ntracks; j++){
                    double dotProductTemp=0;
                    double trackDispl=0;
                    bool pass=Passes_cut(sr->nd.minerva.ixn[i].tracks[j], start_pos.x,end_pos.x,start_pos.y,end_pos.y,start_pos.z, end_pos.z, dotProductTemp, trackDispl);
                    if (!pass) continue;
                    dotProductTemp=abs(dotProductTemp);
                    if (mcOnly && abs(dotProductTemp-dotProductFromCAF)>0.001) continue;
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
                      mx2IntCandidate=i;
                      mx2IdxCandidate=j;

                      if (mcOnly) {
                        maxPartMinerva =
                            sr->nd.minerva.ixn[i].tracks[j].truth[0].part;
                        maxTypeMinerva =
                            sr->nd.minerva.ixn[i].tracks[j].truth[0].type;
                        maxIxnMinerva =
                            sr->nd.minerva.ixn[i].tracks[j].truth[0].ixn;
                      }

                    }
                  }
                        
          }   
            if (dotProductDS > maxDotProductDS) {
                maxDotProductDS = dotProductDS;
                mx2IntCandidateVector.push_back(mx2IntCandidate);
                mx2IdxCandidateVector.push_back(mx2IdxCandidate);
                ndlarTrkCandidateVector.push_back(npart);

                
              }


        }
          }}} }// particles
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
      auto start_pos = sr->common.ixn.dlp[nixn].part.dlp[npart].start;
      auto end_pos = sr->common.ixn.dlp[nixn].part.dlp[npart].end;


      double dX = (end_pos.x - start_pos.x);
      double dY = (end_pos.y - start_pos.y);
      double dZ = (end_pos.z - start_pos.z);
      double length = TMath::Sqrt(dX * dX + dY * dY + dZ * dZ);
      double dirX = dX / length;
      double dirY = dY / length;
      double dirZ = dZ / length;
      
      double cosLReco =
          dirX * beam_x + dirY * beam_y + dirZ * beam_z;

      double trueCosL = -999;

        
        int trueMatchMult = 0;
        int covarIndex=-1;
                Double_t totalWeight[1000];
                if (biggestMatchIndex>-1){
        int covarIndex=rw.FindBinIndex(sr->mc.nu[biggestMatchIndex].pdg,sr->mc.nu[biggestMatchIndex].E);


        if (covarIndex>-1 && covarIndex<151)
        rw.FillAllThrowsForBin(covarIndex, totalWeight);  

                
                };
       if (goodInteraction || goodOutOfSigDef) {

          if (sr->mc.nu[biggestMatchIndex].iscc &&
              abs(sr->mc.nu[biggestMatchIndex].pdg) == 14) {

            for (int primaries = 0;
                 primaries < sr->mc.nu[biggestMatchIndex].prim.size();
                 primaries++) {
              double Elep = sr->mc.nu[biggestMatchIndex].prim[primaries].p.E;
              auto start_pos =
                  sr->mc.nu[biggestMatchIndex].prim[primaries].start_pos;
              auto end_pos =
                  sr->mc.nu[biggestMatchIndex].prim[primaries].end_pos;
              auto p = sr->mc.nu[biggestMatchIndex].prim[primaries].p;
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
              trueCosL = cosL;
              //  std::cout<<cosL<<std::endl;
            }

            if (goodInteraction) {   
              if (trueCosL>0.91){
              responseCosL->Fill(cosLReco, trueCosL);
              for (int index=0;  index<1000; index++){  if (totalWeight[index]>0.01 &&  totalWeight[index]<20){ responseCosLVec.at(index).Fill(cosLReco, trueCosL, totalWeight[index]);
                   backtrackCosLVec.at(index).Fill(trueCosL,totalWeight[index]);}
             else { responseCosLVec.at(index).Fill(cosLReco, trueCosL, 1);
                   backtrackCosLVec.at(index).Fill(trueCosL,1);
              
              
              }
              }
            }

                  
            }}
            // Fill the response matrix for the track multiplicity
            if (cosLReco > 0.9 && trueCosL > 0.9 &&
                goodInteraction) {
              responseMult->Fill(trackMult, trueMatchMult);
              backtrackTrueMult->Fill(trueMatchMult);
               for (int index=0;  index<1000; index++){ if (totalWeight[index]>0.01 &&  totalWeight[index]<20){ responseMultVec.at(index).Fill(trackMult, trueMatchMult, totalWeight[index]);
                    backtrackMultVec.at(index).Fill(trueMatchMult,totalWeight[index]);}
                else{
                backtrackMultVec.at(index).Fill(trueMatchMult,1); responseMultVec.at(index).Fill(trackMult, trueMatchMult, 1);



                }
                                                      
                }
            }
          }
        
        // Fill the stack histograms of the track multiplicity by truth info
        if (cosLReco > 0.9 && (!goodInteraction || trueCosL < 0.9)) {
          track_multBad->Fill(trackMult);
            for (int index=0;  index<1000; index++){
            if (biggestMatchIndex>-1 && totalWeight[index]>0.01 &&  totalWeight[index]<20)  badRecoMultVec.at(index).Fill(trackMult,totalWeight[index]);
            else             badRecoMultVec.at(index).Fill(trackMult,1);}
        }

        if (!goodInteraction || trueCosL < 0.91){
          badRecoCosL->Fill(cosLReco);
            for (int index=0;  index<1000; index++){
            if (biggestMatchIndex>-1 && totalWeight[index]>0.01 &&  totalWeight[index]<20)
            badRecoCosLVec.at(index).Fill(cosLReco,totalWeight[index]);
            
            else      badRecoCosLVec.at(index).Fill(cosLReco,1);}
            }
        if (cosLReco > 0.9) {
          if (goodInteraction && trueCosL > 0.9 &&
              sr->mc.nu[biggestMatchIndex].iscc) {
            track_multGood->Fill(trackMult);
            

        }
        }
      // Fill the reco information for track multiplicity
      if (cosLReco > 0.9) {
        track_mult->Fill(trackMult);
        for (int index=0;  index<1000; index++){
        
        if (biggestMatchIndex>-1 && totalWeight[index]>0.01 &&  totalWeight[index]<20)
        recoMultVec.at(index).Fill(trackMult,totalWeight[index]);
        
        else       recoMultVec.at(index).Fill(trackMult,1);
        }
        
      } 

          recoCosL->Fill(cosLReco);
        for (int index=0;  index<1000; index++){
                double Phi_u   = rw.GetIntegratedFluxThrows()[index];
            recoCosLVec.at(index).SetTitle(Form("%0.08f",Phi_u));
        if (biggestMatchIndex>-1 && totalWeight[index]>0.01 &&  totalWeight[index]<20)


        recoCosLVec.at(index).Fill(cosLReco,totalWeight[index]);
        
        else       recoCosLVec.at(index).Fill(cosLReco,1);
        }

        
  }}// end for interaction
    }   // end for spills
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

  TProfile* trueCosLProfile=new TProfile("trueCosLProfile","trueCosLProfile",6,edges,"S");
    TProfile* trueMultProfile=new TProfile("trueMultProfile","trueMultProfile",binsMult,edgesMult,"S");

for (int index=0;  index<1000; index++){
trueCosLVec.at(index).Write();
true_multVec.at(index).Write();
recoCosLVec.at(index).Write();
recoMultVec.at(index).Write();
responseCosLVec.at(index).Write();
responseMultVec.at(index).Write();
backtrackCosLVec.at(index).Write();
backtrackMultVec.at(index).Write();
badRecoCosLVec.at(index).Write();
badRecoMultVec.at(index).Write();
for (int i=0; i<7; i++){
    double binCenter=trueCosLProfile->GetXaxis()->GetBinCenter(i);   
    trueCosLProfile->Fill(binCenter,trueCosLVec.at(index).GetBinContent(i));
    binCenter=trueMultProfile->GetXaxis()->GetBinCenter(i);   
    trueMultProfile->Fill(binCenter,true_multVec.at(index).GetBinContent(i));
}
}
  trueMultProfile->Write();
  trueCosLProfile->Write();


  track_mult->Write();
  track_multGood->Write();
  track_multBad->Write();

  responseCosL->Write();

  true_mult->Write();
  true_multTrkOnly->Write();



  responseMult->Write();


  recoCosL->Write();
  badRecoCosL->Write();
  trueCosL->Write();

    
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
