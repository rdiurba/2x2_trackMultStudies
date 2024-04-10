#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TEfficiency.h"
#include <iostream>
#include <fstream>
#include <string>
#include "duneanaobj/StandardRecord/StandardRecord.h" //Ideally, this should be SRProxy.h, but there is an include error for that now. Alternatively, you can use SetBranchStatus function in TreeLoader, but it does not work for the common branch (to do)


double signed3dDistance(double firstPoint1, double firstPoint2, double firstPoint3, double secondPoint1, double secondPoint2, double secondPoint3, TVector3 point){

double denominator = sqrt( (secondPoint2-firstPoint2)*(secondPoint2-firstPoint2) + (secondPoint1-firstPoint1)*(secondPoint1-firstPoint1)+ (secondPoint3-firstPoint3)*(secondPoint3-firstPoint3));

double X1=point.X()-firstPoint1;
double Y1=point.Y()-firstPoint2;
double Z1=point.Z()-firstPoint3;

double X2=point.X()-secondPoint1;
double Y2=point.Y()-secondPoint2;
double Z2=point.Z()-secondPoint3;

double numerator=(X1*Y2-Y1*X2)-(X1*Z2-X2*Z1)+(Y1*Z2-Z1*Y2);

return numerator/denominator;




}

int caf_plotter(std::string input_file_list, std::string output_rootfile){
  double minLength=0; int goodIntNum=0; int badIntNum=0; int goodIntInFidVol=0; int badIntInFidVol=0;
int goodMINERvAMatch=0; int totalMINERvAMatch=0; 
 double offset=0;

  //Define histograms
  TH1D *part_energy_hist = new TH1D("recpart_energy", "Reco particle energy in GeV", 100, 0, 1);

  TH1D *track_length= new TH1D("track_length","track_length",100,0,100);
  TH1D *part_multTrkOnly = new TH1D("part_multTrkOnly","part_multTrkOnly",20,0,20);
  TH1D *part_mult=new TH1D("part_mult","part_mult",20,0,20);
   TH1D* bad_origin=new TH1D("rock","rock",2,0,2);
  TH1D *track_mult= new TH1D("track_mult", "track_mult", 20, 0, 20);
  TH1D *track_multMode= new TH1D("track_multMode", "track_multMode", 20, 0, 20);
  TH1D *track_multQE= new TH1D("track_multQE", "track_multQE", 20, 0, 20);
  TH1D *track_multMEC= new TH1D("track_multMEC", "track_multMEC", 20, 0, 20);
  TH1D *track_multRES= new TH1D("track_multRES", "track_multRES", 20, 0, 20);
  TH1D *track_multCOH= new TH1D("track_multCOH", "track_multCOH", 20, 0, 20);
  TH1D *track_multDIS= new TH1D("track_multDIS", "track_multDIS", 20, 0, 20);
  TH1D *track_multNC= new TH1D("track_multNC", "track_multNC", 20, 0, 20);
  TH1D *track_multRock= new TH1D("track_multRock", "track_multRock", 20, 0, 20);
  TH1D *track_multSec= new TH1D("track_multSec", "track_multSec", 20, 0, 20);

  TH1D *track_multBad=new TH1D("track_multBad","track_multBad",20,0,20);
  TH1D *track_multGood=new TH1D("track_multGood","track_multGood",20,0,20);
  TH1D *trackCorrectness=new TH1D("track_correctness","track_correctness",4,0,4);  
  TH1D *showerCorrectness=new TH1D("shower_corrrectness","shower_correctness",4,0,4);
  TH1D *recoHistVertexY= new TH1D("recoHistVertexY","recoHistVertexY",70,-338,-198); 
TH1D *recoHistVertexX= new TH1D("recoHistVertexX","recoHistVertexX",70,-70,70); 
TH1D *recoHistVertexZ= new TH1D("recoHistVertexZ","recoHistVertexZ",70,1230,1370); 
  TH1D *diffPip=new TH1D("diffPip","diffPip",8,-4,4);
  TH1D *diffPim=new TH1D("diffPim","diffPim",8,-4,4);
  TH1D *diffP=new TH1D("diffP","diffP",8,-4,4);
  TH1D *diffN=new TH1D("diffN","diff",8,-4,4);
  TH1D *recoBacktrackCCAr=new TH1D("recoBacktrackCCAr","recoBacktrackCCAr",2,0,2);
  TH1D *recoBacktrackCCRock=new TH1D("recoBacktrackCCRock","recoBacktrackCCRock",2,0,2);
  TH1D *recoBacktrackCCSec=new TH1D("recoBacktrackCCSec","recoBacktrackCCSec",2,0,2);

  TH1D *recoBacktrackPDGAr=new TH1D("recoBacktrackPDGAr","recoBacktrackPDGAr",40,-20,20);
  TH1D *recoBacktrackPDGRock=new TH1D("recoBacktrackPDGRock","recoBacktrackPDGRock",40,-20,20);
  TH1D *recoBacktrackPDGSec=new TH1D("recoBacktrackPDGSec","recoBacktrackPDGSec",40,-20,20);

 TH1D* recoBacktrackElAr=new TH1D("recoBacktrackElAr","recoBacktrackElAr",50,0,20);
 TH1D* recoBacktrackCoslAr=new TH1D("recoBacktrackCoslAr","recoBackTrackCoslAr",100,-1,1);

 TH1D* selPionE=new TH1D("selPionE","selPionE",20,0,0.1);
 TH1D* selProtonE=new TH1D("selProtonE","selProtonE",20,0,0.1);
 TH1D* selKaonE=new TH1D("selKaonE","selKaonE",20,0,0.1);

 TH1D* selPionDirX=new TH1D("selPionDirX","selPionDirX",20,-1,1);
 TH1D* selPionDirY=new TH1D("selPionDirY","selPionDirY",20,-1,1);
 TH1D* selPionDirZ=new TH1D("selPionDirZ","selPionDirZ",20,-1,1);

 TH1D* selProtonDirX=new TH1D("selProtonDirX","selProtonDirX",20,-1,1);
 TH1D* selProtonDirY=new TH1D("selProtonDirY","selProtonDirY",20,-1,1);
 TH1D* selProtonDirZ=new TH1D("selProtonDirZ","selProtonDirZ",20,-1,1);

 TH1D* selKaonDirX=new TH1D("selKaonDirX","selKaonDirX",20,-1,1);
 TH1D* selKaonDirY=new TH1D("selKaonDirY","selKaonDirY",20,-1,1);
 TH1D* selKaonDirZ=new TH1D("selKaonDirZ","selKaonDirZ",20,-1,1);


 TH1D* selProtonXZ=new TH1D("selProtonXZ","selProtonXZ",20,-1,1);
 TH1D* selProtonYZ=new TH1D("selProtonYZ","selProtonYZ",20,-1,1);
 TH1D* selProtonXY=new TH1D("selProtonXY","selProtonXY",20,-1,1);

 TH1D* selPionXZ=new TH1D("selPionXZ","selPionXZ",20,-1,1);
 TH1D* selPionYZ=new TH1D("selPionYZ","selPionYZ",20,-1,1);
 TH1D* selPionXY=new TH1D("selPionXY","selPionXY",20,-1,1);

 TH1D* selKaonXZ=new TH1D("selKaonXZ","selKaonXZ",20,-1,1);
 TH1D* selKaonYZ=new TH1D("selKaonYZ","selKaonYZ",20,-1,1);
 TH1D* selKaonXY=new TH1D("selKaonXY","selKaonXY",20,-1,1);

 TH1D* truePionEWithRecoInt=new TH1D("truePionEWithRecoInt","truePionEWithRecoInt",20,0,0.1);
 TH1D* trueProtonEWithRecoInt=new TH1D("trueProtonEWithRecoInt","trueProtonEWithRecoInt",20,0,0.1);
 TH1D* trueKaonEWithRecoInt=new TH1D("trueKaonEWithRecoInt","trueKaonEWithRecoInt",20,0,0.1);

 TH1D* trueProtonWithRecoIntDirX=new TH1D("trueProtonWithRecoIntDirX","trueProtonWithRecoIntDirX",20,-1,1);
 TH1D* trueProtonWithRecoIntDirY=new TH1D("trueProtonWithRecoIntDirY","trueProtonWithRecoIntDirY",20,-1,1);
 TH1D* trueProtonWithRecoIntDirZ=new TH1D("trueProtonWithRecoIntDirZ","trueProtonWithRecoIntDirZ",20,-1,1);

 TH1D* truePionWithRecoIntDirX=new TH1D("truePionWithRecoIntDirX","truePionWithRecoIntDirX",20,-1,1);
 TH1D* truePionWithRecoIntDirY=new TH1D("truePionWithRecoIntDirY","truePionWithRecoIntDirY",20,-1,1);
 TH1D* truePionWithRecoIntDirZ=new TH1D("truePionWithRecoIntDirZ","truePionWithRecoIntDirZ",20,-1,1);


 TH1D* trueKaonWithRecoIntDirX=new TH1D("trueKaonWithRecoIntDirX","trueKaonWithRecoIntDirX",20,-1,1);
 TH1D* trueKaonWithRecoIntDirY=new TH1D("trueKaonWithRecoIntDirY","trueKaonWithRecoIntDirY",20,-1,1);
 TH1D* trueKaonWithRecoIntDirZ=new TH1D("trueKaonWithRecoIntDirZ","trueKaonWithRecoIntDirZ",20,-1,1);

 TH1D* trueProtonWithRecoIntXY=new TH1D("trueProtonWithRecoIntXY","trueProtonWithRecoIntXY",20,-1,1);
 TH1D* trueProtonWithRecoIntXZ=new TH1D("trueProtonWithRecoIntXZ","trueProtonWithRecoIntXZ",20,-1,1);
 TH1D* trueProtonWithRecoIntYZ=new TH1D("trueProtonWithRecoIntYZ","trueProtonWithRecoIntYZ",20,-1,1);

 TH1D* truePionWithRecoIntXY=new TH1D("truePionWithRecoIntXY","truePionWithRecoIntXY",20,-1,1);
 TH1D* truePionWithRecoIntXZ=new TH1D("truePionWithRecoIntXZ","truePionWithRecoIntXZ",20,-1,1);
 TH1D* truePionWithRecoIntYZ=new TH1D("truePionWithRecoIntYZ","truePionWithRecoIntYZ",20,-1,1);

 TH1D* trueKaonWithRecoIntXY=new TH1D("trueKaonWithRecoIntXY","trueKaonWithRecoIntXY",20,-1,1);
 TH1D* trueKaonWithRecoIntXZ=new TH1D("trueKaonWithRecoIntXZ","trueKaonWithRecoIntXZ",20,-1,1);
 TH1D* trueKaonWithRecoIntYZ=new TH1D("trueKaonWithRecoIntYZ","trueKaonWithRecoIntYZ",20,-1,1);



  TH1D *true_mult=new TH1D("true_mult","true_mult",20,0,20);
  TH1D *true_multTrkOnly=new TH1D("true_multTrkOnly","true_multTrkOnly",20,0,20);
  TH1D *true_multGENIE=new TH1D("true_multGENIE","true_multGENIE",20,0,20);

  TH1D *matchTrue_mult=new TH1D("matchTrue_mult","matchTrue_mult",20,0,20);
  TH1D *matchTrue_multTrkOnly=new TH1D("matchTrue_multTrkOnly","matchTrue_multTrkOnly",20,0,20);
  TH1D *matchTrue_multGENIE=new TH1D("matchTrue_multGENIE","matchTrue_multGENIE",20,0,20);
  TH1D *matchHistEl=new TH1D("matchHistEl","matchHistEl",50,0,20);
  TH1D *matchHistCosl=new TH1D("matchHistCosl","matchHistCosl",20,-1,1);

  TH1D *histEl=new TH1D("histEl","histEl",50,0,20);
  TH1D *histCosl=new TH1D("histCosl","histCosl",100,-1,1);

  TH2D *recoVertex2DNoCuts=new TH2D("recoVertex2DNoCuts","recoVertex2DNoCuts",70,-70,70,70,-70,70);
  TH2D *recoVertex2D=new TH2D("recoVertex2D","recoVertex2D",60,-60,60,60,-60,60);
  TH2D *recoVertex2DBadYZ=new TH2D("recoVertex2DBadYZ","recoVertex2DBadYZ",70,-70,70,70,-70,70);


  TH2D *recoVertex2DBad=new TH2D("recoVertex2DBad","recoVertex2DBad",70,-70,70,70,-70,70);
  TH2D *trueVertex2DBad=new TH2D("trueVertex2DBad","trueVertex2DBad",60,-60,60,60,-60,60);



  TH1D *longest_track=new TH1D("longest_track","longest_track",100,0,100);


  //Give an input list
  std::ifstream caf_list(input_file_list.c_str());

  //Check if input list is present
  if(!caf_list.is_open()){
	std::cerr << Form("File %s not found", input_file_list.c_str()) << std::endl;
	return 1;
  }

  //Add files to CAF chain from input list
  std::string tmp;
  TChain *caf_chain = new TChain("cafTree");

  while(caf_list >> tmp){
	caf_chain->Add(tmp.c_str());
	std::cout << Form("Adding File %s", tmp.c_str()) << std::endl;
  }

  //Check if CAF tree is present
  if(!caf_chain){
	std::cerr << Form("There is no tree in %s", tmp.c_str()) << std::endl;
	return 1;
  }

  long Nentries = caf_chain->GetEntries();
  std::cout << Form("Total number of spills = %ld", Nentries) << std::endl;

  //Define Standard Record and link it to the CAF tree branch "rec"
  auto sr = new caf::StandardRecord;
  caf_chain->SetBranchAddress("rec", &sr);

  std::vector<double> goodEvents={1,7,24,30,50,122,131,142,143,163,210,226,267,272,275,276,339,378,395,442,497,502,517,530,534,597,599,617,655,699,707,787,845,849,898,921,942,948,966,968,986,998,1048,1056,1059,1076,1077,1079,1154,1168,1207};


  bool skipEvent=false;
  for(long n = 0; n < Nentries; ++n){ //Loop over spills/triggers
//if(std::find(goodEvents.begin(), goodEvents.end(), n) != goodEvents.end()){}
//else continue;
        double trueVtxX=-9999; double trueVtxY=-999; double trueVtxZ=-9999;
        skipEvent=false;
	if(n%10000 == 0) std::cout << Form("Processing trigger %ld of %ld", n, Nentries) << std::endl;
	caf_chain->GetEntry(n); //Get spill from tree
       bool hasANeutrino=false;

       for(long unsigned ntrue=0; ntrue<sr->mc.nu.size(); ntrue++){ 
       auto vertex=sr->mc.nu[ntrue].vtx;
       trueVtxX=vertex.x; trueVtxY=vertex.y; trueVtxZ=vertex.z;
      
       auto truePrimary=sr->mc.nu[ntrue].prim;
       int truePart=0; int truePartTrkOnly=0; int truePartNoG4=0;
       if (sr->mc.nu[ntrue].targetPDG!=1000180400) continue;
       if (sr->mc.nu[ntrue].iscc==false) continue;
       if (abs(sr->mc.nu[ntrue].pdg)!=14) continue; 
       if (abs(sr->mc.nu[ntrue].vtx.x)>53 || abs(sr->mc.nu[ntrue].vtx.x)<12) continue;
       if (abs(sr->mc.nu[ntrue].vtx.z)>53 || abs(sr->mc.nu[ntrue].vtx.z)<12) continue;
       if (abs(sr->mc.nu[ntrue].vtx.y)>53) continue;    
  //std::cout<<"True Vertex Location of Argon Interaction: "<<trueVtxX<<","<<trueVtxY<<","<<trueVtxZ<<std::endl;
       int npip=sr->mc.nu[ntrue].npip;
       int npim=sr->mc.nu[ntrue].npim;
       int np=sr->mc.nu[ntrue].nproton;
       int nneutrons=sr->mc.nu[ntrue].nneutron;

       int npipPrimaries=0;
       int npimPrimaries=0;
       int npPrimaries=0;
       int nneutronsPrimaries=0;
       for (int k=0; k<sr->mc.nu[ntrue].prim.size(); k++){ if(sr->mc.nu[ntrue].prim[k].pdg==211) npipPrimaries++; if(sr->mc.nu[ntrue].prim[k].pdg==2212) npPrimaries++;  if(sr->mc.nu[ntrue].prim[k].pdg==-211) npimPrimaries++;  if(sr->mc.nu[ntrue].prim[k].pdg==2112) nneutronsPrimaries++;}

       diffPip->Fill(npipPrimaries-npip);
       diffPim->Fill(npimPrimaries-npim);
       diffP->Fill(npPrimaries-np);
       diffN->Fill(nneutronsPrimaries-nneutrons);

      double  Elep=-999; double cosL=-999;
       for(long unsigned primaries=0; primaries<truePrimary.size(); primaries++){
       int pdg=truePrimary[primaries].pdg;
       double E=truePrimary[primaries].p.E;
       double px=truePrimary[primaries].p.px; double py=truePrimary[primaries].p.py; double pz=truePrimary[primaries].p.pz;
       double totP=TMath::Sqrt(px*px+py*py+pz*pz);

       if ((abs(pdg)==13 || abs(pdg)==2212 || abs(pdg)==211 || abs(pdg)==321)){
        truePartNoG4++;
       
       auto start_pos=sr->mc.nu[ntrue].prim[primaries].start_pos;
       auto end_pos=sr->mc.nu[ntrue].prim[primaries].end_pos;
       if (std::isnan(start_pos.z)){  continue;} 
       double dX=abs(end_pos.x-start_pos.x);
       double dY=abs(end_pos.y-start_pos.y);
       double dZ=abs(end_pos.z-start_pos.z);
       double length=TMath::Sqrt(dX*dX+dY*dY+dZ*dZ);
       if (sr->mc.nu[ntrue].iscc && abs(sr->mc.nu[ntrue].pdg)==14 && abs(pdg)==13){
     
       Elep=sr->mc.nu[ntrue].prim[primaries].p.E;
       cosL=dZ/length;
       histEl->Fill(Elep);
       histCosl->Fill(cosL);


       }
       //if (cosL<0.9) continue;
       if (length>5){ truePartTrkOnly++;
 

       }
       }
       truePart++;
       }
       true_mult->Fill(truePart); if (truePartTrkOnly>0) true_multTrkOnly->Fill(truePartTrkOnly);    true_multGENIE->Fill(truePartNoG4);
       hasANeutrino=true;  
  }
    //for(long unsigned nixn=0; nixn<sr  

  //if (hasANeutrino==false) continue;

       //if (sr->common.ixn.dlp.size()>1) continue; 
        int rock=0;
        bool goodInteraction=false;   
        std::vector<int> trueInteractionIndex;
        std::vector<int> primaryTrkIndex;
	
	for(long unsigned nixn = 0; nixn < sr->common.ixn.dlp.size(); nixn++){
      	  bool oneContained=false; bool oneNotContained=false; goodInteraction=false;
          double biggestMatch=-999; int biggestMatchIndex=-999; double maxDotProductDS=-999; double maxDotProductUS=-999;
            int maxEventPar=-999; int maxEventTyp=-9999; int maxEventIxn=-999;
           for (int ntruth=0; ntruth<sr->common.ixn.dlp[nixn].truth.size(); ntruth++){
          
          if (biggestMatch<sr->common.ixn.dlp[nixn].truthOverlap.at(ntruth)){
          
          biggestMatch=sr->common.ixn.dlp[nixn].truthOverlap.at(ntruth);
          biggestMatchIndex=sr->common.ixn.dlp[nixn].truth.at(ntruth);
           //std::cout<<"Biggest Match: "<<biggestMatch<<std::endl;


		}
		}
          //std::cout<<biggestMatchIndex<<","<<sr->common.ixn.dlp[nixn].truth.at(ntruth)<<std::endl;
          if (sr->mc.nu[biggestMatchIndex].id>1E9){ rock=1; 
	
	}
	
          double recoVertexX=sr->common.ixn.dlp[nixn].vtx.x; double recoVertexY=sr->common.ixn.dlp[nixn].vtx.y; double recoVertexZ=sr->common.ixn.dlp[nixn].vtx.z;
        


	  trueVtxX=sr->mc.nu[biggestMatchIndex].vtx.x; trueVtxY=sr->mc.nu[biggestMatchIndex].vtx.y; trueVtxZ=sr->mc.nu[biggestMatchIndex].vtx.z;
           //if (rock==1) std::cout<<" Checking rock vertex is now a number: "<<trueVtxZ<<std::endl;
//        if (rock==0) std::cout<<recoVertexX<<","<<recoVertexY<<","<<recoVertexZ<<","<<trueVtxX<<","<<trueVtxY<<","<<trueVtxZ<<std::endl;
	  if (abs(sr->mc.nu[biggestMatchIndex].pdg)==14 && abs(sr->mc.nu[biggestMatchIndex].iscc)==true  && sr->mc.nu[biggestMatchIndex].targetPDG==1000180400 &&  abs(recoVertexX-trueVtxX)<5 && abs(recoVertexY-trueVtxY)<5 && abs(recoVertexZ-trueVtxZ)<5 && abs(trueVtxX)>12 && abs(trueVtxX)<53 && abs(trueVtxZ)>12 && abs(trueVtxZ)<53 && abs(trueVtxY)<53) goodInteraction=true;
         recoHistVertexY->Fill(sr->common.ixn.dlp[nixn].vtx.y); recoHistVertexX->Fill(sr->common.ixn.dlp[nixn].vtx.x); recoHistVertexZ->Fill(sr->common.ixn.dlp[nixn].vtx.z);
         if (goodInteraction==true) goodIntNum++; else badIntNum++;
	 if (abs(abs(sr->common.ixn.dlp[nixn].vtx.x)-33)<1.5 || abs(sr->common.ixn.dlp[nixn].vtx.x)>53 || abs(sr->common.ixn.dlp[nixn].vtx.x)<12 || abs(sr->common.ixn.dlp[nixn].vtx.y)>53 || abs(sr->common.ixn.dlp[nixn].vtx.z)<12 || abs(sr->common.ixn.dlp[nixn].vtx.z)>53)  continue; 
          trueInteractionIndex.push_back(biggestMatchIndex);

          if (goodInteraction==true) goodIntInFidVol++; else badIntInFidVol++;
                recoVertex2D->Fill(sr->common.ixn.dlp[nixn].vtx.x,sr->common.ixn.dlp[nixn].vtx.z);

 
       //std::cout<<"counting protons and kaons"<<std::endl;

      //std::cout<<"Done counting"<<std::endl;
       




                int partMult=0;
                int partMultTrkOnly=0;
                bool hasAtLeastOneProton=false;
                bool hasAtLeastOneMuon=false;
                bool hasAtLeastOne=false;
                //std::cout<<"Interaction: "<<nixn<<std::endl;
                int correctTrack=2;
                int correctShower=2;
		 double longestTrk=-9999;
    int trackMult=0;  int trackMultExit=0; int minervaTracks=0; int minervaThrough=0; 
    	for(long unsigned npart=0; npart < sr->common.ixn.dlp[nixn].part.dlp.size(); npart++){ //loop over particles
 	//std::cout<<"Containment variable: "<<sr->common.ixn.dlp[nixn].part.dlp[npart].contained<<std::endl;
                if (!sr->common.ixn.dlp[nixn].part.dlp[npart].primary) continue; 
               int pdg=sr->common.ixn.dlp[nixn].part.dlp[npart].pdg;
	if ( (abs(pdg)==2212 || abs(pdg)==13 || abs(pdg)==211 || abs(pdg)==321)){
	if (sr->common.ixn.dlp[nixn].part.dlp[npart].primary==true) trackCorrectness->Fill(correctTrack);

       auto start_pos=sr->common.ixn.dlp[nixn].part.dlp[npart].start;
       auto end_pos=sr->common.ixn.dlp[nixn].part.dlp[npart].end;

	double diffVertexdZ=abs(start_pos.z-sr->common.ixn.dlp[nixn].vtx.z);
	double diffVertexdX=abs(start_pos.x-sr->common.ixn.dlp[nixn].vtx.x);
	double diffVertexdY=abs(start_pos.y-sr->common.ixn.dlp[nixn].vtx.y);
	double diffVertex=TMath::Sqrt(diffVertexdZ*diffVertexdZ+diffVertexdY*diffVertexdY+diffVertexdX*diffVertexdX);
	if (diffVertex>5) continue;
      
               double dX=(end_pos.x-start_pos.x);
               double dY=(end_pos.y-start_pos.y);
               double dZ=(end_pos.z-start_pos.z);
               double length=TMath::Sqrt(dX*dX+dY*dY+dZ*dZ);
                double dirX=dX/length; double dirY=dY/length; double dirZ=dZ/length;

                if (dirZ<0){ dirZ=-dirZ; dirX=-dirX; dirY=-dirY; auto temp=start_pos; end_pos=start_pos; end_pos=temp;}
        if (std::isnan(start_pos.z)) length=-999;
	if (length>longestTrk) longestTrk=length;
        if (sr->common.ixn.dlp[nixn].part.dlp[npart].primary==true && length>5){   
		partMult++;
                trackMult++;
           if    ((start_pos.z)>58 || (end_pos.z)>58) trackMultExit++; 
		int maxPartMinerva=-999; int maxTypeMinerva=-999;  
                  int maxIxnMinerva=-999;
		int maxPartMinervaUS=-999; int maxTypeMinervaUS=-999;      
		if ((abs(start_pos.z)>58 || abs(end_pos.z)>58) ){ int minervaPass=0;
		double dotProductDS=-999; double deltaExtrapYUS=-999; double deltaExtrapY=-999; double dotProductUS=-999; double deltaExtrapX=-999; double deltaExtrapXUS=-999;	
	 	for(int i=0; i<sr->nd.minerva.ixn.size(); i++){

		for (int j=0; j<sr->nd.minerva.ixn[i].ntracks; j++){
		    double dir_z=sr->nd.minerva.ixn[i].tracks[j].dir.z;
		double end_z=sr->nd.minerva.ixn[i].tracks[j].end.z;
		double start_z=sr->nd.minerva.ixn[i].tracks[j].start.z;
		double end_x=sr->nd.minerva.ixn[i].tracks[j].end.x;
		double start_x=sr->nd.minerva.ixn[i].tracks[j].start.x;

		double end_y=sr->nd.minerva.ixn[i].tracks[j].end.y;
		double start_y=sr->nd.minerva.ixn[i].tracks[j].start.y;
                



		if (start_z>0 && ((start_pos.z)>58 || (end_pos.z)>58) ){
		int truthPart=sr->nd.minerva.ixn[i].tracks[j].truth[0].part;
		double dXMnv=(sr->nd.minerva.ixn[i].tracks[j].end.x-sr->nd.minerva.ixn[i].tracks[j].start.x);
		double dYMnv=(sr->nd.minerva.ixn[i].tracks[j].end.y-sr->nd.minerva.ixn[i].tracks[j].start.y);
		double dZMnv=(sr->nd.minerva.ixn[i].tracks[j].end.z-sr->nd.minerva.ixn[i].tracks[j].start.z);
		double lengthMinerva=TMath::Sqrt(dXMnv*dXMnv+dYMnv*dYMnv+dZMnv*dZMnv);
	        if (lengthMinerva<10) continue;
          	double dirXMinerva=dXMnv/lengthMinerva;
		double dirYMinerva=dYMnv/lengthMinerva;
		double dirZMinerva=dZMnv/lengthMinerva;
		double dotProduct=dirXMinerva*dirX+dirYMinerva*dirY+dirZ*dirZMinerva;
		double extrapdZ=start_z-end_pos.z;
		double extrapY=dirY/dirZ*(extrapdZ)+end_pos.y-start_y;
		double extrapX=dirX/dirZ*(extrapdZ)+end_pos.x-start_x;
		double diffExtrap=TMath::Sqrt(TMath::Power(extrapY-start_y,2));


		if (dotProductDS<dotProduct && abs(extrapY)<10 && abs(TMath::ACos(dirXMinerva)-TMath::ACos(dirX))<0.06 && abs(TMath::ACos(dirYMinerva)-TMath::ACos(dirY))<0.06 && abs(extrapX)<10){ dotProductDS=dotProduct;
		deltaExtrapY=extrapY;
		deltaExtrapX=extrapX;
		maxPartMinerva=sr->nd.minerva.ixn[i].tracks[j].truth[0].part;
		maxTypeMinerva=sr->nd.minerva.ixn[i].tracks[j].truth[0].type;
                maxIxnMinerva=sr->nd.minerva.ixn[i].tracks[j].truth[0].ixn;	
		if (end_z>300) minervaPass=1;
                }}

		if (start_z<0 && end_z>0 && ( (start_pos.z<-58 && end_pos.z>58) || (start_pos.z>58 && end_pos.z<-58))){
		int truthPart=sr->nd.minerva.ixn[i].tracks[j].truth[0].part;
		double dXMnv=(sr->nd.minerva.ixn[i].tracks[j].end.x-sr->nd.minerva.ixn[i].tracks[j].start.x);
		double dYMnv=(sr->nd.minerva.ixn[i].tracks[j].end.y-sr->nd.minerva.ixn[i].tracks[j].start.y);
		double dZMnv=(sr->nd.minerva.ixn[i].tracks[j].end.z-sr->nd.minerva.ixn[i].tracks[j].start.z);
		double lengthMinerva=TMath::Sqrt(dXMnv*dXMnv+dYMnv*dYMnv+dZMnv*dZMnv);
		double dirXMinerva=dXMnv/lengthMinerva;
		double dirYMinerva=dYMnv/lengthMinerva;
		double dirZMinerva=dZMnv/lengthMinerva;
		double dotProduct=dirXMinerva*dirX+dirYMinerva*dirY+dirZ*dirZMinerva;

		double extrapdZUS=end_z-end_pos.z;
		double extrapYUS=dirY/dirZ*(extrapdZUS)+end_pos.y-end_y;
		double extrapXUS=dirX/dirZ*(extrapdZUS)+end_pos.x-end_x;
		
                //double diffExtrap=TMath::Sqrt(TMath::Power(extrapY-end_y,2));
		if (dotProductUS<dotProduct  && abs(extrapYUS)<10 && abs(extrapXUS)<10 && abs(TMath::ACos(dirXMinerva)-TMath::ACos(dirX))<0.06 && abs(TMath::ACos(dirYMinerva)-TMath::ACos(dirY))<0.06) dotProductUS=dotProduct;
		deltaExtrapYUS=extrapYUS;
                deltaExtrapXUS=extrapXUS;
		maxPartMinervaUS=sr->nd.minerva.ixn[i].tracks[j].truth[0].part;
		maxTypeMinervaUS=sr->nd.minerva.ixn[i].tracks[j].truth[0].type;	
//if (maxPartMinervaUS==maxPartNumber) std::cout<<extrapdZ<<","<<dirY/dirZ*(extrapdZ)<<","<<double(dirX/dirZ*(extrapdZ))<<","<<end_y<<","<<end_x<<std::endl;


		}

		
		}} // Minerva
	if (dotProductDS>maxDotProductDS){ maxDotProductDS=dotProductDS;
        maxEventPar=maxPartMinerva;
        maxEventTyp=maxTypeMinerva;
        maxEventIxn=maxIxnMinerva;
	}
	 if (dotProductDS>0.9975){  minervaTracks++; if (minervaPass==1){ minervaThrough++; 
/*		if (maxPartNumber==maxPartMinerva && maxTypeMinerva==maxTypeNumber){ goodMINERvAMatch++; deltaXGood->Fill(deltaExtrapX); deltaYGood->Fill(deltaExtrapY);}
  else{ deltaYBad->Fill(deltaExtrapY); deltaXBad->Fill(deltaExtrapX);}
	      totalMINERvAMatch++;       
*/	
          //std::cout<<maxIxnMinerva<<","<<biggestMatchIndex<<","<<maxTypeMinerva<<","<<maxPartMinerva<<","<<maxTypeNumber<<","<<maxPartNumber<<std::endl;
         if (maxTypeMinerva<3 && maxPartNumber==maxPartMinerva && maxTypeMinerva==maxTypeNumber) std::cout<<sr->mc.nu[maxIxnMinerva].prim[maxPartMinerva].pdg<<std::endl;


	 }}
	if (dotProductUS>maxDotProductUS) maxDotProductUS=dotProductUS;
	 if (dotProductUS>0.9975){   
/*		if (maxPartNumber==maxPartMinervaUS && maxTypeMinervaUS==maxTypeNumber){ deltaYGoodUS->Fill(deltaExtrapYUS); deltaXGoodUS->Fill(deltaExtrapXUS); goodMINERvAMatchUS++;} 
              else{ deltaYBadUS->Fill(deltaExtrapYUS); deltaXBadUS->Fill(deltaExtrapXUS); }
       totalMINERvAMatchUS++;       
*/	 }
  //oneNotContained=true;
		} // exiting particles
//else oneContained=true;
		} // primary of particle greater than 5
		  } } // particles
	if (minervaTracks>1 || minervaThrough<1 ||  maxDotProductDS<0.9975) continue;
//	if (maxDotProductUS>0.9999) continue;


 //std::cout<<"counting protons and kaons"<<std::endl;
      if (sr->mc.nu[biggestMatchIndex].prim.size()<500 && sr->mc.nu[biggestMatchIndex].targetPDG==1000180400){ 
       for (int primaries=0; primaries<sr->mc.nu[biggestMatchIndex].prim.size(); primaries++){ 
      
      auto start_pos=sr->mc.nu[biggestMatchIndex].prim[primaries].start_pos;
      auto end_pos=sr->mc.nu[biggestMatchIndex].prim[primaries].end_pos;
      if (std::isnan(start_pos.z)) continue;
      //std::cout<<"Got a position"<<std::endl; 
	      
               double dX=(end_pos.x-start_pos.x);
               double dY=(end_pos.y-start_pos.y);
               double dZ=(end_pos.z-start_pos.z);
               double length=TMath::Sqrt(dX*dX+dY*dY+dZ*dZ);
                double dirX=dX/length; double dirY=dY/length; double dirZ=dZ/length;
      if(abs(sr->mc.nu[biggestMatchIndex].prim[primaries].pdg)==2212){ trueProtonEWithRecoInt->Fill(sr->mc.nu[biggestMatchIndex].prim[primaries].p.E-0.938272); trueProtonWithRecoIntDirZ->Fill(dirZ); trueProtonWithRecoIntDirX->Fill(dirX); trueProtonWithRecoIntDirY->Fill(dirY);} 
          if(abs(sr->mc.nu[biggestMatchIndex].prim[primaries].pdg)==211){ truePionEWithRecoInt->Fill(sr->mc.nu[biggestMatchIndex].prim[primaries].p.E-0.13957);
truePionWithRecoIntDirZ->Fill(dirZ); truePionWithRecoIntDirX->Fill(dirX); truePionWithRecoIntDirY->Fill(dirY);


}
      if(abs(sr->mc.nu[biggestMatchIndex].prim[primaries].pdg)==321){ trueKaonEWithRecoInt->Fill(sr->mc.nu[biggestMatchIndex].prim[primaries].p.E-0.493677);
      trueKaonWithRecoIntDirZ->Fill(dirZ); trueKaonWithRecoIntDirX->Fill(dirX); trueKaonWithRecoIntDirY->Fill(dirY);}
      }}


    		for(long unsigned npart=0; npart < sr->common.ixn.dlp[nixn].part.dlp.size(); npart++){ //loop over particles
 	//std::cout<<"Containment variable: "<<sr->common.ixn.dlp[nixn].part.dlp[npart].contained<<std::endl;
                if (!sr->common.ixn.dlp[nixn].part.dlp[npart].primary) continue; 
               int pdg=sr->common.ixn.dlp[nixn].part.dlp[npart].pdg;
if ( abs(pdg)==22 || abs(pdg)==11 || abs(pdg)==111){
auto truthSize=sr->common.ixn.dlp[nixn].part.dlp[npart].truth.size();
double maxPartTruthOverlap=0;

for (int backTrack=0; backTrack<truthSize; backTrack++){
int parType=sr->common.ixn.dlp[nixn].part.dlp[npart].truth[backTrack].type;
int partNumber=sr->common.ixn.dlp[nixn].part.dlp[npart].truth[backTrack].part;
int interactionNumber=sr->common.ixn.dlp[nixn].part.dlp[npart].truth[backTrack].ixn;

double  partTruthOverlap=sr->common.ixn.dlp[nixn].part.dlp[npart].truthOverlap[backTrack];
if (maxPartTruthOverlap<partTruthOverlap){maxPartTruthOverlap=partTruthOverlap;}
}

for (int backTrack=0; backTrack<truthSize; backTrack++){
int parType=sr->common.ixn.dlp[nixn].part.dlp[npart].truth[backTrack].type;
int partNumber=sr->common.ixn.dlp[nixn].part.dlp[npart].truth[backTrack].part;
int interactionNumber=sr->common.ixn.dlp[nixn].part.dlp[npart].truth[backTrack].ixn;
double  partTruthOverlap=sr->common.ixn.dlp[nixn].part.dlp[npart].truthOverlap[backTrack];

int backtracked=-9999;
if (partTruthOverlap>maxPartTruthOverlap){ 
if (parType<3){ backtracked=sr->mc.nu[interactionNumber].prim[partNumber].pdg; primaryTrkIndex.push_back(partNumber);}
else backtracked=sr->mc.nu[interactionNumber].sec[partNumber].pdg;
correctShower=0;
if (parType>2) correctShower=3;
else if (backtracked==11 || backtracked==22 || backtracked==111 || backtracked==2112){
correctShower=1;
}

}
else backtracked=-1;
}
if (sr->common.ixn.dlp[nixn].part.dlp[npart].primary==true) showerCorrectness->Fill(correctShower);







}
int maxPartNumber=-999; int maxTypeNumber=-999;
if ( (abs(pdg)==2212 || abs(pdg)==13 || abs(pdg)==211 || abs(pdg)==321)){
auto truthSize=sr->common.ixn.dlp[nixn].part.dlp[npart].truth.size();
double maxPartTruthOverlap=0;
for (int backTrack=0; backTrack<truthSize; backTrack++){
int parType=sr->common.ixn.dlp[nixn].part.dlp[npart].truth[backTrack].type;
int partNumber=sr->common.ixn.dlp[nixn].part.dlp[npart].truth[backTrack].part;
int interactionNumber=sr->common.ixn.dlp[nixn].part.dlp[npart].truth[backTrack].ixn;
double  partTruthOverlap=sr->common.ixn.dlp[nixn].part.dlp[npart].truthOverlap[backTrack];
//std::cout<<"Truth has overlap of: "<<partTruthOverlap<<std::endl;
if (maxPartTruthOverlap<partTruthOverlap){maxPartTruthOverlap=partTruthOverlap;}
}

for (int backTrack=0; backTrack<truthSize; backTrack++){
int parType=sr->common.ixn.dlp[nixn].part.dlp[npart].truth[backTrack].type;
int partNumber=sr->common.ixn.dlp[nixn].part.dlp[npart].truth[backTrack].part;
int interactionNumber=sr->common.ixn.dlp[nixn].part.dlp[npart].truth[backTrack].ixn;
double  partTruthOverlap=sr->common.ixn.dlp[nixn].part.dlp[npart].truthOverlap[backTrack];

int backtracked=-9999;
if (partTruthOverlap>0.5){  
maxPartNumber=partNumber;
maxTypeNumber=parType;
if (parType<3) backtracked=sr->mc.nu[interactionNumber].prim[partNumber].pdg;
else backtracked=sr->mc.nu[interactionNumber].sec[partNumber].pdg;
if (parType<3){
      auto start_pos=sr->mc.nu[interactionNumber].prim[partNumber].start_pos;
      if (std::isnan(start_pos.z)) continue;

      auto end_pos=sr->mc.nu[interactionNumber].prim[partNumber].end_pos;
	      
               double dX=(end_pos.x-start_pos.x);
               double dY=(end_pos.y-start_pos.y);
               double dZ=(end_pos.z-start_pos.z);
               double length=TMath::Sqrt(dX*dX+dY*dY+dZ*dZ);
                double dirX=dX/length; double dirY=dY/length; double dirZ=dZ/length;

if (sr->mc.nu[biggestMatchIndex].targetPDG==1000180400 && backtracked==2212 && parType<3){ selProtonE->Fill(sr->mc.nu[interactionNumber].prim[partNumber].p.E-0.938272);
selProtonDirX->Fill(dirX); selProtonDirZ->Fill(dirZ); selProtonDirY->Fill(dirY);

}
if (sr->mc.nu[biggestMatchIndex].targetPDG==1000180400 && abs(backtracked)==211 && parType<3){ selPionE->Fill(sr->mc.nu[interactionNumber].prim[partNumber].p.E-0.13957039);
selPionDirX->Fill(dirX); selPionDirZ->Fill(dirZ); selPionDirY->Fill(dirY);

}
if (sr->mc.nu[biggestMatchIndex].targetPDG==1000180400 && abs(backtracked)==321 && parType<3){ selKaonE->Fill(sr->mc.nu[interactionNumber].prim[partNumber].p.E-0.493677);
selKaonDirX->Fill(dirX); selKaonDirZ->Fill(dirZ); selKaonDirY->Fill(dirY);

}
}
correctTrack=0;
//std::cout<<"Backtrackilng for "<<parType<<" : "<<pdg<<","<<backtracked<<std::endl;
if (parType>2) correctTrack=3;
else if (backtracked==11 || backtracked==22 || backtracked==111 || backtracked==2112){
correctTrack=1;
}



}
if (partTruthOverlap<0.5) correctTrack=2;
else backtracked=-1; 
}}
}

       if (trackMult>0 && trackMultExit>0 /*&& oneContained && oneNotContained*/){  
       if (goodInteraction){ recoBacktrackCCAr->Fill(sr->mc.nu[biggestMatchIndex].iscc);
recoBacktrackPDGAr->Fill(sr->mc.nu[biggestMatchIndex].pdg);
       if (sr->mc.nu[biggestMatchIndex].iscc && abs(sr->mc.nu[biggestMatchIndex].pdg)==14){
       for (int primaries=0; primaries<sr->mc.nu[biggestMatchIndex].prim.size(); primaries++){
       if (abs(sr->mc.nu[biggestMatchIndex].prim[primaries].pdg)!=13) continue;
       double Elep=sr->mc.nu[biggestMatchIndex].prim[primaries].p.E;
            auto start_pos=sr->mc.nu[biggestMatchIndex].prim[primaries].start_pos;
       auto end_pos=sr->mc.nu[biggestMatchIndex].prim[primaries].end_pos;
       /*if (std::isnan(start_pos.z)) continue;
       double dX=abs(end_pos.x-start_pos.x);
       double dY=abs(end_pos.y-start_pos.y);
       double dZ=abs(end_pos.z-start_pos.z);
       double length=TMath::Sqrt(dX*dX+dY*dY+dZ*dZ);
       double cosL=dZ/length;
//if (goodInteraction)       std::cout<<biggestMatch<<","<<sr->common.ixn.dlp[nixn].vtx.z<<","<<sr->mc.nu[biggestMatchIndex].vtx.z<<","<<start_pos.z<<","<<end_pos.z<<std::endl;
       recoBacktrackElAr->Fill(Elep);
       recoBacktrackCoslAr->Fill(cosL);*/
       }
      
       }










       if (maxEventTyp<3){
       auto start_pos=sr->mc.nu[maxEventIxn].prim[maxEventPar].start_pos;
       auto end_pos=sr->mc.nu[maxEventIxn].prim[maxEventPar].end_pos;
       if (!std::isnan(start_pos.z)){
       double dX=abs(end_pos.x-start_pos.x);
       double dY=abs(end_pos.y-start_pos.y);
       double dZ=abs(end_pos.z-start_pos.z);
       double length=TMath::Sqrt(dX*dX+dY*dY+dZ*dZ);
       double cosL=dZ/length;
       //std::cout<<end_pos.z<<","<<start_pos.z<<std::endl;
	recoBacktrackCoslAr->Fill(cosL);
       }}
  
}
if (!goodInteraction && sr->mc.nu[biggestMatchIndex].id>1E9){  recoBacktrackCCRock->Fill(sr->mc.nu[biggestMatchIndex].iscc);
recoBacktrackPDGRock->Fill(sr->mc.nu[biggestMatchIndex].pdg);


}
if (!goodInteraction && sr->mc.nu[biggestMatchIndex].id<1E9){ recoBacktrackCCSec->Fill(sr->mc.nu[biggestMatchIndex].iscc);
recoBacktrackPDGSec->Fill(sr->mc.nu[biggestMatchIndex].pdg);


}
		track_mult->Fill(trackMult);

              if (goodInteraction){ track_multGood->Fill(trackMult);
        int cc=sr->mc.nu[biggestMatchIndex].iscc;
        int mode=sr->mc.nu[biggestMatchIndex].mode;
		track_multMode->Fill(mode);
	if (cc==1){ 
	   if (mode==1 || mode==1001) track_multQE->Fill(trackMult);
           if (mode==10)     track_multMEC->Fill(trackMult);
         if (mode==3)       track_multDIS->Fill(trackMult);
         if (mode==4) track_multRES->Fill(trackMult);
         if (mode==5) track_multCOH->Fill(trackMult);
	
           }
        else track_multNC->Fill(trackMult);
		}
              if (!goodInteraction){
        if (rock==1) track_multRock->Fill(trackMult); else track_multSec->Fill(trackMult);	
        bad_origin->Fill(rock);      track_multBad->Fill(trackMult);
         recoVertex2DBad->Fill(sr->common.ixn.dlp[nixn].vtx.x,sr->common.ixn.dlp[nixn].vtx.z);  trueVertex2DBad->Fill(trueVtxX,trueVtxZ);
          recoVertex2DBadYZ->Fill(sr->common.ixn.dlp[nixn].vtx.y,sr->common.ixn.dlp[nixn].vtx.z);
        } 
	}
    
    part_mult->Fill(partMult);      
 
	
	} //end for interaction
  }// end for spills
std::cout<<"Total: "<<true_mult->GetEntries()<<std::endl;
std::cout<<"All Interactions No Additional Selection: "<<double(goodIntNum)/true_mult->GetEntries()<<","<<float(goodIntNum)/float(goodIntNum+badIntNum)<<std::endl;
std::cout<<"In fiducial volume: "<<double(goodIntInFidVol)/true_mult->GetEntries()<<","<<float(goodIntInFidVol)/float(goodIntInFidVol+badIntInFidVol)<<std::endl;
std::cout<<"All events with MINERvA match: "<<track_multGood->GetEntries()<<"/"<<track_mult->GetEntries()<<","<<double(track_multGood->GetEntries())/true_mult->GetEntries()<<","<<double(track_multGood->GetEntries())/track_mult->GetEntries()<<std::endl;  
std::cout<<"All events with MINERvA match and <=2 tracks: "<<track_multGood->Integral(3,50)<<"/"<<track_mult->Integral(3,50)<<","<<double(track_multGood->Integral(3,50))/true_mult->Integral(3,50)<<","<<double(track_multGood->Integral(3,50))/track_mult->Integral(3,50)<<std::endl; 

//std::cout<<"Minerva Match Purity: "<<goodMINERvAMatch<<"/"<<totalMINERvAMatch<<std::endl;
//Create output file and write your histograms
  TFile *caf_out_file = new TFile(output_rootfile.c_str(), "recreate");
  matchHistCosl->Write();
  matchHistEl->Write();
  histCosl->Write();
  histEl->Write();
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
recoBacktrackCCAr->Write(); recoBacktrackCCRock->Write();
   recoBacktrackCCSec->Write();
recoBacktrackPDGSec->Write(); recoBacktrackPDGAr->Write();
recoBacktrackPDGRock->Write();

recoBacktrackElAr->Write(); recoBacktrackCoslAr->Write();


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
  recoVertex2DNoCuts->Write();
  longest_track->Write();
  showerCorrectness->Write();
  trackCorrectness->Write();
  selPionE->Write(); selKaonE->Write(); selProtonE->Write();

 truePionEWithRecoInt->Write();
 trueProtonEWithRecoInt->Write();
 trueKaonEWithRecoInt->Write();


  diffPip->Write(); diffPim->Write(); diffP->Write(); diffN->Write();
  recoVertex2DBadYZ->Write();  
  recoVertex2DBad->Write(); trueVertex2DBad->Write();

  auto protonKEEff=new TEfficiency(*selProtonE,*trueProtonEWithRecoInt);
  auto kaonKEEff=new TEfficiency(*selKaonE,*trueKaonEWithRecoInt);
  auto pionKEEff=new TEfficiency(*selPionE,*truePionEWithRecoInt);

  auto protonDirX=new TEfficiency(*selProtonDirX,*trueProtonWithRecoIntDirX);
  auto protonDirY=new TEfficiency(*selProtonDirY,*trueProtonWithRecoIntDirY);
  auto protonDirZ=new TEfficiency(*selProtonDirZ,*trueProtonWithRecoIntDirZ);

  auto pionDirX=new TEfficiency(*selPionDirX,*truePionWithRecoIntDirX);
  auto pionDirY=new TEfficiency(*selPionDirY,*truePionWithRecoIntDirY);
  auto pionDirZ=new TEfficiency(*selPionDirZ,*truePionWithRecoIntDirZ);

  auto kaonDirX=new TEfficiency(*selKaonDirX,*trueKaonWithRecoIntDirX);
  auto kaonDirY=new TEfficiency(*selKaonDirY,*trueKaonWithRecoIntDirY);
  auto kaonDirZ=new TEfficiency(*selKaonDirZ,*trueKaonWithRecoIntDirZ);

  protonDirX->Write(); protonDirY->Write(); protonDirZ->Write();
  pionDirX->Write(); pionDirY->Write(); pionDirZ->Write();
  kaonDirX->Write(); kaonDirY->Write(); kaonDirZ->Write();

  protonKEEff->Write();
  kaonKEEff->Write();
  pionKEEff->Write();


  caf_out_file->Close();
  
  return 1;  
}

int main(int argc, char** argv){

  if(argc!=3){
    std::cout << "\n USAGE: " << argv[0] << "input_caf_file_list output_root_file\n" << std::endl;
    return 1;
  }

  std::string input_file_list = argv[1];
  std::string output_rootfile = argv[2];

  caf_plotter(input_file_list, output_rootfile);

  return 0;
}
