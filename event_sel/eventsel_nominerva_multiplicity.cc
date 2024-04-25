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
double minTrkLength=1;
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
  TH1D *track_multBad=new TH1D("track_multBad","track_multBad",20,0,20);
  TH1D *track_multGood=new TH1D("track_multGood","track_multGood",20,0,20);


  TH1D *track_multFirst= new TH1D("track_multFirst", "track_multFirst", 20, 0, 20);
  TH1D *track_multBadFirst=new TH1D("track_multBadFirst","track_multBadFirst",20,0,20);
  TH1D *track_multGoodFirst=new TH1D("track_multGoodFirst","track_multGoodFirst",20,0,20);
  TH1D *bad_originFirst=new TH1D("rockFirst","rockFirst",2,0,2);


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
 TH1D* recoBacktrackCoslAr=new TH1D("recoBacktrackCoslAr","recoBackTrackCoslAr",30,0.9,1);

 TH1D* selPionE=new TH1D("selPionE","selPionE",20,0,0.1);
 TH1D* selProtonE=new TH1D("selProtonE","selProtonE",20,0,0.1);
 TH1D* selKaonE=new TH1D("selKaonE","selKaonE",20,0,0.1);

 TH1D* truePionEWithRecoInt=new TH1D("truePionEWithRecoInt","truePionEWithRecoInt",20,0,0.1);
 TH1D* trueProtonEWithRecoInt=new TH1D("trueProtonEWithRecoInt","trueProtonEWithRecoInt",20,0,0.1);
 TH1D* trueKaonEWithRecoInt=new TH1D("trueKaonEWithRecoInt","trueKaonEWithRecoInt",20,0,0.1);



  TH1D* recoCosL=new TH1D("recoCosL","recoCosL",50,-1,1);
  TH1D *true_mult=new TH1D("true_mult","true_mult",20,0,20);
  TH1D *true_multFirst=new TH1D("true_multFirst","true_multFirst",20,0,20);
  TH1D *true_multTrkOnly=new TH1D("true_multTrkOnly","true_multTrkOnly",20,0,20);
  TH1D *true_multGENIE=new TH1D("true_multGENIE","true_multGENIE",20,0,20);

  TH1D *matchTrue_mult=new TH1D("matchTrue_mult","matchTrue_mult",20,0,20);
  TH1D *matchTrue_multTrkOnly=new TH1D("matchTrue_multTrkOnly","matchTrue_multTrkOnly",20,0,20);
  TH1D *matchTrue_multGENIE=new TH1D("matchTrue_multGENIE","matchTrue_multGENIE",20,0,20);
  TH1D *matchHistEl=new TH1D("matchHistEl","matchHistEl",50,0,20);
  TH1D *matchHistCosl=new TH1D("matchHistCosl","matchHistCosl",50,-1,1);

  TH1D *histEl=new TH1D("histEl","histEl",50,0,20);
  TH1D *histCosl=new TH1D("histCosl","histCosl",50,-1,1);

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
       if (length>minTrkLength){ truePartTrkOnly++;
 

       }
       }
       truePart++;
       }
       true_mult->Fill(truePart); if (truePartTrkOnly>0) true_multTrkOnly->Fill(truePartTrkOnly);    true_multGENIE->Fill(truePartNoG4); if (trueVtxZ<0) true_multFirst->Fill(truePart);
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
        double dirZExiting=-999;   
       for (int ntruth=0; ntruth<sr->common.ixn.dlp[nixn].truth.size(); ntruth++){
          if (biggestMatch<sr->common.ixn.dlp[nixn].truthOverlap.at(ntruth)){
          biggestMatch=sr->common.ixn.dlp[nixn].truthOverlap.at(ntruth);
          biggestMatchIndex=sr->common.ixn.dlp[nixn].truth.at(ntruth);}}
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
      if (sr->mc.nu[biggestMatchIndex].prim.size()<500 && sr->mc.nu[biggestMatchIndex].targetPDG==1000180400){ 
       for (int primaries=0; primaries<sr->mc.nu[biggestMatchIndex].prim.size(); primaries++){ 
      
      auto start_pos=sr->mc.nu[biggestMatchIndex].prim[primaries].start_pos;
      if (std::isnan(start_pos.z)) continue;
      //std::cout<<"Got a position"<<std::endl; 
      if(abs(sr->mc.nu[biggestMatchIndex].prim[primaries].pdg)==2212) trueProtonEWithRecoInt->Fill(sr->mc.nu[biggestMatchIndex].prim[primaries].p.E-0.938272);
          if(abs(sr->mc.nu[biggestMatchIndex].prim[primaries].pdg)==211) truePionEWithRecoInt->Fill(sr->mc.nu[biggestMatchIndex].prim[primaries].p.E-0.13957);
      if(abs(sr->mc.nu[biggestMatchIndex].prim[primaries].pdg)==321) trueKaonEWithRecoInt->Fill(sr->mc.nu[biggestMatchIndex].prim[primaries].p.E-0.493677);
      
      }}
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
    int trackMult=0;  int trackMultExit=0; 
    		for(long unsigned npart=0; npart < sr->common.ixn.dlp[nixn].part.dlp.size(); npart++){ //loop over particles
 	//std::cout<<"Containment variable: "<<sr->common.ixn.dlp[nixn].part.dlp[npart].contained<<std::endl;
                if (!sr->common.ixn.dlp[nixn].part.dlp[npart].primary) continue; 
               int pdg=sr->common.ixn.dlp[nixn].part.dlp[npart].pdg;
		if ( (abs(pdg)==2212 || abs(pdg)==13 || abs(pdg)==211 || abs(pdg)==321)){

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
        if (sr->common.ixn.dlp[nixn].part.dlp[npart].primary==true && length>minTrkLength){   
		partMult++;
                trackMult++;
           if    ((start_pos.z)>58 || (end_pos.z)>58){ trackMultExit++; 

               if(dirZ>dirZExiting) dirZExiting=dirZ;


		}
		int maxPartMinerva=-999; int maxTypeMinerva=-999;  
                  int maxIxnMinerva=-999;
		int maxPartMinervaUS=-999; int maxTypeMinervaUS=-999;      
		
		} // primary of particle greater than 5
		   }} // particles

         longest_track->Fill(longestTrk);
       if (trackMult>0 && trackMultExit>0){  

recoCosL->Fill(dirZExiting);

if (goodInteraction){ recoBacktrackCCAr->Fill(sr->mc.nu[biggestMatchIndex].iscc);
recoBacktrackPDGAr->Fill(sr->mc.nu[biggestMatchIndex].pdg);
       if (sr->mc.nu[biggestMatchIndex].iscc && abs(sr->mc.nu[biggestMatchIndex].pdg)==14){
       for (int primaries=0; primaries<sr->mc.nu[biggestMatchIndex].prim.size(); primaries++){
       if (abs(sr->mc.nu[biggestMatchIndex].prim[primaries].pdg)!=13) continue;
       double Elep=sr->mc.nu[biggestMatchIndex].prim[primaries].p.E;
            auto start_pos=sr->mc.nu[biggestMatchIndex].prim[primaries].start_pos;
       auto end_pos=sr->mc.nu[biggestMatchIndex].prim[primaries].end_pos;
       if (std::isnan(start_pos.z)) continue;
       double dX=(end_pos.x-start_pos.x);
       double dY=(end_pos.y-start_pos.y);
       double dZ=(end_pos.z-start_pos.z);
       double length=TMath::Sqrt(dX*dX+dY*dY+dZ*dZ);
       double cosL=dZ/length;
       recoBacktrackElAr->Fill(Elep);
       recoBacktrackCoslAr->Fill(cosL);
       }

       }




}
if (!goodInteraction && sr->mc.nu[biggestMatchIndex].id>1E9){  recoBacktrackCCRock->Fill(sr->mc.nu[biggestMatchIndex].iscc);
recoBacktrackPDGRock->Fill(sr->mc.nu[biggestMatchIndex].pdg);


}
if (!goodInteraction && sr->mc.nu[biggestMatchIndex].id<1E9){ recoBacktrackCCSec->Fill(sr->mc.nu[biggestMatchIndex].iscc);
recoBacktrackPDGSec->Fill(sr->mc.nu[biggestMatchIndex].pdg);


}
		track_mult->Fill(trackMult); if (sr->common.ixn.dlp[nixn].vtx.z<0) track_multFirst->Fill(trackMult);
              if (goodInteraction){ track_multGood->Fill(trackMult);	 if (sr->common.ixn.dlp[nixn].vtx.z<0) track_multGoodFirst->Fill(trackMult);			}
              if (!goodInteraction){ bad_origin->Fill(rock);      track_multBad->Fill(trackMult);  if (sr->common.ixn.dlp[nixn].vtx.z<0){ bad_originFirst->Fill(rock); track_multBadFirst->Fill(trackMult);}
         recoVertex2DBad->Fill(sr->common.ixn.dlp[nixn].vtx.x,sr->common.ixn.dlp[nixn].vtx.z);  trueVertex2DBad->Fill(trueVtxX,trueVtxZ);
          recoVertex2DBadYZ->Fill(sr->common.ixn.dlp[nixn].vtx.y,sr->common.ixn.dlp[nixn].vtx.z);
        } 
	}
    
    part_mult->Fill(partMult);      
 
	
	} //end for interactions

 	for(long unsigned ntrue=0; ntrue<sr->mc.nu.size(); ntrue++){ 

        int matchTruePartTrkOnly=0; int matchTruePartNoG4=0;
       auto vertex=sr->mc.nu[ntrue].vtx;
       trueVtxX=vertex.x; trueVtxY=vertex.y; trueVtxZ=vertex.z;
       auto truePrimary=sr->mc.nu[ntrue].prim;
       int truePart=0; int truePartTrkOnly=0; int truePartNoG4=0;
       if (sr->mc.nu[ntrue].targetPDG!=1000180400) continue;
       if(std::find(trueInteractionIndex.begin(), trueInteractionIndex.end(), ntrue) == trueInteractionIndex.end()) continue;
	double missE=sr->mc.nu[ntrue].E;
       for(long unsigned primaries=0; primaries<truePrimary.size(); primaries++){
       
       int pdg=truePrimary[primaries].pdg;
       double E=truePrimary[primaries].p.E;
       double px=truePrimary[primaries].p.px; double py=truePrimary[primaries].p.py; double pz=truePrimary[primaries].p.pz;
       double totP=TMath::Sqrt(px*px+py*py+pz*pz);
       //recoTrueCC->Fill(sr->mc.nu[ntrue].iscc
       if ((abs(pdg)==13 || abs(pdg)==2212 || abs(pdg)==211 || abs(pdg)==321)){
        matchTruePartNoG4++;
       
       auto start_pos=sr->mc.nu[ntrue].prim[primaries].start_pos;
       auto end_pos=sr->mc.nu[ntrue].prim[primaries].end_pos;
       if (std::isnan(start_pos.z)) continue;
       double dX=abs(end_pos.x-start_pos.x);
       double dY=abs(end_pos.y-start_pos.y);
       double dZ=abs(end_pos.z-start_pos.z);
       double length=TMath::Sqrt(dX*dX+dY*dY+dZ*dZ);
       if (sr->mc.nu[ntrue].iscc && abs(sr->mc.nu[ntrue].pdg)==14 && abs(pdg)==13){
     
       double Elep=sr->mc.nu[ntrue].prim[primaries].p.E;
       double cosL=dZ/length;
       matchHistEl->Fill(Elep);
       matchHistCosl->Fill(cosL);


       }
       if (length>0){ matchTruePartTrkOnly++;
 

       }
       
	}
	}

	matchTrue_multTrkOnly->Fill(matchTruePartTrkOnly);
        matchTrue_multGENIE->Fill(matchTruePartNoG4);
        


        }
   


  }// end for spills
std::cout<<"Total: "<<true_mult->GetEntries()<<","<<true_multFirst->GetEntries()<<std::endl;
std::cout<<goodIntNum<<"/"<<goodIntNum+badIntNum<<","<<double(goodIntNum)/true_mult->GetEntries()<<","<<double(goodIntNum)/double(goodIntNum+badIntNum)<<std::endl;
std::cout<<goodIntInFidVol<<","<<goodIntInFidVol+badIntInFidVol<<","<<double(goodIntInFidVol)/true_mult->GetEntries()<<","<<double(goodIntInFidVol)/double(goodIntInFidVol+badIntInFidVol)<<std::endl;
std::cout<<"All modules: "<<track_multGood->GetEntries()<<"/"<<track_mult->GetEntries()<<","<<double(track_multGood->GetEntries())/true_mult->GetEntries()<<","<<double(track_multGood->GetEntries())/track_mult->GetEntries()<<std::endl;  
std::cout<<"All modules more than 2 tracks: "<<track_multGood->Integral(3,50)<<"/"<<track_mult->Integral(3,50)<<","<<double(track_multGood->Integral(3,50))/true_mult->Integral(3,50)<<","<<double(track_multGood->Integral(3,50))/track_mult->Integral(3,50)<<std::endl;  


std::cout<<"First modules: "<<track_multGoodFirst->GetEntries()<<"/"<<track_multFirst->GetEntries()<<","<<double(track_multGoodFirst->GetEntries())/true_mult->GetEntries()<<","<<double(track_multGoodFirst->GetEntries())/track_multFirst->GetEntries()<<std::endl;  
std::cout<<"First modules more than 2 tracks: "<<track_multGoodFirst->Integral(3,50)<<"/"<<track_multFirst->Integral(3,50)<<","<<double(track_multGoodFirst->Integral(3,50))/true_mult->Integral(0,50)<<","<<double(track_multGoodFirst->Integral(3,50))/track_multFirst->Integral(3,50)<<std::endl; 
std::cout<<"Minerva Match Purity: "<<goodMINERvAMatch<<"/"<<totalMINERvAMatch<<std::endl;
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

  track_multFirst->Write();
  track_multGoodFirst->Write();
  track_multBadFirst->Write();
  recoCosL->Write();
  bad_originFirst->Write();


  diffPip->Write(); diffPim->Write(); diffP->Write(); diffN->Write();
  recoVertex2DBadYZ->Write();  
  recoVertex2DBad->Write(); trueVertex2DBad->Write();

  auto protonKEEff=new TEfficiency(*selProtonE,*trueProtonEWithRecoInt);
  auto kaonKEEff=new TEfficiency(*selKaonE,*trueKaonEWithRecoInt);
  auto pionKEEff=new TEfficiency(*selPionE,*truePionEWithRecoInt);
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
