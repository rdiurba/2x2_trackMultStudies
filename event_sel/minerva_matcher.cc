#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
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

int caf_plotter(std::string input_file_list, std::string output_rootfile, bool mcOnly){

   int goodMINERvAMatch=0; int totalMINERvAMatch=0; int goodMINERvAMatchUS=0; int totalMINERvAMatchUS=0;
  //Give an input list
  std::ifstream caf_list(input_file_list.c_str());
  TH1D* dotProduct=new TH1D("dotProduct","dotProduct",30,0.9975,1);
  TH1D* dotProductGood=new TH1D("dotProductGood","dotProductGood",30,0.9975,1);
  TH1D* dotProductPlotUS=new TH1D("dotProductUS","dotProductUS",30,0.9975,1);
  TH1D* deltaX=new TH1D("deltaX","deltaX",60,-30,30);
  TH1D* deltaY=new TH1D("deltaY","deltaY",60,-30,30);
  TH1D* deltaYGood=new TH1D("deltaYGood","deltaYGood",50,-25,25);
  TH1D* deltaYBad=new TH1D("deltaYBad","deltaYBad",50,-25,25);
  TH1D* deltaXGood=new TH1D("deltaXGood","deltaXGood",50,-25,25);
  TH1D* deltaXBad=new TH1D("deltaXBad","deltaXBad",50,-25,25);

  TH1D* deltaYGoodUS=new TH1D("deltaYGoodUS","deltaYGoodUS",50,-25,25);
  TH1D* deltaYBadUS=new TH1D("deltaYBadUS","deltaYBadUS",50,-25,25);
  TH1D* deltaXGoodUS=new TH1D("deltaXGoodUS","deltaXGoodUS",50,-25,25);
  TH1D* deltaXBadUS=new TH1D("deltaXBadUS","deltaXBadUS",50,-25,25);

  TH1D* deltaXUS=new TH1D("deltaXUS","deltaXUS",60,-30,30);
  TH1D* deltaYUS=new TH1D("deltaYUS","deltaYUS",60,-30,30);
    TH1D* deltaXUSFront=new TH1D("deltaXUSFront","deltaXUSFront",60,-30,30);
  TH1D* deltaYUSFront=new TH1D("deltaYUSFront","deltaYUSFront",60,-30,30);

  TH1D* trackStartXUS= new TH1D("trackStartXUS","trackStartXUS",60,-60,60);
  TH1D* trackStartYUS= new TH1D("trackStartYUS","trackStartYUS",60,-60,60);
  TH2D* trackStartUS=new TH2D("trackStartUS","trackStartUS",60,-60,60,60,-60,60);



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

  for(long n = 0; n < Nentries; ++n){ //Loop over spills/triggers
//if(std::find(goodEvents.begin(), goodEvents.end(), n) != goodEvents.end()){}
//else continue;

	if(n%10000 == 0) std::cout << Form("Processing trigger %ld of %ld", n, Nentries) << std::endl;
	caf_chain->GetEntry(n); //Get spill from tree
       double mnvOffsetX=-10; double mnvOffsetY=5;
       if (mcOnly){ mnvOffsetX=0; mnvOffsetY=0;}
	
	for(long unsigned nixn = 0; nixn < sr->common.ixn.dlp.size(); nixn++){
          double biggestMatch=-999; int biggestMatchIndex=-999; double maxDotProductDS=-999; double maxDotProductUS=-999;   
            int maxEventPar=-999; int maxEventTyp=-9999; int maxEventIxn=-999;
           if (mcOnly){
           for (int ntruth=0; ntruth<sr->common.ixn.dlp[nixn].truth.size(); ntruth++){
          
          if (biggestMatch<sr->common.ixn.dlp[nixn].truthOverlap.at(ntruth)){
          
          biggestMatch=sr->common.ixn.dlp[nixn].truthOverlap.at(ntruth);
          biggestMatchIndex=sr->common.ixn.dlp[nixn].truth.at(ntruth);
           //std::cout<<"Biggest Match: "<<biggestMatch<<std::endl;


		}
		}
 		}
    	for(long unsigned npart=0; npart < sr->common.ixn.dlp[nixn].part.dlp.size(); npart++){ //loop over particles
                if (!sr->common.ixn.dlp[nixn].part.dlp[npart].primary) continue; 
               int pdg=sr->common.ixn.dlp[nixn].part.dlp[npart].pdg;
		int maxIxnNumber=-9999; int maxPartNumber=-999; int maxTypeNumber=-999; int correctTrack=2;
		if ( (abs(pdg)==2212 || abs(pdg)==13 || abs(pdg)==211 || abs(pdg)==321)){
			if (mcOnly){
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
                        maxIxnNumber=interactionNumber;
			if (parType<3){ backtracked=sr->mc.nu[interactionNumber].prim[partNumber].pdg;
                        //std::cout<<sr->common.ixn.dlp[nixn].part.dlp[npart].start.z<<","<<sr->mc.nu[maxIxnNumber].prim[maxPartNumber].start_pos.z<<std::endl;

			}
			else backtracked=sr->mc.nu[interactionNumber].sec[partNumber].pdg;
			correctTrack=0;
			//std::cout<<"Backtracking for "<<parType<<" : "<<pdg<<","<<backtracked<<std::endl;
			if (parType>2) correctTrack=3;
			else if (backtracked==11 || backtracked==22 || backtracked==111 || backtracked==2112){
			correctTrack=1;
			}



			}
			if (partTruthOverlap<0.5) correctTrack=2;
			else backtracked=-1; 
			}
			}
                //if (maxPartTruthOverlap<0.5) continue;

	       auto start_pos=sr->common.ixn.dlp[nixn].part.dlp[npart].start;
	       auto end_pos=sr->common.ixn.dlp[nixn].part.dlp[npart].end;

		double diffVertexdZ=abs(start_pos.z-sr->common.ixn.dlp[nixn].vtx.z);
		double diffVertexdX=abs(start_pos.x-sr->common.ixn.dlp[nixn].vtx.x);
		double diffVertexdY=abs(start_pos.y-sr->common.ixn.dlp[nixn].vtx.y);
		double diffVertex=TMath::Sqrt(diffVertexdZ*diffVertexdZ+diffVertexdY*diffVertexdY+diffVertexdX*diffVertexdX);
		if (diffVertex>5) continue;
	       if (std::isnan(start_pos.z)) continue;
	       double dX=(end_pos.x-start_pos.x);
	       double dY=(end_pos.y-start_pos.y);
	       double dZ=(end_pos.z-start_pos.z);
	       double length=TMath::Sqrt(dX*dX+dY*dY+dZ*dZ);
		double dirX=dX/length; double dirY=dY/length; double dirZ=dZ/length;

                if (dirZ<0){ dirZ=-dirZ; dirX=-dirX; dirY=-dirY;}

		if (std::isnan(start_pos.z)) length=-999;

		if (sr->common.ixn.dlp[nixn].part.dlp[npart].primary==true && length>1){
  
		int maxPartMinerva=-999; int maxTypeMinerva=-999;    int maxIxnMinerva=-999;   int maxIxnMinervaUS=-999;
  
		int maxPartMinervaUS=-999; int maxTypeMinervaUS=-999;
                

if ( (start_pos.z<-60 && end_pos.z>60) || (start_pos.z>60 && end_pos.z<-60)){
		trackStartXUS->Fill(start_pos.x); trackStartYUS->Fill(start_pos.y);
          trackStartUS->Fill(start_pos.x,start_pos.y);      }
		if ((abs(start_pos.z)>58 || abs(end_pos.z)>58) ){
		double dotProductDS=-999; double deltaExtrapYUS=-999; double deltaExtrapY=-999; double dotProductUS=-999; double deltaExtrapX=-999; double deltaExtrapXUS=-999;	
double deltaExtrapXUSFront=-999; double deltaExtrapYUSFront=-999;
	 	for(int i=0; i<sr->nd.minerva.ixn.size(); i++){

		for (int j=0; j<sr->nd.minerva.ixn[i].ntracks; j++){
		    double dir_z=sr->nd.minerva.ixn[i].tracks[j].dir.z;
		double end_z=sr->nd.minerva.ixn[i].tracks[j].end.z;
		double start_z=sr->nd.minerva.ixn[i].tracks[j].start.z;
		double end_x=sr->nd.minerva.ixn[i].tracks[j].end.x;
		double start_x=sr->nd.minerva.ixn[i].tracks[j].start.x;

		double end_y=sr->nd.minerva.ixn[i].tracks[j].end.y;
		double start_y=sr->nd.minerva.ixn[i].tracks[j].start.y;
                



		if (start_z>0 && end_z>0 && ((start_pos.z)>60 || (end_pos.z)>60) ){


          if (/*abs(abs(sr->common.ixn.dlp[nixn].vtx.x)-33)<1 ||*/ abs(sr->common.ixn.dlp[nixn].vtx.x)>59 || abs(sr->common.ixn.dlp[nixn].vtx.x)<5 || abs(sr->common.ixn.dlp[nixn].vtx.y)>57 || abs(sr->common.ixn.dlp[nixn].vtx.z)<5 || abs(sr->common.ixn.dlp[nixn].vtx.z)>59.5)  continue;


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


		if (dotProductDS<dotProduct && abs(extrapY-mnvOffsetY)<15 && abs(extrapX-mnvOffsetX)<15  &&  abs(TMath::ATan(dirXMinerva/dirZMinerva)-TMath::ATan(dirX/dirZ))<0.06 && abs(TMath::ATan(dirYMinerva/dirZMinerva)-TMath::ATan(dirY/dirZ))<0.06){ dotProductDS=dotProduct;
		deltaExtrapY=extrapY;
		deltaExtrapX=extrapX;
		if (mcOnly){
		maxPartMinerva=sr->nd.minerva.ixn[i].tracks[j].truth[0].part;
		maxTypeMinerva=sr->nd.minerva.ixn[i].tracks[j].truth[0].type;
		maxIxnMinerva=sr->nd.minerva.ixn[i].tracks[j].truth[0].ixn;
		}

		}}

		if (start_z<0 && end_z>0 && ( (start_pos.z<-60 && end_pos.z>60) || (start_pos.z>60 && end_pos.z<-60))){
	       if (start_z>-1000 && start_z<-900)	std::cout<<start_z<<std::endl;
	int truthPart=sr->nd.minerva.ixn[i].tracks[j].truth[0].part;
		double dXMnv=(sr->nd.minerva.ixn[i].tracks[j].end.x-sr->nd.minerva.ixn[i].tracks[j].start.x);
		double dYMnv=(sr->nd.minerva.ixn[i].tracks[j].end.y-sr->nd.minerva.ixn[i].tracks[j].start.y);
		double dZMnv=(sr->nd.minerva.ixn[i].tracks[j].end.z-sr->nd.minerva.ixn[i].tracks[j].start.z);
		double lengthMinerva=TMath::Sqrt(dXMnv*dXMnv+dYMnv*dYMnv+dZMnv*dZMnv);
		double dirXMinerva=dXMnv/lengthMinerva;
		double dirYMinerva=dYMnv/lengthMinerva;
		double dirZMinerva=dZMnv/lengthMinerva;
		double dotProductTemp=dirXMinerva*dirX+dirYMinerva*dirY+dirZ*dirZMinerva;

		double extrapdZUS=end_z-end_pos.z;
		double extrapYUS=dirY/dirZ*(extrapdZUS)+end_pos.y-end_y;
		double extrapXUS=dirX/dirZ*(extrapdZUS)+end_pos.x-end_x;
		double extrapdZUSFront=start_z-start_pos.z;
	        double extrapYUSFront=dirY/dirZ*(extrapdZUSFront)+start_pos.y-start_y;
                double extrapXUSFront=dirX/dirZ*(extrapdZUSFront)+start_pos.x-start_x;	
                //double diffExtrap=TMath::Sqrt(TMath::Power(extrapY-end_y,2));
		if (dotProductUS<dotProductTemp  && abs(extrapYUS-mnvOffsetY)<15 && abs(extrapXUS-mnvOffsetX)<15 && abs(TMath::ATan(dirXMinerva/dirZMinerva)-TMath::ATan(dirX/dirZ))<0.06 && abs(TMath::ATan(dirYMinerva/dirZMinerva)-TMath::ATan(dirY/dirZ))<0.06){ dotProductUS=dotProductTemp;
		deltaExtrapYUS=extrapYUS;
                deltaExtrapXUS=extrapXUS;
                deltaExtrapXUSFront=extrapXUSFront;
                deltaExtrapYUSFront=extrapYUSFront;
		if (mcOnly){
                maxPartMinervaUS=sr->nd.minerva.ixn[i].tracks[j].truth[0].part;
		maxTypeMinervaUS=sr->nd.minerva.ixn[i].tracks[j].truth[0].type;	
		maxIxnMinervaUS=sr->nd.minerva.ixn[i].tracks[j].truth[0].ixn;
                }}


//if (maxPartMinervaUS==maxPartNumber) std::cout<<extrapdZ<<","<<dirY/dirZ*(extrapdZ)<<","<<double(dirX/dirZ*(extrapdZ))<<","<<end_y<<","<<end_x<<std::endl;
		}

		
		}} // Minerva
	if (dotProductDS>maxDotProductDS) maxDotProductDS=dotProductDS;
	 if (dotProductDS>0.99){   
         dotProduct->Fill(dotProductDS);
         deltaX->Fill(deltaExtrapX); deltaY->Fill(deltaExtrapY);	
	if (mcOnly && maxIxnMinerva==maxIxnNumber && maxPartNumber==maxPartMinerva && maxTypeMinerva==maxTypeNumber){ goodMINERvAMatch++; deltaXGood->Fill(deltaExtrapX); deltaYGood->Fill(deltaExtrapY); dotProductGood->Fill(dotProductDS);
          if (maxTypeMinerva==1){
               auto start_posMnv=sr->mc.nu[maxIxnMinerva].prim[maxPartMinerva].start_pos;
       auto end_posMnv=sr->mc.nu[maxIxnMinerva].prim[maxPartMinerva].end_pos;
       if (!std::isnan(start_posMnv.z)){
       double dX=abs(end_posMnv.x-start_posMnv.x);
       double dY=abs(end_posMnv.y-start_posMnv.y);
       double dZ=abs(end_posMnv.z-start_posMnv.z);
       auto p=sr->mc.nu[maxIxnMinerva].prim[maxPartMinerva].p;
       double length=TMath::Sqrt(p.px*p.px+p.py*p.py+p.pz*p.pz);
       double cosL=p.pz/length;
       
      
       double lengthDist=TMath::Sqrt(dX*dX+dY*dY+dZ*dZ);
       double cosLDist=dZ/lengthDist;




       //if (abs(start_pos.z-start_posMnv.z)<5)   minervaMatchCos->Fill(cosL);
       //if (cosLDist<0.9) std::cout<<start_pos.z<<","<<end_pos.z<<","<<start_posMnv.z<<","<<end_posMnv.z<<","<<cosL<<","<<cosLDist<<","<<sr->mc.nu[maxIxnMinerva].prim[maxPartMinerva].pdg<<std::endl;

         } 
	}}
  else{ deltaYBad->Fill(deltaExtrapY); deltaXBad->Fill(deltaExtrapX);}
	      totalMINERvAMatch++;       
	 }
	if (dotProductUS>maxDotProductUS) maxDotProductUS=dotProductUS;
	 if (dotProductUS>0.99){   
	  dotProductPlotUS->Fill(dotProductUS);
	 deltaXUSFront->Fill(deltaExtrapXUSFront); deltaYUSFront->Fill(deltaExtrapYUSFront);	
         deltaXUS->Fill(deltaExtrapXUS); deltaYUS->Fill(deltaExtrapYUS);
         /*trackStartXUS->Fill(start_pos.x); trackStartYUS->Fill(start_pos.y);
          trackStartUS->Fill(start_pos.x,start_pos.y);
		*/

 if (mcOnly && maxIxnMinervaUS==maxIxnNumber && maxPartNumber==maxPartMinervaUS && maxTypeMinervaUS==maxTypeNumber){ deltaYGoodUS->Fill(deltaExtrapYUS); deltaXGoodUS->Fill(deltaExtrapXUS); goodMINERvAMatchUS++;} 
              else{ deltaYBadUS->Fill(deltaExtrapYUS); deltaXBadUS->Fill(deltaExtrapXUS); }
       totalMINERvAMatchUS++;       
	 }
		} // exiting particles
		} // primary of particle greater than 5
		   }} // particles


	} // interactions


  }// end for spills
  
std::cout<<"Minerva Match Purity US: "<<goodMINERvAMatchUS<<"/"<<totalMINERvAMatchUS<<std::endl;
std::cout<<"Minerva Match Purity DS: "<<goodMINERvAMatch<<"/"<<totalMINERvAMatch<<std::endl;
//Create output file and write your histograms
  TFile *caf_out_file = new TFile(output_rootfile.c_str(), "recreate");  
  deltaX->Write(); deltaY->Write();
  deltaYGood->Write();
  deltaYBad->Write();
  deltaXGood->Write();
  deltaXBad->Write();
  dotProduct->Write(); 
  dotProductGood->Write();
  dotProductPlotUS->Write();
  deltaYGoodUS->Write();
  deltaYBadUS->Write();
  deltaXGoodUS->Write();
  deltaXBadUS->Write();
  deltaXUS->Write();
  deltaYUS->Write(); 
  trackStartXUS->Write();
  trackStartYUS->Write();
  trackStartUS->Write();
  deltaXUSFront->Write();
  deltaYUSFront->Write();



  caf_out_file->Write();
  caf_out_file->Close();

    return 1;  
}

int main(int argc, char** argv){

  if(argc!=4){
    std::cout << "\n USAGE: " << argv[0] << "input_caf_file_list output_root_file\n" << std::endl;
    return 1;
  }

  std::string input_file_list = argv[1];
  std::string output_rootfile = argv[2];
  std::string mcOnlyString=argv[3];
  bool mcOnly=true;
  if (mcOnlyString=="0") mcOnly=false;
  std::cout<<mcOnly<<","<<argv[3]<<std::endl;

  caf_plotter(input_file_list, output_rootfile,mcOnly);

  return 0;
}
