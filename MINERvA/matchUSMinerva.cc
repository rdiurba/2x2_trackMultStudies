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

int matchDSMinerva(std::string input_file_list, std::string output_rootfile){
 
   int goodMINERvAMatch=0; int totalMINERvAMatch=0; int goodMINERvAMatchUS=0; int totalMINERvAMatchUS=0;
  //Give an input list
  std::ifstream caf_list(input_file_list.c_str());
  TH1D* dotProduct=new TH1D("dotProduct","dotProduct",30,0.9975,1);
  TH1D* dotProductGood=new TH1D("dotProductGood","dotProductGood",30,0.9975,1);
  TH1D* deltaX=new TH1D("deltaX","deltaX",50,-25,25);
  TH1D* deltaY=new TH1D("deltaY","deltaY",50,-25,25);
  TH1D* deltaYGood=new TH1D("deltaYGood","deltaYGood",50,-25,25);
  TH1D* deltaYBad=new TH1D("deltaYBad","deltaYBad",50,-25,25);
  TH1D* deltaXGood=new TH1D("deltaXGood","deltaXGood",50,-25,25);
  TH1D* deltaXBad=new TH1D("deltaXBad","deltaXBad",50,-25,25);

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
 	for(long unsigned nixn=0; nixn < sr->common.ixn.dlp.size(); nixn++){ 

    	for(long unsigned npart=0; npart < sr->common.ixn.dlp[nixn].part.dlp.size(); npart++){ //loop over particles
              //  if (!sr->common.ixn.dlp[nixn].part.dlp[npart].primary) continue;  Consider reco primaries only, you can toggle this if you want to
               int pdg=sr->common.ixn.dlp[nixn].part.dlp[npart].pdg;
int maxIxnNumber=-999; int maxPartNumber=-999; int maxTypeNumber=-999; int correctTrack=2;
		if ( (abs(pdg)==2212 || abs(pdg)==13 || abs(pdg)==211 || abs(pdg)==321)){  // Only look at tracks, only consider the best mach truth particle
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
			if (partTruthOverlap==maxPartTruthOverlap){  
			maxPartNumber=partNumber;
			maxTypeNumber=parType;
                        maxIxnNumber=interactionNumber;
			if (parType<3) backtracked=sr->mc.nu[interactionNumber].prim[partNumber].pdg;
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


	       auto tpcTrk_start=sr->common.ixn.dlp[nixn].part.dlp[npart].start;
	       auto tpcTrk_end=sr->common.ixn.dlp[nixn].part.dlp[npart].end;
	       if (std::isnan(tpcTrk_start.z)) continue;
		if (tpcTrk_end.z<tpcTrk_start.z){
		auto tpcTrk_start_temp=tpcTrk_end;
		tpcTrk_end=tpcTrk_start;
		tpcTrk_start=tpcTrk_start_temp;

		}
	       double dX=(tpcTrk_end.x-tpcTrk_start.x);
	       double dY=(tpcTrk_end.y-tpcTrk_start.y);
	       double dZ=(tpcTrk_end.z-tpcTrk_start.z);
	       double length=TMath::Sqrt(dX*dX+dY*dY+dZ*dZ);
		double dirX=dX/length; double dirY=dY/length; double dirZ=dZ/length;

                if (dirZ<0){ dirZ=-dirZ; dirX=-dirX; dirY=-dirY; auto temp=tpcTrk_start; tpcTrk_end=tpcTrk_start; tpcTrk_end=temp;}  // Flip the track so end is greater than start in Z

		if (std::isnan(tpcTrk_start.z)) length=-999;

		if (length>5){
		if ((tpcTrk_start.z<-58 || tpcTrk_end.z<-58) ){
                int punchthrough=-1;
		double dotProductDS=-999; double deltaExtrapYUS=-999; double deltaExtrapY=-999; double dotProductUS=-999; double deltaExtrapX=-999; double deltaExtrapXUS=-999;	
	 	int maxTypeMinerva=-999; int maxPartMinerva=-999; int maxIxnMinerva=-999;
                for(int i=0; i<sr->nd.minerva.ixn.size(); i++){    // Loop over MINERvA tracks

		for (int j=0; j<sr->nd.minerva.ixn[i].ntracks; j++){
		    double dir_z=sr->nd.minerva.ixn[i].tracks[j].dir.z;
		double end_z=sr->nd.minerva.ixn[i].tracks[j].end.z;
		double start_z=sr->nd.minerva.ixn[i].tracks[j].start.z;
		double end_x=sr->nd.minerva.ixn[i].tracks[j].end.x;
		double start_x=sr->nd.minerva.ixn[i].tracks[j].start.x;

		double end_y=sr->nd.minerva.ixn[i].tracks[j].end.y;
		double start_y=sr->nd.minerva.ixn[i].tracks[j].start.y;
                
		if (start_z<0 && end_z<0){  // Minerva track needs to be throughgoing
		int truthPart=sr->nd.minerva.ixn[i].tracks[j].truth[0].part;
		double dXMnv=(sr->nd.minerva.ixn[i].tracks[j].end.x-sr->nd.minerva.ixn[i].tracks[j].start.x);   // Find the direction
		double dYMnv=(sr->nd.minerva.ixn[i].tracks[j].end.y-sr->nd.minerva.ixn[i].tracks[j].start.y);
		double dZMnv=(sr->nd.minerva.ixn[i].tracks[j].end.z-sr->nd.minerva.ixn[i].tracks[j].start.z);
		double lengthMinerva=TMath::Sqrt(dXMnv*dXMnv+dYMnv*dYMnv+dZMnv*dZMnv);
	        if (lengthMinerva<10) continue;		
          	double dirXMinerva=dXMnv/lengthMinerva;
		double dirYMinerva=dYMnv/lengthMinerva;
		double dirZMinerva=dZMnv/lengthMinerva;
		double dotProduct=dirXMinerva*dirX+dirYMinerva*dirY+dirZ*dirZMinerva;  // Calculate dot product
		double extrapdZ=end_z-tpcTrk_end.z;
		double extrapY=dirY/dirZ*(extrapdZ)+tpcTrk_end.y-end_y;      // Calculate extrapolated distortions
		double extrapX=dirX/dirZ*(extrapdZ)+tpcTrk_end.x-end_x;
		double diffExtrap=TMath::Sqrt(TMath::Power(extrapY-start_y,2));


		if (dotProductUS<dotProduct && abs(extrapY)<10){ dotProductUS=dotProduct;  // Select the track that is most parallel
		deltaExtrapY=extrapY;
		deltaExtrapX=extrapX;
		maxPartMinerva=sr->nd.minerva.ixn[i].tracks[j].truth[0].part;   // Save info of best match
		maxTypeMinerva=sr->nd.minerva.ixn[i].tracks[j].truth[0].type;
		maxIxnMinerva=sr->nd.minerva.ixn[i].tracks[j].truth[0].ixn;
                
		}}
		
		}} // Minerva
	 if (dotProductUS>0.9975 ){    // Sanity check on dot product 
         dotProduct->Fill(dotProductUS);
         deltaX->Fill(deltaExtrapX); deltaY->Fill(deltaExtrapY);	
	if (maxIxnMinerva==maxIxnNumber && maxPartNumber==maxPartMinerva && maxTypeMinerva==maxTypeNumber){ goodMINERvAMatch++; deltaXGood->Fill(deltaExtrapX); deltaYGood->Fill(deltaExtrapY); dotProductGood->Fill(dotProductUS);}
  else{ deltaYBad->Fill(deltaExtrapY); deltaXBad->Fill(deltaExtrapX);}
	      totalMINERvAMatch++;       
	 } // done saving the information
		} // exiting particles
		} // primary of particle greater than 5
		   }} // particles


	} // interactions


  }// end for spills
 
std::cout<<"Minerva Match Purity US: "<<goodMINERvAMatch<<"/"<<totalMINERvAMatch<<std::endl;
//Create output file and write your histograms
  TFile *caf_out_file = new TFile(output_rootfile.c_str(), "recreate");  
  deltaX->Write(); deltaY->Write();
  deltaYGood->Write();
  deltaYBad->Write();
  deltaXGood->Write();
  deltaXBad->Write();
  dotProduct->Write(); 
  dotProductGood->Write();

  





  caf_out_file->Write();
  caf_out_file->Close();

    return 1;  
}

int main(int argc, char** argv){

  if(argc<3){
    std::cout << "\n USAGE: " << argv[0] << "input_caf_file_list output_root_file  \n" << std::endl;
    return 1;
  }

  std::string input_file_list = argv[1];
  std::string output_rootfile = argv[2];
  
  matchDSMinerva(input_file_list, output_rootfile);

  return 0;
}
