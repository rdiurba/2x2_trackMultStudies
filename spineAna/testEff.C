#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <iomanip>
#include <TH2D.h>
#include <TMatrixD.h>



double binomialProb(int N, int k, double delta)
{
    if (k < 0 || k > N) return 0.0;

    double logC =
        std::lgamma(N + 1)
      - std::lgamma(k + 1)
      - std::lgamma(N - k + 1);

    return std::exp(
        logC
      + k * std::log(1.0 - delta)
      + (N - k) * std::log(delta)
    );
}
TH2D* convertCovCountsToXS(
    const TMatrixD& covCounts,
    double flux,                     // neutrinos per m^2 
    const std::vector<double>& binWidths)  // width of each multiplicity bin
{
    float nTargets=2.824206e+28;
    
    int n = covCounts.GetNrows();
    TMatrixD covXS(n, n);

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            covXS(i,j) = covCounts(i,j) / 
                         (flux * flux * nTargets * nTargets * binWidths[i] * binWidths[j]);
        }
    }


    TH2D* h = new TH2D(
        "covMatXS",
        "covMatXS",
        n, 0, n,
        n, 0, n
    );

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            // ROOT bins start at 1
            h->SetBinContent(i+1, j+1, covXS(i,j));
        }
    }

    return h;

}
TMatrixD computeCovarianceBinned(TH2D* hTrueLS, double deltaLong, double deltaShort,
                                 int binsMult, Double_t* edgesMult)
{
    TMatrixD cov(binsMult, binsMult);
    cov.Zero();

    int nBinsX = hTrueLS->GetNbinsX();
    int nBinsY = hTrueLS->GetNbinsY();

    // Loop over true (N_L, N_S)
    for (int ix = 0; ix <= nBinsX; ++ix)
    {
        for (int iy = 0; iy <= nBinsY; ++iy)
        {
            double nEvents = hTrueLS->GetBinContent(ix, iy);
            if (nEvents == 0.0) continue;

            int N_S = hTrueLS->GetXaxis()->GetBinCenter(ix);
            int N_L = hTrueLS->GetYaxis()->GetBinCenter(iy);


            // Compute raw probabilities for each nObs
            std::vector<double> pObs(N_L + N_S + 1, 0.0);
            for (int nL = 0; nL <= N_L; ++nL)
            {
                for (int nS = 0; nS <= N_S; ++nS)
                {

                    
             double effLong = ((nL == 1 && nS==0) || (N_L==1 && N_S==0)) ? 1E-20 : deltaLong;
            double effShort = ((nL == 0 && nS==1 )|| (N_L==0 && N_S==1)) ? 1E-20 : deltaShort;
                    double pS = binomialProb(N_S, nS, effShort);
                    double pL = binomialProb(N_L, nL, effLong);
                    int nObs = nL + nS;
                    pObs[nObs] += pL * pS;
                }
            }

            // Normalize probabilities
            double sumProb = 0.0;
            for (double p : pObs) sumProb += p;
            for (double &p : pObs) p /= sumProb;

            // Map nObs to multiplicity bins
            std::vector<double> pBin(binsMult, 0.0);
            for (size_t nObs = 0; nObs < pObs.size(); ++nObs)
            {
                for (int b = 0; b < binsMult; ++b)
                {
                    if (nObs >= edgesMult[b] && nObs < edgesMult[b+1])
                    {
                        pBin[b] += pObs[nObs];
                        break;
                    }
                }
            }

            // Fill covariance matrix in multiplicity bins
            for (int i = 0; i < binsMult; ++i)
            {
                for (int j = 0; j < binsMult; ++j)
                {
                    double term = (i == j ? pBin[i] : 0.0) - pBin[i]*pBin[j];
                    cov(i,j) += nEvents * term;
                }
            }
        }
    }

    return cov;
}

TH1D* migrateMultiplicity(
    TH2D* hTrueLS,
    double deltaLong,
    double deltaShort
)
{




int maxMult = 10;
    // Output histogram: observed total multiplicity
    TH1D* hObs = new TH1D(
        "hObservedMult",
        "Observed multiplicity",
        20,0,20);


    int nBinsX = hTrueLS->GetNbinsX();
    int nBinsY = hTrueLS->GetNbinsY();

    for (int ix = 0; ix <= nBinsX; ++ix){
        for (int iy = 0; iy <= nBinsY; ++iy){
            double nEvents = hTrueLS->GetBinContent(ix, iy);
            if (nEvents == 0.0) continue;

            int N_L = hTrueLS->GetYaxis()->GetBinCenter(iy);
            int N_S = hTrueLS->GetXaxis()->GetBinCenter(ix);

                double totalProb = 0.0;
            for (int nL = 0; nL <= N_L; ++nL)
            {

                for (int nS = 0; nS <= N_S; ++nS)
                {
             double effLong = ((nL == 1 && nS==0) || (N_L==1 && N_S==0)) ? 1E-20 : deltaLong;
            double effShort = ((nL == 0 && nS==1 )|| (N_L==0 && N_S==1)) ? 1E-20 : deltaShort;
                    double pS = binomialProb(N_S, nS, effShort);
                    double pL = binomialProb(N_L, nL, effLong);

                    totalProb += pL * pS;
                }
            }





            // Loop over observed long/short tracks
            for (int nS = 0; nS <= N_S; ++nS)
                {

            
            
            for (int nL = 0; nL <= N_L; ++nL)
            {
             double effLong = ((nL == 1 && nS==0) || (N_L==1 && N_S==0)) ? 1E-20 : deltaLong;
            double effShort = ((nL == 0 && nS==1 )|| (N_L==0 && N_S==1)) ? 1E-20 : deltaShort;
                    double pS = binomialProb(N_S, nS, effShort);
                    double pL = binomialProb(N_L, nL, effLong);
  
                    int nObs = nL + nS;
                    double weight = nEvents * pL * pS/totalProb;
                    hObs->Fill(nObs, weight);
                }
            }
        }
    }

    return hObs;
}

int testEff() {
double deltaLong  = 0.03;
double deltaShort = 0.15;



TFile fPurity("testSPINEPart1.root");
TH2D* true_multResp=(TH2D*)fPurity.Get("true_multResp");
TH1D* totalPOT=(TH1D*)fPurity.Get("totalPOT");
totalPOT->SetName("totalPOTMC");
TH1D* true_mult=(TH1D*)fPurity.Get("true_multTrkOnly");
double flux=totalPOT->Integral()*stod(true_mult->GetTitle());
TH1D* hObsInt =
    migrateMultiplicity(true_multResp, deltaLong, deltaShort);


std::cout<<true_mult->GetEntries()<<","<<true_multResp->GetEntries()<<std::endl;
    int binsMult = 7;
Double_t edgesMult[8] = {1, 2, 3, 4, 5, 6, 8, 10};
TH1D* hObsFinal = (TH1D*)hObsInt->Rebin(
    binsMult,
    "hObsFinal",
    edgesMult
);
    
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Bin   Nominal   Shifted   RelChange(%)\n";
    std::cout << "--------------------------------------\n";

    for (int i = 0; i < 7; ++i)
    {
        double nominal = true_mult->GetBinContent(i+1);
        double shifted = hObsFinal->GetBinContent(i+1);

        double relChange =
            100.0 * (shifted - nominal) / nominal;

        std::cout << std::setw(3) << (i + 1)
                  << std::setw(10) << nominal
                  << std::setw(10) << shifted
                  << std::setw(14) << relChange << std::endl;
    }

   std::cout<<hObsInt->Integral()
<<","<<true_mult->Integral()<<","<<hObsInt->GetEntries()
<<","<<true_mult->GetEntries()<<std::endl;


    TMatrixD cov = computeCovarianceBinned(true_multResp, deltaLong, deltaShort,
                                       binsMult, edgesMult);
    std::vector<double> binWidths={1,1,1,1,1,2,2};
 TH2D* covXS=(TH2D*) convertCovCountsToXS( cov, flux, binWidths);
 
    for (int i=0; i<binsMult; i++) std::cout<<TMath::Sqrt(cov(i,i))<<std::endl;
   TFile output("output.root","UPDATE");
    covXS->SetName("covMatComb"); covXS->SetTitle("covMatComb");
    covXS->Write(); output.Write(); output.Close();
    return 0;
}


