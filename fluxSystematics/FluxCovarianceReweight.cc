#include "FluxCovarianceReweight.h"

#include "TFile.h"
#include "TH2D.h"
#include "TRandom3.h"
#include "TDecompChol.h"
#include "TMatrixDSymEigen.h"
#include "TH1D.h"

#include <fstream>
#include <iostream>
#include <cmath>

// =====================================================
// Load binning
// =====================================================
bool FluxCovarianceReweight::LoadBinning(const char* binfile)
{
    std::ifstream fin("flux_covariance_binning_NuMI_GeV.txt");
    if (!fin.is_open()) {
        std::cerr << "ERROR: cannot open binning file\n";
        return false;
    }

    int isRHC, pdg;
    double elo, ehi;

    while (fin >> isRHC >> pdg >> elo >> ehi) {
        FluxBin b;
        b.pdg = pdg;
        b.elo = elo;
        b.ehi = ehi;
        fBins.push_back(b);
    }

    std::cout << "Loaded " << fBins.size() << " flux bins\n";
    return true;
}

// =====================================================
// Load covariance matrix
// =====================================================
bool FluxCovarianceReweight::LoadCovariance()
{
    TFile* f = TFile::Open("numiFluxSyst.root", "READ");
    if (!f || f->IsZombie()) return false;

    TH2D* h = dynamic_cast<TH2D*>(f->Get("covariance_matrices/hadron/total/hcov_total"));


    if (!h) return false;

    int n = h->GetNbinsX();
    if ((int)fBins.size() != n) {
        std::cerr << "ERROR: binning size does not match covariance\n";
        return false;
    }

    fCov.ResizeTo(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            fCov(i,j) = h->GetBinContent(i+1, j+1);



    return true;
}

// =====================================================
// Generate correlated throws
// =====================================================
bool FluxCovarianceReweight::GenerateThrows(int nThrows, int seed)
{
    TRandom3 rng(seed);
    AddDiagonalEpsilon();
    TDecompChol chol(fCov);
    if (!chol.Decompose()) {
        std::cerr << "ERROR: Cholesky failed\n";
        return false;
    }

    TMatrixD* U = (TMatrixD*)chol.GetU().Clone();
    U->Transpose(chol.GetU()); // L = U^T

    int n = fCov.GetNrows();
    fThrows.clear();
    fThrows.reserve(nThrows);

    for (int t = 0; t < nThrows; ++t) {
        TVectorD z(n);
        for (int i = 0; i < n; ++i)
            z(i) = rng.Gaus(0,1);

        TVectorD delta = *U * z;

        // Convert to reweights: w = 1 + delta
        for (int i = 0; i < n; ++i)
            delta(i) = 1.0 + delta(i);

        fThrows.push_back(delta);
    }

    std::cout << "Generated " << nThrows << " throws\n";
    ComputeIntegratedFluxes();

    return true;
}

// =====================================================
// Find bin index
// =====================================================
int FluxCovarianceReweight::FindBinIndex(int pdg, double Enu) const
{
    for (size_t i = 0; i < fBins.size(); ++i)
        if (fBins[i].Contains(pdg, Enu))
            return i;

    return -1;
}

// =====================================================
// Public lookup
// =====================================================
double FluxCovarianceReweight::GetWeight(int ithrow,
                                          int pdg,
                                          double Enu) const
{
    if (ithrow < 0 || ithrow >= (int)fThrows.size())
        return 1.0;

    int idx = FindBinIndex(pdg, Enu);
    if (idx < 0)
        return 1.0;

    return fThrows[ithrow](idx);
}
bool FluxCovarianceReweight::GenerateThrows2(int nThrows, int seed)
{
    TRandom3 rng(seed);

    const int n = fCov.GetNrows();

    // Eigen decomposition
    TMatrixDSymEigen eig(fCov);
    TVectorD eigenValues = eig.GetEigenValues();
    TMatrixD eigenVectors = eig.GetEigenVectors();

    // Build sqrt(D) with protection
    TMatrixD sqrtD(n, n);
    sqrtD.Zero();

    int nClipped = 0;
    for (int i = 0; i < n; ++i) {
        if (eigenValues(i) > 0) {
            sqrtD(i,i) = std::sqrt(eigenValues(i));
        } else {
            sqrtD(i,i) = 0.0;  // clip negative modes
            nClipped++;
        }
    }

    if (nClipped > 0) {
        std::cout << "WARNING: clipped " << nClipped
                  << " non-positive eigenvalues\n";
    }

    // Transformation matrix
    TMatrixD A = eigenVectors * sqrtD;

    fThrows.clear();
    fThrows.reserve(nThrows);

    for (int t = 0; t < nThrows; ++t) {

        // Uncorrelated standard normals
        TVectorD z(n);
        for (int i = 0; i < n; ++i)
            z(i) = rng.Gaus(0,1);

        // Correlated variation
        TVectorD delta = A * z;

        // Convert to reweights
        for (int i = 0; i < n; ++i)
            delta(i) = 1.0 + delta(i);

        fThrows.push_back(delta);
    }
       ComputeIntegratedFluxes();

    std::cout << "Generated " << nThrows << " throws using eigen decomposition\n";
    return true;
}
const TVectorD& FluxCovarianceReweight::GetAllWeights(int pdg,
                                                       double Enu) const
{
    static TVectorD unity;  // fallback

    int idx = FindBinIndex(pdg, Enu);
    if (idx < 0) {
        unity.ResizeTo(fThrows.size());
        unity = 1.0;
        return unity;
    }

    static TVectorD weights;
    weights.ResizeTo(fThrows.size());

    for (size_t i = 0; i < fThrows.size(); ++i)
        weights(i) = fThrows[i](idx);

    return weights;
}
void FluxCovarianceReweight::AddDiagonalEpsilon(double rel)
{
    double maxDiag = 0.0;
    for (int i = 0; i < fCov.GetNrows(); ++i)
        maxDiag = std::max(maxDiag, fCov(i,i));

    double eps = rel * maxDiag;

    for (int i = 0; i < fCov.GetNrows(); ++i)
        fCov(i,i) += eps;

    std::cout << "Added diagonal epsilon = " << eps << "\n";
}
int FluxCovarianceReweight::GetBinIndex(int pdg, double Enu) const
{
    return FindBinIndex(pdg, Enu);
}

double FluxCovarianceReweight::GetWeightFast(int ithrow,
                                               int binIndex) const
{
    if (ithrow < 0 || ithrow >= (int)fThrows.size())
        return 1.0;

    if (binIndex < 0 || binIndex >= fThrows[ithrow].GetNrows())
        return 1.0;

    return fThrows[ithrow](binIndex);
}
void FluxCovarianceReweight::FillAllThrowsForBin(int binIndex,
                                                  Double_t* out) const
{
    const int nThrows = fThrows.size();

    for (int ithrow = 0; ithrow < nThrows; ++ithrow)
        out[ithrow] = 1.0;

    if (binIndex < 0)
        return;

    for (int ithrow = 0; ithrow < nThrows; ++ithrow) {
        if (binIndex < fThrows[ithrow].GetNrows())
            out[ithrow] = fThrows[ithrow](binIndex);
    }
}


std::vector<double> FluxCovarianceReweight::BuildFluxBinWidthVector(
    TH1D* h_nue,
    TH1D* h_nuebar,
    TH1D* h_numu,
    TH1D* h_numubar
) const {
    std::vector<double> w;
    w.reserve(150);

    for (int i = 1; i <= h_nue->GetNbinsX(); ++i)
        w.push_back(h_nue->GetBinWidth(i));

    for (int i = 1; i <= h_nuebar->GetNbinsX(); ++i)
        w.push_back(h_nuebar->GetBinWidth(i));

    for (int i = 1; i <= h_numu->GetNbinsX(); ++i)
        w.push_back(h_numu->GetBinWidth(i));

    for (int i = 1; i <= h_numubar->GetNbinsX(); ++i)
        w.push_back(h_numubar->GetBinWidth(i));

    return w;
}

TVectorD FluxCovarianceReweight::BuildNominalFluxVector(
    TH1D* h_nue,
    TH1D* h_nuebar,
    TH1D* h_numu,
    TH1D* h_numubar
) const{
    TVectorD fluxNom(150);
    int bin = 0;

    for (int i = 1; i <= h_nue->GetNbinsX(); ++i)
        fluxNom[bin++] = h_nue->GetBinContent(i);

    for (int i = 1; i <= h_nuebar->GetNbinsX(); ++i)
        fluxNom[bin++] = h_nuebar->GetBinContent(i);

    for (int i = 1; i <= h_numu->GetNbinsX(); ++i)
        fluxNom[bin++] = h_numu->GetBinContent(i);

    for (int i = 1; i <= h_numubar->GetNbinsX(); ++i)
        fluxNom[bin++] = h_numubar->GetBinContent(i);

    return fluxNom;
}
double FluxCovarianceReweight::GetFluxWeight(
    int binIndex,
    int ithrow
) const {
    return fThrows[ithrow](binIndex) / fNominalFlux[binIndex];
}
void FluxCovarianceReweight::LoadFluxHistograms(
    TH1D* h_nue,
    TH1D* h_nuebar,
    TH1D* h_numu,
    TH1D* h_numubar
) {
    const int NBINS = 150;

    fNominalFlux.ResizeTo(NBINS);
    fFluxBinWidth.clear();
    fFluxBinWidth.reserve(NBINS);

    int bin = 0;

    // nue (25)
    for (int i = 1; i <= h_nue->GetNbinsX(); ++i, ++bin) {
        fNominalFlux[bin] = h_nue->GetBinContent(i);
        fFluxBinWidth.push_back(h_nue->GetBinWidth(i));
    }

    // nuebar (25)
    for (int i = 1; i <= h_nuebar->GetNbinsX(); ++i, ++bin) {
        fNominalFlux[bin] = h_nuebar->GetBinContent(i);
        fFluxBinWidth.push_back(h_nuebar->GetBinWidth(i));
    }

    // numu (50)
    for (int i = 1; i <= h_numu->GetNbinsX(); ++i, ++bin) {
        fNominalFlux[bin] = h_numu->GetBinContent(i);
        fFluxBinWidth.push_back(h_numu->GetBinWidth(i));
    }

    // numubar (50)
    for (int i = 1; i <= h_numubar->GetNbinsX(); ++i, ++bin) {
        fNominalFlux[bin] = h_numubar->GetBinContent(i);
        fFluxBinWidth.push_back(h_numubar->GetBinWidth(i));
    }
}
void FluxCovarianceReweight::ComputeIntegratedFluxes()
{
    const int nUniv = (int)fThrows.size();
    const int nBins = fThrows[1].GetNrows();

    // Nominal
    fIntegratedFluxNominal = 0.0;
    for (int b = 50; b < nBins; ++b)
        fIntegratedFluxNominal +=
            fNominalFlux[b];// * fFluxBinWidth[b];

    // Universes
    fIntegratedFluxThrows.assign(nUniv, 0.0);

    for (int u = 0; u < nUniv; ++u) {
        for (int b = 50; b < nBins; ++b) {
            fIntegratedFluxThrows[u] +=
                fThrows[u](b)* fNominalFlux[b];// * fFluxBinWidth[b];
        }
    }
}