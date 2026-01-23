#ifndef FLUXCOVARIANCEREWEIGHT_H
#define FLUXCOVARIANCEREWEIGHT_H

#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TH1D.h"

#include <vector>

// -------------------------------------
struct FluxBin {
    int pdg;
    double elo;
    double ehi;

    bool Contains(int p, double e) const {
        return (p == pdg && e >= elo && e < ehi);
    }
};

// -------------------------------------
class FluxCovarianceReweight {
public:
    FluxCovarianceReweight() {}

    bool LoadBinning(const char* binfile);
    bool LoadCovariance();
    void AddDiagonalEpsilon(double rel = 1e-8);
    double GetWeightFast(int ithrow, int binIndex) const;
    int    GetBinIndex(int pdg, double Enu) const;
    void FillAllThrowsForBin(int binIndex, Double_t* out) const;
      void LoadFluxHistograms(
        TH1D* h_nue,
        TH1D* h_nuebar,
        TH1D* h_numu,
        TH1D* h_numubar
      );
    std::vector<double> BuildFluxBinWidthVector(
    TH1D* h_nue,
    TH1D* h_nuebar,
    TH1D* h_numu,
    TH1D* h_numubar) const;
    TVectorD BuildNominalFluxVector(
    TH1D* h_nue,
    TH1D* h_nuebar,
    TH1D* h_numu,
    TH1D* h_numubar) const;
    bool GenerateThrows(int nThrows, int seed = 0);
    bool GenerateThrows2(int nThrows, int seed = 0);
    const TVectorD& GetAllWeights(int pdg, double Enu) const;
    int FindBinIndex(int pdg, double Enu) const;
    double GetWeight(int ithrow, int pdg, double Enu) const;

  const TVectorD& GetNominalFlux() const { return fNominalFlux; }
  const std::vector<double>& GetFluxBinWidth() const { return fFluxBinWidth; }

  double GetIntegratedFluxNominal() const { return fIntegratedFluxNominal; }
  const std::vector<double>& GetIntegratedFluxThrows() const {
    return fIntegratedFluxThrows;
  }

  double GetFluxWeight(int fluxBin, int universe) const;
  void ComputeIntegratedFluxes();
private:

    std::vector<FluxBin> fBins;
    std::vector<TVectorD> fThrows;
    TMatrixDSym fCov;

      // Nominal flux (150 bins)
      TVectorD fNominalFlux;
    
      // Bin widths (150 bins)
      std::vector<double> fFluxBinWidth;
    
      // Integrated flux
      double fIntegratedFluxNominal;
    
      // Integrated flux per universe
      std::vector<double> fIntegratedFluxThrows;
};

#endif
