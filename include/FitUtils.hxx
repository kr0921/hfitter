#ifndef FitUtils_h
#define FitUtils_h

#include <iostream>
#include <vector>
#include <map>

// psyche
#include "BasicTypes.hxx"

// ROOT
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>

#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooWorkspace.h>
#include <RooFitResult.h>
#include <RooAbsPdf.h>
#include <RooAbsData.h>
#include <RooMCStudy.h>

#include <RooStats/HypoTestInverterResult.h>

namespace FitUtils{

  /// Get POT from file
  Float_t GetPOT(std::string file, bool dumpPOT = false);

  /// Check if tree exist in file
  bool CheckTree(TString InputFile, TString InputTree);

  /// Check if histogram has entries in one bin only 
  bool IsSingleBinHisto(TH1* histo);

  /// Check if two vectors have the same size
  template<typename vecT> 
  bool CheckVectorSize(const std::vector<vecT> vec1, const std::vector<vecT> vec2){ return (vec1.size() == vec2.size()); }

  /// Print best-fit for each POI, error and error/best-fit
  void PrintFractionalError(RooArgList paramsfit);

  /// Print (theta - theta0)DTheta for each nuisance parameter
  void PrintNuisanceParamsImpact(RooArgList paramsinit, RooArgList paramsfit);

  /// Get systematic histogram containing up and down events
  TH1* GetSystematicHistoFromNominal(TH1* nominal, TH1* syst, TString name); 

  /// Split string into substrings separated by sep
  std::vector<std::string> SplitString(std::string input, std::string sep);

  /// Set interpolation code
  void SetInterpolationCode(RooWorkspace* ws, Int_t code);

  /// Save workspace snapshot
  void SaveSnapshot(RooWorkspace* ws, TString snapshotname);

  /// Load workspace snapshot
  bool LoadSnapshot(RooWorkspace* ws, TString snapshotname);

  /// Save workspace in ROOT file
  void SaveWorkspace(RooWorkspace* ws, TString outFileName, bool recreate=true);

  /// Get fit covariance matrix
  TH2D* GetFitCovariance(RooFitResult* result);

  /// Get fit correlation matrix
  TH2* GetFitCorrelationMatrix(RooFitResult* result);

  /// Get fit correlation matrix and print highly correlated parameters
  TH2* PrintHighCorrelations(RooFitResult* result, Double_t threshold);

  /// Fit data. Results are stored in RooFitResult container
  RooFitResult* FitData(RooWorkspace* ws, std::string minimizer = "Minuit2", Int_t fitstrategy = 1, bool minoserror = false);

  /// Fit a given dataset. Results are stored in RooFitResult container
  RooFitResult* FitToyData(RooWorkspace* ws, RooAbsData* obsdata, std::string minimizer = "Minuit2", Int_t fitstrategy = 1, bool minoserror = false);

  /// Fit Asimov data. Results are stored in RooFitResult container
  RooFitResult* FitAsimovData(RooWorkspace* ws, std::string minimizer = "Minuit2", Int_t fitstrategy = 1, bool minoserror = false);

  /// Generate a toy MC
  RooAbsData* GenerateToyMC(RooWorkspace* ws, bool datanorm);

  /// Generate and fit a number ot toy experiments. Fit results are stored in tree 
  TTree* GenerateAndFit(RooWorkspace* ws, Int_t nexp, std::string minimizer = "Minuit2", Int_t fitstrategy = 1, bool minoserror = false, Double_t muVal = -1.0);

  /// Redo fit by setting all nuisance parameters to constant and shifting those in params by +/- 1 sigma one each time
  TTree* FitDataWithShiftedNuisanceParameters(RooWorkspace* ws, TString wsSnapshot, std::vector<TString> params);

  /// Save results of a MC study container in a tree
  TTree* RooMCStudyToTTree(RooMCStudy* mc);

  /// Write fit results in a tree
  TTree* RooFitResultToTTree(RooWorkspace* ws, RooFitResult* res);

  /// Write a RooDataSet in a tree
  TTree* RooDataSetToTTree(RooAbsData* data, TString treename);

  /// Set parameters in paramsPair to a new value and make them constant
  void MakeParametersConstant(RooWorkspace* ws, std::map<TString, Double_t> paramsPair);

  /// Reset values of the uncertainties in parList to nominal 
  void ResetValues(RooWorkspace* ws, const RooArgList& parList);

  /// Reset values of nominal parameters in parSet
  void ResetValuesToNominal(RooWorkspace* ws, const RooArgSet& parSet);

  /// Reset error on parameters in parList
  void ResetError(RooWorkspace* ws, const RooArgList& parList);

  // Get list of floating parameters
  RooArgList GetFloatingParsList(RooWorkspace* ws);

  // Get set of global observables
  const RooArgSet* GetGlobalObservablesSet(RooWorkspace* ws);

  /// Vector of plots with data and pdfs
  std::vector<TCanvas*> PlotDatasetsAndPdfs(RooWorkspace *work, TString name, TString error, TString plottodraw, TString binNames, RooAbsData* data=NULL, RooFitResult* result=NULL);

  /// Remove empty bins from RooPlot
  void RemoveEmptyBins(RooPlot* frame);

  /// Vector of plots for the NLL
  std::vector<TCanvas*> PlotNLL(RooWorkspace *work, TString name, RooFitResult* result, bool plotPLL=false);

  /// Plot nuisance parameters
  TCanvas* PlotNuisanceParameters(TTree* tree, RooWorkspace* ws);

}


#endif
