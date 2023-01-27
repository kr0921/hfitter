#include "FitUtils.hxx"

// ROOT
#include <TFile.h>
#include <TTree.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLine.h>
#include <THStack.h>

#include <RooDataSet.h>
#include <RooSimultaneous.h>
#include <RooArgSet.h>
#include <RooCategory.h>
#include <RooRealIntegral.h>
#include <RooPlot.h>
#include <RooHist.h>
#include <Roo1DTable.h>
#include <RooRealSumPdf.h>
#include <RooProduct.h>

#include <RooStats/HistFactory/PiecewiseInterpolation.h>
#include <RooStats/ProfileLikelihoodTestStat.h>
#include <RooStats/ToyMCSampler.h>
#include <RooStats/ModelConfig.h>

// psyche
#include "Header.hxx"
//#include "SelectionBase.hxx"
//#include "ND280AnalysisUtils.hxx"
//#include "ConfigurationManager.hxx"


//********************************************************************
Float_t FitUtils::GetPOT(std::string file, bool dumpPOT){
  //********************************************************************

  TFile *_file    = new TFile(file.c_str(),"READ");
  if(!_file) return -1.0;
  if(_file->GetNkeys() <= 0){ 
    std::cerr << "ERROR::Empty file " << file.c_str() << ". POT not calculated." << std::endl;
    _file->Close();
    return -1.0;
  }

  TTree *_pot     = (TTree*)_file->Get("header");

  Header header1;
  header1.AddHeader(_pot);
  Float_t mmpot = header1.GetPOT();
  if(dumpPOT)
    header1.DumpPOT();

  //POTDef def200  = k200KA;
  //POTDef def250  = k250KA;
  //POTDef defm250 = km250KA;

  _file->Close();

  return mmpot;

}

//********************************************************************
bool FitUtils::CheckTree(TString InputFile, TString InputTree){
  //********************************************************************

  TFile* f = new TFile(InputFile.Data(), "READ");
  bool checkttree =  f->GetListOfKeys()->Contains(InputTree.Data());
  f->Close();

  return checkttree;

}

//********************************************************************
bool FitUtils::IsSingleBinHisto(TH1* histo){
  //********************************************************************

  Int_t nbinswithevents = 0;
  for(Int_t j = 1; j <= histo->GetNbinsX(); j++){
    if(histo->GetBinContent(j) > 0.0)
      nbinswithevents++;
  }

  if(nbinswithevents == 1) return true;

  return false;

}

//********************************************************************
void FitUtils::PrintFractionalError(RooArgList paramsfit){
  //********************************************************************

  std::cout << std::endl;
  std::cout << "INFO::Printing fractional error for each fit parameter:-" << std::endl;
  TIterator* iter = paramsfit.createIterator();
  RooRealVar* var(0);
  for(Int_t i=0; (var = (RooRealVar*)iter->Next()); ++i){
    TString varname(var->GetName());
    if(varname.Contains("alpha_") || varname.Contains("gamma_")) continue;
    std::cout << varname.Data() << ": Fit result = " << var->getVal() << " +/- " << var->getError() << ". Fractional error =  " <<  var->getError()/var->getVal() << std::endl;
  }

  delete iter;

}

//********************************************************************
void FitUtils::PrintNuisanceParamsImpact(RooArgList paramsinit, RooArgList paramsfit){
  //********************************************************************

  std::cout << std::endl;
  std::cout << "INFO::Printing (theta - theta0)/Dtheta:-" << std::endl;
  TIterator* Inititer = paramsinit.createIterator();
  RooRealVar* Initvar(0);
  for(Int_t i=0; (Initvar = (RooRealVar*)Inititer->Next()); ++i){
    TString Initvarname(Initvar->GetName());
    if(!Initvarname.Contains("alpha_")) continue;
    Double_t Initvarerror = Initvar->getError();
    // Protection
    if(Initvarerror == 0)
      Initvarerror = 1.0;

    TIterator* iter = paramsfit.createIterator();
    RooRealVar* var(0);
    for(Int_t j=0; (var = (RooRealVar*)iter->Next()); ++j){
      TString varname(var->GetName());
      if(varname != Initvarname) continue;

      std::cout << varname.Data() << " = " << (var->getVal() - Initvar->getVal())/Initvarerror << std::endl;
    }
    delete iter;
  }
}

//********************************************************************
TH1* FitUtils::GetSystematicHistoFromNominal(TH1* nominal, TH1* syst, TString name){
  //********************************************************************
  
  TString histoname = TString(syst->GetName()) + TString("_") + name;
  TH1F* histo = new TH1F(histoname.Data(), histoname.Data(), nominal->GetNbinsX(), nominal->GetBinLowEdge(1), nominal->GetBinLowEdge(nominal->GetNbinsX()+1));
  if(nominal->GetNbinsX() != syst->GetNbinsX()){
    std::cout << "WARNING::Nominal and systematic histograms does not have the same number of bins. No systematic is applied!" << std::endl;
    return histo;
  }

  for(Int_t i=0; i < nominal->GetNbinsX(); i++){
    Double_t nom = nominal->GetBinContent(i+1);
    Double_t sys = syst->GetBinContent(i+1);

    if(nom <= 0.){
      std::cout << "WARNING::No events in histogram " << nominal->GetName() << " for bin " << i << ". Skip systematics propagation for this bin and systematic histogram " << syst->GetName() << std::endl;
      histo->SetBinContent(i+1, 0.0);
      continue;
    }

    if(sys < 0.)
      std::cout << "WARNING::Negative fractional error for histogram " << syst->GetName() << " and bin " << i << std::endl;

    // Avoid large fluctuation in the systematics. Could be caused if the number of events in a bin is very small
    if(sys > 1.0){
      std::cout << "WARNING::Fractional error for histogram " << syst->GetName() << " and bin " << i << " is larger than 100%. Set it to 100%." << std::endl;
      sys = 1.0;
    }

    Float_t w = 1.0;
    if(name == "up" || name == "Up" || name == "UP" || name == "UP2" ){
      w = nom + nom*sys;
    }
    else if(name == "down" || name == "Down" || name == "DOWN" || name == "DOWN2"){
      w = nom - nom*sys;
    }

    if(w < 0.) w = 0.;

    histo->SetBinContent(i+1, w);
  }

  return histo;

}

//********************************************************************
std::vector<std::string> FitUtils::SplitString(std::string input, std::string sep){
  //********************************************************************

  std::vector<std::string> string_vector;

  while(input.find_first_of(sep) < input.length()){
    string_vector.push_back( input.substr(0, input.find_first_of(sep)) );
    input = input.substr(input.find_first_of(sep)+1, input.length());
  }
  string_vector.push_back(input);

  return string_vector;

}

//********************************************************************
void FitUtils::SetInterpolationCode(RooWorkspace* ws, Int_t code){
  //********************************************************************

  if(!ws){
    std::cout << "ERROR::NULL workspace. No interpolation code added." << std::endl;
    return;
  }

  RooArgSet fun = ws->allFunctions();
  TIterator* iter = fun.createIterator();

  RooAbsArg* arg(0);
  while( (arg=(RooAbsArg*)iter->Next()) ) {
    if(arg->ClassName()!=TString("PiecewiseInterpolation") ) continue;
    
    PiecewiseInterpolation* p = (PiecewiseInterpolation*)ws->function(arg->GetName() );
    

   // p->setAllInterpCodes(code);
  
  }

  delete iter;

}

//********************************************************************
void FitUtils::SaveSnapshot(RooWorkspace* ws, TString snapshotname){
  //********************************************************************

  if(!ws){
    std::cerr << "ERROR::Workspace not found. Will not save snapshot." << std::endl;
    return;
  }

  RooSimultaneous* simpdf = (RooSimultaneous*)ws->pdf("simPdf");
  if(!simpdf)
    simpdf = (RooSimultaneous*)ws->pdf("combPdf");
  if(!simpdf){
    std::cerr << "ERROR::Pdf not found in workspace. Will not save snapshot." << std::endl;
    return;
  }

  RooAbsData* obsdata = (RooAbsData*)ws->data("obsData");
  RooArgSet* parameters = (RooArgSet*)simpdf->getParameters(*obsdata);
  if(!ws->loadSnapshot(snapshotname.Data())){
    ws->saveSnapshot(snapshotname.Data(),*parameters);
  }
  else{
    std::cout << "WARNING::Snapshot " << snapshotname.Data() << " found in workspace. Will not overwrite it." << std::endl;
  }

  //gDirectory->Add(ws);

}

//********************************************************************
bool FitUtils::LoadSnapshot(RooWorkspace* ws, TString snapshotname){
  //********************************************************************

  if(!ws){
    std::cerr << "ERROR::Workspace not found. Can't find snapshot." << std::endl;
    return false;
  }

  if(!ws->loadSnapshot(snapshotname)){
    std::cerr << "ERROR::Snapshot not found in workspace." << std::endl;
    return false;
  }

  ws->loadSnapshot(snapshotname);
  std::cout << "INFO::Workspace snapshot " << snapshotname.Data() << " loaded." << std::endl;

  return true;

}

//********************************************************************
void FitUtils::SaveWorkspace(RooWorkspace* ws, TString outFileName, bool recreate){
  //********************************************************************

  if(!ws){
    std::cerr << "ERROR::Workspace not found. Won't save." << std::endl;
    return;
  }

  ws->writeToFile(outFileName.Data(), recreate);
  std::cout << "INFO::Workspace written to file " << outFileName.Data() << std::endl;

  return;

}

//********************************************************************
TH2D* FitUtils::GetFitCovariance(RooFitResult* result){
  //********************************************************************

  RooArgList floatParsFinal = result->floatParsFinal();
  const TMatrixDSym& covmatrix = result->covarianceMatrix();
  Int_t n = covmatrix.GetNcols();
  TH2D* CovarianceHisto = new TH2D("FitCovariance", "FitCovariance", n,0,n, n,0,n);
  CovarianceHisto->GetXaxis()->SetLabelSize(0.01);
  CovarianceHisto->GetYaxis()->SetLabelSize(0.01);
  for(Int_t j = 0 ; j<n ; j++) {
    for(Int_t k = 0 ; k<n; k++) {
      CovarianceHisto->Fill(j+0.5,n-k-0.5, covmatrix(j,k));
    }
    CovarianceHisto->GetXaxis()->SetBinLabel(j+1, floatParsFinal.at(j)->GetName());
    CovarianceHisto->GetYaxis()->SetBinLabel(n-j, floatParsFinal.at(j)->GetName());
  }
  CovarianceHisto->SetMinimum(-1);
  CovarianceHisto->SetMaximum(+1);

  return CovarianceHisto;

}

//********************************************************************
TH2* FitUtils::GetFitCorrelationMatrix(RooFitResult* result){
  //********************************************************************

  TH2* corrhisto = result->correlationHist();
  corrhisto->GetXaxis()->SetLabelSize(0.01);
  corrhisto->GetYaxis()->SetLabelSize(0.01);

  return corrhisto;

}

//********************************************************************
TH2* FitUtils::PrintHighCorrelations(RooFitResult* result, Double_t threshold){
  //********************************************************************

  TH2* corrhisto = FitUtils::GetFitCorrelationMatrix(result);
  for(Int_t i = 1; i <= corrhisto->GetNbinsY(); i++){
    for(Int_t j = 1; j <= corrhisto->GetNbinsX() && j <= (corrhisto->GetNbinsX()-i); j++){
      if(fabs(corrhisto->GetBinContent(i,j)) > threshold)
	std::cout << "INFO::High correlation coefficient " << corrhisto->GetBinContent(i,j)
		  << " , for parameters " << corrhisto->GetXaxis()->GetBinLabel(j)
		  << " - " << corrhisto->GetYaxis()->GetBinLabel(i)  << std::endl;
    }
  }

  return corrhisto;

}

//********************************************************************
RooFitResult* FitUtils::FitData(RooWorkspace* ws, std::string minimizer, Int_t fitstrategy, bool minoserror){
  //********************************************************************

  if(!ws){
    std::cerr << "ERROR::NULL Workspace. Return empty fit!!" << std::endl;
    return NULL;
  }

  // Get the configuration model from workspace
  RooStats::ModelConfig *combined_config = (RooStats::ModelConfig*)ws->obj("ModelConfig");
  if(!combined_config){
    std::cerr << "ERROR::No model config " << "ModelConfig " << " in workspace. Return empty fit!!" << std::endl;
    return NULL;
  }

  // Silence the output
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

  // Dataset
  RooAbsData *obsdata = ws->data("obsData");

  RooFitResult* fresult = combined_config->GetPdf()->fitTo(*obsdata, RooFit::Minimizer(minimizer.c_str()), RooFit::Save(), RooFit::Strategy(fitstrategy), RooFit::Minos(minoserror), RooFit::PrintLevel(1), RooFit::Constrain(*combined_config->GetPdf()->getParameters(*obsdata)), RooFit::Extended(), RooFit::Offset(true) );

  if(fresult){
    fresult->SetName("FitData_result");
    ws->import(*fresult);
    //gDirectory->Add(fresult);
  }

  // Save post-fit snapshot
  RooArgSet* params = (RooArgSet*)combined_config->GetPdf()->getParameters(*obsdata);
  ws->saveSnapshot("snapshot_datafit",*params);

  // Reset the verbosity
  RooMsgService::instance().reset();

  return fresult;

}

//********************************************************************
RooFitResult* FitUtils::FitAsimovData(RooWorkspace* ws, std::string minimizer, Int_t fitstrategy, bool minoserror){
  //********************************************************************

  if(!ws){
    std::cerr << "ERROR::NULL Workspace. Return empty fit!!" << std::endl;
    return NULL;
  }

  // Get the configuration model from workspace
  RooStats::ModelConfig *combined_config = (RooStats::ModelConfig*)ws->obj("ModelConfig");
  if(!combined_config){
    std::cerr << "ERROR::No model config " << "ModelConfig " << " in workspace. Return empty fit!!" << std::endl;
    return NULL;
  }

  // Dataset
  RooAbsData *obsdata = ws->data("asimovData");

  RooFitResult* fresult = combined_config->GetPdf()->fitTo(*obsdata, RooFit::Minimizer(minimizer.c_str()), RooFit::Save(), RooFit::Strategy(fitstrategy), RooFit::Minos(minoserror), RooFit::Constrain(*combined_config->GetPdf()->getParameters(*obsdata)), RooFit::Extended(), RooFit::Offset(true));

  if(fresult){
    fresult->SetName("FitAsimovData_result");
    ws->import(*fresult);
    //gDirectory->Add(fresult);
  }

  // Save post-fit snapshot
  RooArgSet* params = (RooArgSet*)combined_config->GetPdf()->getParameters(*obsdata);
  ws->saveSnapshot("snapshot_asimovfit",*params);

  return fresult;

}

//********************************************************************
RooFitResult* FitUtils::FitToyData(RooWorkspace* ws, RooAbsData* obsdata, std::string minimizer, Int_t fitstrategy, bool minoserror){
  //********************************************************************

  if(!ws){
    std::cerr << "ERROR::NULL Workspace. Return empty fit!!" << std::endl;
    return NULL;
  }

  // Get the configuration model from workspace
  RooStats::ModelConfig *combined_config = (RooStats::ModelConfig*)ws->obj("ModelConfig");
  if(!combined_config){
    std::cerr << "ERROR::No model config " << "ModelConfig " << " in workspace. Return empty fit!!" << std::endl;
    return NULL;
  }

  RooFitResult* fresult = combined_config->GetPdf()->fitTo(*obsdata, RooFit::Minimizer(minimizer.c_str()), RooFit::Save(), RooFit::Strategy(fitstrategy), RooFit::Minos(minoserror), RooFit::PrintLevel(-1), RooFit::Constrain(*combined_config->GetPdf()->getParameters(*obsdata)), RooFit::Extended(), RooFit::Offset(true) );

  return fresult;

}

//********************************************************************
RooAbsData* FitUtils::GenerateToyMC(RooWorkspace* ws, bool datanorm){
  //********************************************************************

  if(!ws){
    std::cerr << "ERROR::NULL workspace. No dataset generated!" << std::endl;
    return NULL;
  }

  // Get the configuration model from workspace
  RooStats::ModelConfig *combined_config = (RooStats::ModelConfig*)ws->obj("ModelConfig");
  if(!combined_config){
    std::cerr << "ERROR::No model config " << " ModelConfig " << " in workspace!" << std::endl;
    return NULL;
  }

  // Get the pdf
  RooAbsPdf* pdf = combined_config->GetPdf();
  if(!pdf){
    std::cerr << "ERROR::No pdf found in ModelConfig. No dataset generated!" << std::endl;
    return NULL;
  }

  RooAbsData* data = ws->data("obsData");
  if(!data){
    std::cerr << "ERROR::No dataset found in workspace. No dataset generated!" << std::endl;
    return NULL;
  }

  const RooArgSet* obs = combined_config->GetObservables();
  if(!obs){
    std::cerr << "ERROR::No observables found in ModelConfig. No dataset generated!" << std::endl;
    return NULL;
  }

  // Silence the output
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

  RooAbsData* mctoy;
  if(datanorm){
    Int_t NEventsToGenerate = (Int_t)(data->sumEntries());
    mctoy = pdf->generate(*obs, RooFit::NumEvents(NEventsToGenerate), RooFit::AutoBinned(false));
  }
  else{
    mctoy = pdf->generate(*obs, RooFit::AutoBinned(false), RooFit::Extended());
  }

  // Reset verbosity
  RooMsgService::instance().reset();

  return mctoy;

}

//********************************************************************
TTree* FitUtils::GenerateAndFit(RooWorkspace* ws, Int_t nexp, std::string minimizer, Int_t fitstrategy, bool minoserror, Double_t muVal){
  //********************************************************************

  if(!ws){
    std::cerr << "ERROR::NULL Workspace. Return empty fit!!" << std::endl;
    return NULL;
  }

  // Get the configuration model from workspace
  RooStats::ModelConfig *combined_config = (RooStats::ModelConfig*)ws->obj("ModelConfig");
  if(!combined_config){
    std::cerr << "ERROR::No model config " << " ModelConfig " << " in workspace. Return empty fit!!" << std::endl;
    return NULL;
  }

  // Silence the output
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

  RooAbsPdf* pdf = combined_config->GetPdf();
  const RooArgSet* obs_set = combined_config->GetObservables();

  // Create mc study. Only used to store the fit result
  RooMCStudy* mcstudy = new RooMCStudy( *pdf, *obs_set, RooFit::FitOptions("r"));

  // change parameter of interest value
  if(muVal >= 0){
    RooRealVar *poiparam = (RooRealVar*)(combined_config->GetParametersOfInterest()->first());
    if(poiparam) poiparam->setVal(muVal);
  }

  // Create test statistics object
  RooStats::ProfileLikelihoodTestStat ts(*combined_config->GetPdf());

  // Create toy mc sampler
  RooStats::ToyMCSampler sampler(ts,nexp);
  sampler.SetPdf(*combined_config->GetPdf());
  sampler.SetObservables(*combined_config->GetObservables());
  sampler.SetGlobalObservables(*combined_config->GetGlobalObservables());
  sampler.SetParametersForTestStat(*combined_config->GetParametersOfInterest());

  RooArgSet poiAndNuisance;
  poiAndNuisance.add(*combined_config->GetParametersOfInterest());
  if(combined_config->GetNuisanceParameters())
    poiAndNuisance.add(*combined_config->GetNuisanceParameters());
  RooArgSet* nullParams = (RooArgSet*) poiAndNuisance.snapshot();
  //nullParams->Print();
  // Will be used as start values of fit
  ws->saveSnapshot("paramsToFit",poiAndNuisance);

  for(Int_t i=0; i < nexp; ++i) {
    if(i%10==0)
      std::cout << "INFO::Running on toy " << i << std::endl;

    // Reset starting values of fit
    ws->loadSnapshot("paramsToFit");

    RooAbsData* toyMC = sampler.GenerateToyData( *nullParams );
    toyMC->Print();

    RooFitResult* fresult = combined_config->GetPdf()->fitTo(*toyMC, RooFit::Minimizer(minimizer.c_str()), RooFit::Save(), RooFit::Strategy(fitstrategy), RooFit::Minos(minoserror), RooFit::Constrain(*combined_config->GetPdf()->getParameters(*toyMC)), RooFit::Extended(), RooFit::Offset(true), RooFit::PrintLevel(-1)); // RooFit::CloneData(kFALSE)

    mcstudy->addFitResult(*fresult);

    delete toyMC;
  }

  // Reset verbosity
  RooMsgService::instance().reset();

  TTree* mcstree = FitUtils::RooMCStudyToTTree(mcstudy);

  delete mcstudy;

  return mcstree;

}

//********************************************************************
TTree* FitUtils::FitDataWithShiftedNuisanceParameters(RooWorkspace* ws, TString wsSnapshot, std::vector<TString> params){
  //********************************************************************

  if(!ws){
    std::cerr << "ERROR::NULL Workspace. Return empty fit!!" << std::endl;
    return NULL;
  }

  // Get the configuration model from workspace
  RooStats::ModelConfig *combined_config = (RooStats::ModelConfig*)ws->obj("ModelConfig");
  if(!combined_config){
    std::cerr << "ERROR::No model config " << "ModelConfig " << " in workspace. Return empty fit!!" << std::endl;
    return NULL;
  }

  // Silence the output
  RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

  if(!FitUtils::LoadSnapshot(ws, wsSnapshot.Data())){
    std::cerr << "ERROR::Snapshot not found in workspace. Return empty fit!!" << std::endl;
    return NULL;
  }

  if(params.empty()){
    std::cout << "WARNING::List of nuisance parameters is empty. Nothing to do here." << std::endl;
    return NULL;
  }

  RooAbsPdf* pdf = combined_config->GetPdf();
  const RooArgSet* obs_set = combined_config->GetObservables();

  // Create mc study. Only used to store the fit result
  RooMCStudy* mcstudy = new RooMCStudy( *pdf, *obs_set, RooFit::FitOptions("r"));

  std::vector<Double_t> valVec; std::vector<TString> strvalVec; std::vector<Double_t> defvalVec;
  for(UInt_t i = 0; i < params.size(); i++){
    bool loadsnapshot = FitUtils::LoadSnapshot(ws, wsSnapshot.Data());

    if(!loadsnapshot){
      std::cout << "WARNING::Still can't find " << wsSnapshot.Data() << ". Check workspace." << std::endl;
      continue;
    }

    TString nuisparam = params[i];
    RooRealVar* var = (RooRealVar*)ws->var(nuisparam);
    if(!var){
      std::cout << "Parameter not found " << nuisparam << ". Skip fit." << std::endl;
      continue;
    }
    Double_t varvalue = var->getVal();
    Double_t valup    = varvalue + var->getErrorHi();
    Double_t valdown  = varvalue + var->getErrorLo();

    // Set all nuisance constants to their post-fit values
    const RooArgSet* nuisanceParams = ( (RooStats::ModelConfig*)ws->obj("ModelConfig") )->GetNuisanceParameters();
    RooStats::SetAllConstant(*nuisanceParams);

    // Dataset
    RooAbsData *obsdata = ws->data("obsData");

    std::cout << "INFO::Fitting data with all nuisance parameters constant and shifting " << nuisparam << " from nominal of " << varvalue << " to up and down values " << valup << " , " << valdown << std::endl;

    // Repeat fit twice, with parameter shifted by +/- 1sigma
    var->setVal(valup);
    RooFitResult* upresult = combined_config->GetPdf()->fitTo(*obsdata, RooFit::Minimizer("Minuit2"), RooFit::Save(), RooFit::Strategy(1), RooFit::Minos(false), RooFit::PrintLevel(-1), RooFit::Extended(), RooFit::Offset(true) );

    mcstudy->addFitResult(*upresult);
    defvalVec.push_back(varvalue);
    valVec.push_back(valup);
    strvalVec.push_back(nuisparam);

    var->setVal(valdown);
    RooFitResult* downresult = combined_config->GetPdf()->fitTo(*obsdata, RooFit::Minimizer("Minuit2"), RooFit::Save(), RooFit::Strategy(1), RooFit::Minos(false), RooFit::PrintLevel(-1), RooFit::Extended(), RooFit::Offset(true) );

    mcstudy->addFitResult(*downresult);
    defvalVec.push_back(varvalue);
    valVec.push_back(valdown);
    strvalVec.push_back(nuisparam);

  }

  // Reset the verbosity
  RooMsgService::instance().reset();

  TTree* mcstree = FitUtils::RooMCStudyToTTree(mcstudy);
  mcstree->SetNameTitle(Form("ShiftNuisanceTree_%s",wsSnapshot.Data()),Form("ShiftNuisanceTree_%s",wsSnapshot.Data()));

  delete mcstudy;

  // Add extra branches to the tree
  Double_t systDefaultValue; Double_t systShiftedValue; TString systName;
  TBranch *bsystDefaultValue = mcstree->Branch("systDefaultValue", &systDefaultValue, "systDefaultValue/D");
  TBranch *bsystShiftedValue = mcstree->Branch("systShiftedValue", &systShiftedValue, "systShiftedValue/D");
  TBranch *bsystName         = mcstree->Branch("systName",         &systName,         "systName/C");

  for(UInt_t i=0; i < mcstree->GetEntries(); i++){
    //mcstree->GetEntry(i);

    if(i > defvalVec.size()){
      std::cerr << "ERROR::TTree from RooMCStudy has more entries than the nuisance parameters vectors. Will stop filling the tree." << std::endl;
      break;
    }

    systDefaultValue = defvalVec[i];
    systShiftedValue = valVec[i];
    systName         = TString(strvalVec[i].Data());

    bsystDefaultValue->Fill();
    bsystShiftedValue->Fill();
    bsystName->Fill();
  }

  return mcstree;

}

//********************************************************************
TTree* FitUtils::RooMCStudyToTTree(RooMCStudy* mc){
  //********************************************************************

  TTree* myTree = new TTree("mcstree","mcstree");

  const RooDataSet& toydata = mc->fitParDataSet();

  // Fit status for Minuit2
  //status = 0 : OK
  //status = 1 : Covariance was made pos defined
  //status = 2 : Hesse is invalid
  //status = 3 : Edm is above max
  //status = 4 : Reached call limit
  //status = 5 : Any other failure

  Int_t covQual=0, status=0, numInvalidNLL=0;
  Float_t minNll=0, edm=0;
  myTree->Branch("covQual",        &covQual,        "covQual/I" );
  myTree->Branch("status",         &status,         "status/I" );
  myTree->Branch("numInvalidNLL",  &numInvalidNLL,  "numInvalidNLL/I" );
  myTree->Branch("minNll",         &minNll,         "minNll/F" );
  myTree->Branch("edm",            &edm,            "edm/F" );

  std::vector<Float_t> varVals;
  const RooArgSet* args = toydata.get();
  varVals.resize( args->getSize(), -999. );

  RooRealVar* var(0);
  TIterator* Itr = args->createIterator();
  for(Int_t i=0; (var = (RooRealVar*)Itr->Next()); ++i){
    TString varName = var->GetName();
    TString varNameF = TString(var->GetName()) + "/F";
    myTree->Branch( varName.Data(), &varVals[i], varNameF.Data() );
  }
  delete Itr;

  // Save correlations between parameters
  std::vector<Float_t> corrVals;
  const RooArgList& fitparams = mc->fitResult(0)->floatParsFinal();
  corrVals.resize( fitparams.getSize()*fitparams.getSize(), -999. );

  Int_t counter = 0;
  for(Int_t i=0; i < fitparams.getSize(); ++i){
    RooRealVar* fvar = (RooRealVar*)fitparams.at(i);
    TString varName = fvar->GetName();
    for(Int_t j=0; j < fitparams.getSize(); ++j){
      if(j == i) continue;
      RooRealVar* fvar2 = (RooRealVar*)fitparams.at(j);
      TString varName2 = varName + TString("_corr_") + TString(fvar2->GetName());
      TString varNameF = varName2 + "/F";
      myTree->Branch( varName2.Data(), &corrVals[counter], varNameF.Data() );
      counter++;
    }
  }

  // Fill the tree by looping over mcstudy
  for(int iToy = 0; iToy < toydata.numEntries(); iToy++){
    const RooFitResult* result = mc->fitResult(iToy);
    edm     = result->edm();
    covQual = result->covQual();
    status  = result->status();
    minNll  = result->minNll();
    numInvalidNLL = result->numInvalidNLL();

    toydata.get(iToy); // reset args to new value
    Itr = args->createIterator();
    for (Int_t i=0; (var=(RooRealVar*)Itr->Next()); ++i) { varVals[i] = var->getVal(); }
    delete Itr;

    const RooArgList& fparams = result->floatParsFinal();
    counter = 0;
    for(Int_t i=0; i < fparams.getSize(); ++i){
      RooRealVar* fvar = (RooRealVar*)fparams.at(i);
      for(Int_t j=0; j < fparams.getSize(); ++j){
        if(j == i) continue;
        RooRealVar* fvar2 = (RooRealVar*)fparams.at(j);
        corrVals[counter] = result->correlation(fvar->GetName(), fvar2->GetName());
        counter++;
      }
    }

    myTree->Fill();
  }

  return myTree;

}

//********************************************************************
TTree* FitUtils::RooFitResultToTTree(RooWorkspace* ws, RooFitResult* res){
  //********************************************************************

  if(!ws){
    std::cerr << "ERROR::NULL Workspace. Return NULL TTree!!" << std::endl;
    return NULL;
  }

  // Get the configuration model from workspace
  RooStats::ModelConfig *combined_config = (RooStats::ModelConfig*)ws->obj("ModelConfig");
  if(!combined_config){
    std::cerr << "ERROR::No model config " << "ModelConfig " << " in workspace. Return empty fit!!" << std::endl;
    return NULL;
  }

  RooAbsPdf* pdf = combined_config->GetPdf();
  const RooArgSet* obs_set = combined_config->GetObservables();

  // Create mc study. Only used to store the fit result
  RooMCStudy* mcstudy = new RooMCStudy( *pdf, *obs_set, RooFit::FitOptions("r"));
  mcstudy->addFitResult(*res);

  TTree* myTree = FitUtils::RooMCStudyToTTree(mcstudy);

  delete mcstudy;

  return myTree;

}

//********************************************************************
TTree* FitUtils::RooDataSetToTTree(RooAbsData* data, TString treename){
  //********************************************************************

  TTree* myTree = new TTree(treename, treename);

  std::vector<Float_t> varVals;
  const RooArgSet* args = data->get();
  varVals.resize( args->getSize(), -999. );

  RooRealVar* var(0);
  TIterator* Itr = args->createIterator();
  for(Int_t i=0; (var = (RooRealVar*)Itr->Next()); ++i) {
    TString varName = var->GetName();
    TString varNameF = TString(var->GetName()) + "/F";
    myTree->Branch( varName.Data(), &varVals[i], varNameF.Data() );
  }
  delete Itr;

  // Fill the tree by looping over mcstudy
  for(int iToy=0; iToy<data->numEntries(); iToy++){
    data->get(iToy); // reset args to new value
    TIterator* Itr2 = args->createIterator();
    for (Int_t i=0; (var=(RooRealVar*)Itr2->Next()); ++i) {varVals[i] = var->getVal();}
    delete Itr2;

    myTree->Fill();
  }

  return myTree;

}

//********************************************************************
void FitUtils::MakeParametersConstant(RooWorkspace* ws, std::map<TString, Double_t> paramsPair){
  //********************************************************************

  if(!ws){
    std::cout << "ERROR::NULL workspace. Will not fix parameters." << std::endl;
    return;
  }

  std::map<TString, Double_t>::iterator exp;
  for(exp = paramsPair.begin(); exp != paramsPair.end(); exp++){
    TString param_name = exp->first;
    Double_t param_val = exp->second;
    RooRealVar* var = (RooRealVar*)ws->var(param_name);
    if(var){
      std::cout << "INFO::Setting parameter " << param_name.Data() << " to constant value " << param_val << std::endl;
      var->setVal(param_val);
      var->setConstant(kTRUE);
    }
    else{
      std::cout << "WARNING::Can't find parameter " << param_name.Data() << " in workspace." << std::endl;
    }
  }

}

//********************************************************************
void FitUtils::ResetValues(RooWorkspace* ws, const RooArgList& parList){
  //********************************************************************

  TIterator* iter = parList.createIterator();
  RooAbsArg* arg;
  while((arg=(RooAbsArg*)iter->Next())){
    TString varname;
    if(arg->InheritsFrom("RooRealVar") && !arg->isConstant()){
      varname = TString(arg->GetName());
    }
    else{
      continue;
    }

    RooRealVar* var = ws->var(varname.Data());
    if(!var){
      std::cout << "WARNING::Can't find parameter " << varname.Data() << " in workspace. Will not reset value." << std::endl;
      continue;
    }

    if(varname == ""){
      std::cout << "WARNING::Empty variable name. Skipping reset value." << std::endl;
      continue;
    }
    else if(varname.Contains("alpha")){
      var->setVal(0.0);
    }
    else if(varname.Contains("gamma_stat")){
      var->setVal(1.0);
    }
    else{
      var->setVal(1.0);
    }
  }

  delete iter;

}

//********************************************************************
void FitUtils::ResetValuesToNominal(RooWorkspace* ws, const RooArgSet& parSet){
  //********************************************************************

  TIterator* iter = parSet.createIterator();
  RooAbsArg* arg;
  while((arg=(RooAbsArg*)iter->Next())){
    TString varname;
    if(arg->InheritsFrom("RooRealVar") && arg->isConstant()){
      varname = TString(arg->GetName());
    }
    else{
      continue;
    }

    RooRealVar* var = ws->var(varname.Data());
    if(!var){
      std::cout << "WARNING::Can't find parameter " << varname.Data() << " in workspace. Will not reset nominal value." << std::endl;
      continue;
    }

    if(varname == ""){
      std::cout << "WARNING::Empty variable name. Skipping reset value to nominal." << std::endl;
      continue;
    }
    else if(varname.BeginsWith("nom_gamma_stat")){
      var->setVal(1.0);
    }
    else if(varname == TString("nominalLumi")){
      var->setVal(1.0);
    }
    else if(varname.BeginsWith("nom")){
      var->setVal(0.0);
    }
    else{
      var->setVal(0.0);
    }
  }

  delete iter;

}

//********************************************************************
void FitUtils::ResetError(RooWorkspace* ws, const RooArgList& parList){
  //********************************************************************

  TIterator* iter = parList.createIterator();
  RooAbsArg* arg;
  while((arg=(RooAbsArg*)iter->Next())){
    TString varname;
    if(arg->InheritsFrom("RooRealVar") && !arg->isConstant()){
      varname = TString(arg->GetName());
    }
    else{
      continue;
    }

    RooRealVar* var = ws->var(varname.Data());
    if(!var){
      std::cout << "WARNING::Can't find parameter " << varname.Data() << " in workspace. Will not reset error." << std::endl;
      continue;
    }

    if(varname == ""){
      std::cout << "WARNING::Empty variable name. Skipping reset error." << std::endl;
      continue;
    }
    else if(varname.Contains("alpha")){
      var->setError(1.0);
      if(var->getMin() < var->getVal() - 6.) var->setMin(var->getVal() - 6.);
      if(var->getMax() > var->getVal() + 6.) var->setMax(var->getVal() + 6.);
    }
    else if(varname.Contains("gamma_stat")){
      // Constraint could be either Gaus or Poisson
      RooAbsReal* constraint = (RooAbsReal*) ws->obj(Form("%s_constraint",varname.Data()));
      if(!constraint){
	std::cout << "WARNING::Constraint for variable " << varname.Data() << " not found. Skip reset error." << std::endl;
        continue;
      }

      TString constraintString = TString(constraint->IsA()->GetName());
      if(constraintString == "") continue;
      else if(constraintString.Contains("RooGaussian")){
        RooAbsReal* ErrorVar = (RooAbsReal*)ws->obj(Form("%s_sigma",varname.Data()));
        if(!ErrorVar){
	  std::cout << "WARNING::Constraint type " << constraintString.Data() << " for variable " << varname.Data() << " not found. Skip reset error." << std::endl;
          continue;
        }
        Double_t err = ErrorVar->getVal();
        var->setError(err);
        if(var->getMin() < var->getVal() - 6*err) var->setMin(var->getVal() - 6*err);
        if(var->getMax() > var->getVal() + 6*err) var->setMax(var->getVal() + 6*err);
      }
      else if(constraintString.Contains("RooPoisson")){
        RooAbsReal* ErrorVar = (RooAbsReal*)ws->obj(Form("nom_%s",varname.Data()));
        if(!ErrorVar){
	  std::cout << "WARNING::Constraint type " << constraintString.Data() << " for variable " << varname.Data() << " not found. Skip reset error." << std::endl;
          continue;
        }
        Double_t err = 1/sqrt(ErrorVar->getVal());
        var->setError(err);
        if(var->getMin() < var->getVal() - 6*err) var->setMin(var->getVal() - 6*err);
        if(var->getMax() > var->getVal() + 6*err) var->setMax(var->getVal() + 6*err);
      }
      else{
	std::cout << "WARNING::Unknown constraint type " << constraintString.Data() << ". Set prefit uncertainty to 0.00001." << std::endl;
        var->setError(0.00001);
      }
    }

  }

  delete iter;

}

//********************************************************************
RooArgList FitUtils::GetFloatingParsList(RooWorkspace* ws){
  //********************************************************************

  RooArgList floatingParsList;

  // Get the configuration model from workspace
  RooStats::ModelConfig *combined_config = (RooStats::ModelConfig*)ws->obj("ModelConfig");
  if(!combined_config){
    std::cerr << "ERROR::No model config " << "ModelConfig " << " in workspace. Return empty list of floating parameters!!" << std::endl;
    return floatingParsList;
  }

  const RooArgSet* obsSet = combined_config->GetObservables();
  if (obsSet == 0) { return floatingParsList; }

  RooAbsPdf* pdf = combined_config->GetPdf();
  if (pdf == 0) { return floatingParsList; }

  const RooArgSet* pars = pdf->getParameters(obsSet);
  if (pars==0) { return floatingParsList; }

  TIterator* iter = pars->createIterator() ;
  RooAbsArg* arg ;
  while( (arg=(RooAbsArg*)iter->Next()) ) {
    if(arg->InheritsFrom("RooRealVar") && !arg->isConstant()){
      floatingParsList.add( *arg );
    }
  }
  delete iter;

  return floatingParsList;

}

//********************************************************************
const RooArgSet* FitUtils::GetGlobalObservablesSet(RooWorkspace* ws){
  //********************************************************************

  // Get the configuration model from workspace
  RooStats::ModelConfig *combined_config = (RooStats::ModelConfig*)ws->obj("ModelConfig");

  if(!combined_config){
    std::cerr << "ERROR::No model config " << "ModelConfig " << " in workspace. Return empty set of global parameters!!" << std::endl;
    return NULL;
  }

  return (combined_config->GetGlobalObservables());
}

//********************************************************************
std::vector<TCanvas*> FitUtils::PlotDatasetsAndPdfs(RooWorkspace *work, TString name, TString error, TString plottodraw, TString binNames, RooAbsData* data, RooFitResult* result){
  //********************************************************************

  std::vector<TCanvas*> rooplots;

  if(!work){
    std::cerr << "ERROR::NULL workspace. Will return empty vector!" << std::endl;
    return rooplots;
  }

  // Silence output
  RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);

  // Get pdf from workspace
  RooSimultaneous* pdf = (RooSimultaneous*)work->pdf("simPdf");
  if(!pdf){
    std::cerr << "ERROR::No pdf found in workspace. Will return empty vector!" << std::endl;
    return rooplots;
  }

  // If not provide, get data from workspace
  if(!data)
    data = work->data("obsData");

  // Get category components
  RooCategory* categories = work->cat("channelCat");

  // Print table with the number of data entries per category (channel)
  data->table(*((RooAbsCategory*)categories))->Print("v");

  std::vector<TString> categoriesName;
  TIterator* iter = categories->typeIterator();
  RooCatType* catType;
  while( (catType = (RooCatType*) iter->Next())) {
    TString catname = catType->GetName();
    categoriesName.push_back(catname);
  }

  for(UInt_t i = 0; i < categoriesName.size(); i++){
    TString catname = categoriesName[i];

    if(categories->setLabel(catname)){
      std::cout << "WARNING::Category " << catname.Data() << " is not a member of channelCat" << std::endl;
      continue;
    }

    RooAbsPdf* subpdf = (RooAbsPdf*)pdf->getPdf(catname.Data());
    if(!subpdf){
      std::cout << "WARNING::Can't find sub-pdf for region " << catname.Data() << std::endl;
      continue;
    }

    TString subdataset_str = Form("channelCat==channelCat::%s",catname.Data());
    RooAbsData* subdataset = (RooAbsData*) data->reduce(subdataset_str.Data());
    if(!subdataset){
      std::cout << "WARNING::Can't find sub-dataset for region " << catname.Data() << std::endl;
      continue;
    }

    RooRealVar* var =(RooRealVar*) ((RooArgSet*) subpdf->getObservables(*subdataset))->find(Form("obs_x_%s", catname.Data()));

    // Define canvas first
    TString canName = Form("canvas_%s_%s", catname.Data(), name.Data());
    TCanvas* c = new TCanvas(canName, canName);

    // Legend
    TLegend* legend = new TLegend(0.10,0.70,0.85,0.90,"");
    legend->SetFillStyle(0);
    legend->SetFillColor(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.036);

    legend->SetNColumns(2);

    TLegendEntry* entry = legend->AddEntry("","Data","pl") ;
    entry->SetMarkerColor(kBlack);
    entry->SetMarkerStyle(20);
    entry=legend->AddEntry("","MC Total","lf") ;
    entry->SetLineColor(kBlack);
    entry->SetFillColor(kBlue);
    entry->SetFillStyle(3444);

    RooPlot* frame = var->frame();
    frame->SetName(Form("frame_%s_%s", catname.Data(), name.Data()));

    // Draw data and pdf on top of error band - change to RooAbsData::SumW2 if data is weighted
    if(error.Contains("SumW2") || error.Contains("sumw2") || error.Contains("sumW2"))
      subdataset->plotOn(frame, RooFit::DataError(RooAbsData::SumW2), RooFit::MarkerColor(1), RooFit::LineColor(1));
    else
      subdataset->plotOn(frame, RooFit::DataError(RooAbsData::Poisson), RooFit::MarkerColor(1), RooFit::LineColor(1));

    TString RRSumPdfName = Form("%s_model",catname.Data());
    RooRealSumPdf* RRSumPdf = (RooRealSumPdf*) subpdf->getComponents()->find(RRSumPdfName);

    RooArgList RRSumComponentsList =  RRSumPdf->funcList();

    TString binWidthName =  Form("binWidth_obs_x_%s_0",catname.Data());
    RooRealVar* regionBinWidth = ((RooRealVar*) RRSumPdf->getVariables()->find(Form("binWidth_obs_x_%s_0",catname.Data()))) ;
    if(regionBinWidth == NULL)
      std::cout << "WARNING::bindWidth variable not found for region(" << catname << "). PLOTTING COMPONENTS WILL BE WRONG!" << std::endl;

    // normalize data to expected number of events
    double normCount = subpdf->expectedEvents(*var);
    if(result)
      std::cout << "INFO::MC events after fit  = " << normCount << " for " << subpdf->GetName() << std::endl;
    else 
      std::cout << "INFO::MC events before fit = " << normCount << " for " << subpdf->GetName() << std::endl;

    RooLinkedListIter iter = RRSumComponentsList.iterator() ;
    RooProduct* component;
    std::vector<TString> compNameVec;
    std::vector <double> compFracVec;
    std::vector<TString> compStackNameVec;
    std::vector <double> compStackFracVec;
    compNameVec.clear();
    compStackNameVec.clear();
    compFracVec.clear();
    compStackFracVec.clear();

    while( (component = (RooProduct*) iter.Next()) ){
      TString  componentName = component->GetName();
      TString stackComponentName = componentName;
      if(!compStackNameVec.empty())
        stackComponentName  = Form("%s,%s",compStackNameVec.back().Data() ,componentName.Data());
      compNameVec.push_back(componentName);
      compStackNameVec.push_back(stackComponentName);

      RooAbsReal*  i_RRSumPdf = ((RooAbsPdf*)work->pdf(RRSumPdfName))->createIntegral(RooArgSet(*var));
      Double_t Int_RRSumPdf = i_RRSumPdf->getVal();
      RooAbsReal*  i_component =   ((RooProduct*)work->obj(componentName))->createIntegral(RooArgSet(*var));
      Double_t Int_component = i_component->getVal();

      Double_t componentFrac = 0.;
      if(Int_RRSumPdf != 0.)
        componentFrac =  Int_component * regionBinWidth->getVal() / Int_RRSumPdf;

      double stackComponentFrac = componentFrac;
      if(!compStackFracVec.empty())
        stackComponentFrac  = compStackFracVec.back() + componentFrac;

      compFracVec.push_back(componentFrac);
      compStackFracVec.push_back(stackComponentFrac);
    }

    for(Int_t iVec = (compFracVec.size()-1); iVec>-1; iVec--){
      Int_t  compPlotColor = iVec;
      if(compPlotColor <= 0) compPlotColor = kMagenta;
      if(compPlotColor == 1) compPlotColor = kGray;
      if(compPlotColor == 6) compPlotColor = 47;
      if(compPlotColor == 10) compPlotColor = 20;
      subpdf->plotOn(frame,RooFit::Components(compStackNameVec[iVec].Data()),RooFit::FillColor(compPlotColor),RooFit::FillStyle(3001),RooFit::DrawOption("F"),RooFit::Normalization(compStackFracVec[iVec]*normCount,RooAbsReal::NumEvent),RooFit::Precision(1e-5));
    }

    if(result)
      subpdf->plotOn(frame, RooFit::Normalization(1,RooAbsReal::RelativeExpected), RooFit::Precision(1e-5), RooFit::VisualizeError(*result), RooFit::FillColor(kBlue), RooFit::FillStyle(3444));

    // Plot again so that it is on top of errors
    subpdf->plotOn(frame, RooFit::Normalization(1,RooAbsReal::RelativeExpected), RooFit::Precision(1e-5), RooFit::LineColor(1));
    if(error.Contains("SumW2") || error.Contains("sumw2") || error.Contains("sumW2"))
      subdataset->plotOn(frame, RooFit::DataError(RooAbsData::SumW2), RooFit::MarkerColor(1), RooFit::LineColor(1));
    else
      subdataset->plotOn(frame, RooFit::DataError(RooAbsData::Poisson), RooFit::MarkerColor(1), RooFit::LineColor(1));

    // Remove empty data bins
    FitUtils::RemoveEmptyBins(frame);

    // These are hardcoded - need to change
    //TString compshnamevec[8] = {"FGD1 #nu-e Elastic", "FGD1 CC-#nu_{e}", "FGD1 Other Backgrounds", "FGD1 LDM Signal", "FGD2 #nu-e Elastic", "FGD2 CC-#nu_{e}", "FGD2 Other Backgrounds", "FGD2 LDM Signal"};

    for(int iComp = (compNameVec.size()-1) ; iComp>-1; iComp--){
      Int_t compPlotColor    = iComp;
      if(compPlotColor <= 0) compPlotColor = kMagenta;
      if(compPlotColor == 1) compPlotColor = kGray;
      if(compPlotColor == 6) compPlotColor = 47;
      if(compPlotColor == 8) compPlotColor = 20;
      //TString  compShortName  = compshnamevec[iComp];
      //TString legName = compShortName;
      //entry=legend->AddEntry("",legName.Data(),"f");
      if     (compNameVec[iComp].Contains("LDM_FGD1")) entry=legend->AddEntry("","FGD1 LDM Signal","f");
      else if(compNameVec[iComp].Contains("LDM_FGD2")) entry=legend->AddEntry("","FGD2 LDM Signal","f");
      else if(compNameVec[iComp].Contains("LDM_FGD_")) entry=legend->AddEntry("","FGD LDM Signal","f");
      else if(compNameVec[iComp].Contains("_LDM_")) entry=legend->AddEntry("","LDM Signal","f");
      else if(compNameVec[iComp].Contains("NuElecElastic_FGD1")) entry=legend->AddEntry("","FGD1 #nu-e Elastic","f");
      else if(compNameVec[iComp].Contains("NuElecElastic_FGD2")) entry=legend->AddEntry("","FGD2 #nu-e Elastic","f");
      else if(compNameVec[iComp].Contains("NuElecElastic_FGD_")) entry=legend->AddEntry("","FGD #nu-e Elastic","f");
      else if(compNameVec[iComp].Contains("NuElecElastic")) entry=legend->AddEntry("","#nu-e Elastic","f");
      else if(compNameVec[iComp].Contains("CCnue_FGD1")) entry=legend->AddEntry("","FGD1 CC-#nu_{e}","f");
      else if(compNameVec[iComp].Contains("CCnue_FGD2")) entry=legend->AddEntry("","FGD2 CC-#nu_{e}","f");
      else if(compNameVec[iComp].Contains("CCnue_FGD_")) entry=legend->AddEntry("","FGD CC-#nu_{e}","f");
      else if(compNameVec[iComp].Contains("CCnue")) entry=legend->AddEntry("","CC-#nu_{e}","f");
      else if(compNameVec[iComp].Contains("Otherbkg_FGD1")) entry=legend->AddEntry("","FGD1 Other Backgrounds","f");
      else if(compNameVec[iComp].Contains("Otherbkg_FGD2")) entry=legend->AddEntry("","FGD2 Other Backgrounds","f");
      else if(compNameVec[iComp].Contains("Otherbkg_FGD_")) entry=legend->AddEntry("","FGD Other Backgrounds","f");
      else if(compNameVec[iComp].Contains("Otherbkg")) entry=legend->AddEntry("","Other Backgrounds","f");
      else entry=legend->AddEntry("","Unknown","f");
      entry->SetLineColor(compPlotColor);
      entry->SetFillColor(compPlotColor);
      entry->SetFillStyle(3001);
    }

    // two pads, one for 'standard' plot, one for data/MC ratio
    float yMinP1=0.305;
    float bottomMarginP1=0.005;
    TPad *pad1 = new TPad(Form("%s_pad1",canName.Data()),Form("%s_pad1",canName.Data()),0.,yMinP1,.99,1);
    pad1->SetBottomMargin(bottomMarginP1);
    pad1->SetFillColor(kWhite);
    pad1->SetTickx();
    pad1->SetTicky();
    TPad *pad2 = new TPad(Form("%s_pad2",canName.Data()),Form("%s_pad2",canName.Data()),0.,0.01,.99,0.295);
    pad2->SetTopMargin(0.005);
    pad2->SetBottomMargin(0.3);
    pad2->SetFillColor(kWhite);

    pad1->Draw();
    pad2->Draw();
    frame->GetXaxis()->SetLabelSize(0.);

    pad1->cd();
    frame->SetTitle("");
    frame->Draw();//Draw("same");
    frame->GetXaxis()->SetTitle("Reco Momentum Bin");
    frame->GetYaxis()->SetTitle("Events");
    legend->Draw();

    pad2->cd();
    RooPlot* frame_dummy =  var->frame();
    subdataset->plotOn(frame_dummy,RooFit::Cut(subdataset_str),RooFit::DataError(RooAbsData::Poisson));

    // Get RooHist of the data - will be used later in the ratio plot
    const char* curvename = 0;
    RooHist* nominal_hist = (RooHist*) frame_dummy->findObject(curvename, RooHist::Class());
    if(!nominal_hist){
      std::cerr << "ERROR::Drawing the data/MC histogram. Can't find nominal histogram." << std::endl;
      return rooplots;
    }

    // normalize pdf to number of expected events, not to number of events in dataset
    subpdf->plotOn(frame_dummy,RooFit::Normalization(1,RooAbsReal::RelativeExpected),RooFit::Precision(1e-5));

    // Get RooCurve of the pdf - will be used later in the ratio plot
    RooCurve* nominal_curve = (RooCurve*) frame_dummy->findObject(curvename, RooCurve::Class());
    if(!nominal_curve){
      std::cerr << "ERROR::Drawing the data/MC histogram. Can't find nominal curve." << std::endl;
      return rooplots;
    }

    RooHist* hratio = NULL;
    RooCurve* ratio_curve = new RooCurve;

    if(plottodraw == "pull"){
      hratio = (RooHist*) frame_dummy->pullHist();
      hratio->SetTitle("Pull Distribution");
    }
    else if(plottodraw == "residuals"){
      hratio = (RooHist*) frame_dummy->residHist();
      hratio->SetTitle("Residual Distribution");
    }
    else if(plottodraw == "ratio"){
      if(result)
        subpdf->plotOn(frame_dummy, RooFit::Normalization(1,RooAbsReal::RelativeExpected), RooFit::Precision(1e-5), RooFit::VisualizeError(*result, 1), RooFit::FillColor(kBlue-5), RooFit::FillStyle(3004));

      // Get error plot
      RooCurve* error_curve = (RooCurve*)frame_dummy->findObject(curvename, RooCurve::Class());
      if(!error_curve){
	std::cerr << "ERROR::Drawing the data/MC histogram. Can't find error curve." << std::endl;
        return rooplots;
      }

      ratio_curve->SetName(Form("%s_ratiocurve",nominal_curve->GetName()));
      ratio_curve->SetLineColor(kBlue-5);
      ratio_curve->SetFillColor(kBlue-5);
      ratio_curve->SetFillStyle(3004);
      ratio_curve->SetLineWidth(1);

      // Fill error curve
      Int_t j = 0;
      bool bottomCurve = false;
      for(Int_t i=1; i < error_curve->GetN()-1; i++){
        Double_t x = 0.;
        Double_t y = 0.;
        error_curve->GetPoint(i,x,y) ;

        if( i >= (nominal_curve->GetN()-1) ) bottomCurve = true;

        Double_t xNom = x;
        Double_t yNom = y;

        if( i == (nominal_curve->GetN() - 1) ||  i == nominal_curve->GetN() ){
          ratio_curve->addPoint(x, 0.);
          continue;
        }

	if(bottomCurve){
          nominal_curve->GetPoint(j,xNom,yNom);
          j--;
        }
        else{
          j++;
          nominal_curve->GetPoint(j,xNom,yNom);
        }

        if(fabs(yNom) > 0.00001){
          ratio_curve->addPoint(x, (y / yNom));
        }
        else{
          ratio_curve->addPoint(x, 0.);
        }
      }

      // Define ratio plot
      hratio = new RooHist(nominal_hist->getNominalBinWidth());

      // Determine range of curve
      Double_t xstart, xstop, y;
      nominal_curve->GetPoint(2,xstart,y);
      nominal_curve->GetPoint(nominal_curve->GetN()-1,xstop,y);

      for(Int_t i=0; i < nominal_hist->GetN(); i++){
        Double_t x,ypoint;
        nominal_hist->GetPoint(i,x,ypoint);

        if(x < xstart || x > xstop) continue;
        if(fabs(ypoint) < 0.000001 ) continue;

        Double_t yy;
        Double_t yerrorl = nominal_hist->GetErrorYlow(i);
        Double_t yerrorh = nominal_hist->GetErrorYhigh(i);

	yy = ypoint / nominal_curve->interpolate(x);
        yerrorl /= nominal_curve->interpolate(x);
        yerrorh /= nominal_curve->interpolate(x);

        hratio->addBinWithError(x,yy,yerrorl,yerrorh);
      }
    }

    if(!hratio){
      std::cerr << "ERROR::Drawing data and MC histogram failed. RooHist is not found." << std::endl;
      return rooplots;
    }

    hratio->SetMarkerColor(kRed);
    hratio->SetLineColor(kRed);

    // Create a new frame to draw the residual distribution and add the distribution to the frame
    RooPlot* frame2 = var->frame();
    if(plottodraw == "ratio" && result) frame2->addPlotable(ratio_curve,"F");
    frame2->addPlotable(hratio,"P");

    // These are hardcoded - need to change
    std::string sbinNames(binNames.Data());
    std::vector<std::string> recoBinName = FitUtils::SplitString(sbinNames, ",");
    for(unsigned int i = 0; i < recoBinName.size(); i++)
      frame2->GetXaxis()->SetBinLabel(i+1, recoBinName[i].c_str());
    
    frame2->GetXaxis()->SetTitle("p [MeV/c]");

    // Cosmetics
    int firstbin = frame_dummy->GetXaxis()->GetFirst();
    int lastbin  = frame_dummy->GetXaxis()->GetLast();
    double xmax  = frame_dummy->GetXaxis()->GetBinUpEdge(lastbin);
    double xmin  = frame_dummy->GetXaxis()->GetBinLowEdge(firstbin);

    if(plottodraw == "pull"){
      TLine* lp1 = new TLine(xmin,1.,xmax,1.);
      TLine* lp2 = new TLine(xmin,2.,xmax,2.);
      TLine* lp3 = new TLine(xmin,3.,xmax,3.);
      TLine* lp4 = new TLine(xmin,4.,xmax,4.);

      TLine* lp6 = new TLine(xmin,-1.,xmax,-1.);
      TLine* lp7 = new TLine(xmin,-2.,xmax,-2.);
      TLine* lp8 = new TLine(xmin,-3.,xmax,-3.);
      TLine* lp9 = new TLine(xmin,-4.,xmax,-4.);

      TLine* lp10 = new TLine(xmin,0,xmax,0);

      lp1->SetLineStyle(3);
      lp2->SetLineStyle(3);
      lp3->SetLineStyle(3);
      lp4->SetLineStyle(3);
      lp6->SetLineStyle(3);
      lp7->SetLineStyle(3);
      lp8->SetLineStyle(3);
      lp9->SetLineStyle(3);
      lp10->SetLineStyle(3);

      frame2->addObject(lp1);
      frame2->addObject(lp2);
      frame2->addObject(lp3);
      frame2->addObject(lp4);
      frame2->addObject(lp6);
      frame2->addObject(lp7);
      frame2->addObject(lp8);
      frame2->addObject(lp9);
      frame2->addObject(lp10);

      frame2->SetMinimum(-4.0);
      frame2->SetMaximum(4.0);

      frame2->GetYaxis()->SetTitle("Pull");
    }
    else if(plottodraw == "residuals"){
      TLine* l = new TLine(xmin,0.,xmax,0.);
      l->SetLineWidth(1);
      l->SetLineStyle(2);
      frame2->addObject(l);
      frame2->GetYaxis()->SetTitle("Residual");
    }
    else if(plottodraw == "ratio"){
      TLine* lp1 = new TLine(xmin,1.,xmax,1.);
      TLine* lp2 = new TLine(xmin,0.5,xmax,0.5);
      TLine* lp3 = new TLine(xmin,1.5,xmax,1.5);
      TLine* lp4 = new TLine(xmin,2.,xmax,2.);
      TLine* lp5 = new TLine(xmin,2.5,xmax,2.5);
      lp1->SetLineWidth(1);
      lp1->SetLineStyle(3);
      lp2->SetLineStyle(3);
      lp3->SetLineStyle(3);
      lp4->SetLineStyle(3);
      lp5->SetLineStyle(3);
      frame2->addObject(lp1);
      frame2->addObject(lp2);
      frame2->addObject(lp3);
      frame2->addObject(lp4);
      frame2->addObject(lp5);

      frame2->SetMinimum(.0);
      frame2->SetMaximum(3.0);

      frame2->GetYaxis()->SetTitle("Data / MC");
    }

    frame2->GetYaxis()->SetLabelSize(0.10);
    frame2->GetYaxis()->SetNdivisions(504);
    frame2->GetXaxis()->SetLabelSize(0.10);
    frame2->GetYaxis()->SetTitleSize(0.10);
    frame2->GetXaxis()->SetTitleSize(0.10);
    frame2->GetYaxis()->SetTitleOffset(0.35);
    frame2->GetXaxis()->SetTitleOffset(1.);
    frame2->GetYaxis()->SetLabelOffset(0.01);
    frame2->GetXaxis()->SetLabelOffset(0.03);
    frame2->GetXaxis()->SetTickLength(0.06);

    frame2->SetTitle("");
    frame2->GetYaxis()->CenterTitle();
    frame2->Draw();

    rooplots.push_back(c);
  }

  return rooplots;

}

//********************************************************************
void FitUtils::RemoveEmptyBins(RooPlot* frame){
  //********************************************************************

  const char* histname = 0;

  // Find histogram
  RooHist* histo = (RooHist*) frame->findObject(histname, RooHist::Class());
  if(!histo) return;

  for(Int_t i=0; i < histo->GetN(); i++){
    Double_t x,y;
    histo->GetPoint(i,x,y);

    if(fabs(y) < 0.00001 && histo->GetErrorYhigh(i) > 0.){
      histo->RemovePoint(i);
      if(i != histo->GetN()) --i;
    }
  }

  return;

}

//********************************************************************
std::vector<TCanvas*> FitUtils::PlotNLL(RooWorkspace *work, TString name, RooFitResult* result, bool plotPLL){
  //********************************************************************

  std::vector<TCanvas*> rooplots;

  if(!work){
    std::cerr << "ERROR::NULL workspace. Will return empty vector!" << std::endl;
    return rooplots;
  }

  if(!result){
    std::cerr << "ERROR::NULL fit result. Will return empty vector!" << std::endl;
    return rooplots;
  }

  RooSimultaneous* pdf = (RooSimultaneous*) work->pdf("simPdf");
  if(!pdf){
    std::cout << "ERROR::No pdf found in workspace. Will return empty vector!" << std::endl;
    return rooplots;
  }

  // Get data from workspace
  RooAbsData* data = work->data("obsData");

  // Get category components
  RooCategory* categories = work->cat("channelCat");
  // Print table with the number of data entries per category (channel)
  data->table(*((RooAbsCategory*)categories))->Print("v");

  RooStats::ModelConfig *combined_config = (RooStats::ModelConfig*)work->obj("ModelConfig");
  if(!combined_config){
    std::cout << "ERROR::No model config " << " ModelConfig " << " in workspace.  Will return empty vector!" << std::endl;
    return rooplots;
  }

  const RooArgSet* obs = combined_config->GetGlobalObservables();
  if(!obs){
    std::cout << "ERROR::No observables found in ModelConfig.  Will return empty vector!" << std::endl;
    return rooplots;
  }

  RooArgList floatParsFinal = result->floatParsFinal();

  // Create NLL
  RooAbsReal* nll = pdf->createNLL(*data, RooFit::NumCPU(2), RooFit::GlobalObservables(*obs), RooFit::Offset(true));

  for(Int_t i = 0; i < floatParsFinal.getSize(); i++){
    RooAbsArg* arg = floatParsFinal.at(i);
    if(!arg->InheritsFrom("RooRealVar")) continue;

    RooRealVar* par = (RooRealVar*) arg;
    TString parName = par->GetName();

    // set parameter range to readable range
    double minRange = par->getMin();
    double maxRange = par->getMax();
    if(minRange < 0.){
      par->setMin(-3.);
      par->setMax(3.);
    }
    else {
      par->setMin(minRange);
      par->setMax(2.);
    }

    RooPlot* frame = par->frame();
    nll->plotOn(frame, RooFit::ShiftToZero());
    frame->SetMinimum(0.);
    // To be able to see the 1/2 sigma
    frame->SetMaximum(2.5);

    RooAbsReal* pll = NULL;
    if(plotPLL) {
      pll = nll->createProfile(*par) ;
      pll->plotOn(frame, RooFit::LineColor(kRed), RooFit::LineStyle(kDashed), RooFit::NumCPU(4));
    }

    TString canName=Form("Canvas_NLL_%s_%s", name.Data(), parName.Data());
    TCanvas* c = new TCanvas(canName,canName,600,600);
    c->cd();
    frame->Draw();

    TLegend* legend = new TLegend(0.55,0.65,0.85,0.95,"");
    legend->SetFillStyle(0);
    legend->SetFillColor(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.036);
    TLegendEntry* entry=legend->AddEntry("","NLL","l") ;
    entry->SetLineColor(kBlue);
    if(plotPLL){
      entry=legend->AddEntry("","PLL","l") ;
      entry->SetLineColor(kRed);
      entry->SetLineStyle(kDashed);
    }
    legend->Draw();

    // reset parameter range to previous values
    par->setMin(minRange);
    par->setMax(maxRange);

    if(pll) delete pll;

    rooplots.push_back(c);
  }

  return rooplots;

}

//********************************************************************
TCanvas* FitUtils::PlotNuisanceParameters(TTree* tree, RooWorkspace* ws){
  //********************************************************************

  if(!tree || !ws){
    std::cerr << "ERROR::No tree or workspace found. Pull plot failed!" << std::endl;
    return NULL;
  }

  RooStats::ModelConfig *combined_config = (RooStats::ModelConfig*)ws->obj("ModelConfig");
  if(!combined_config){
    std::cerr << "ERROR::No model config " << " ModelConfig " << " in workspace. No nuisance parameters plots!" << std::endl;
    return NULL;
  }

  Int_t status;
  tree->SetBranchAddress("status", &status);
  
  RooRealVar* var(0);
  std::vector<Float_t> varValsNuisance;
  const RooArgSet* nuisanceParsList = combined_config->GetNuisanceParameters();
  if(!nuisanceParsList){
    std::cerr << "ERROR::No nuisance parameters found. Nuisance parameters plot failed!" << std::endl;
    return NULL;
  }
  const Int_t n = nuisanceParsList->getSize();
  varValsNuisance.resize( nuisanceParsList->getSize(), -999. );

  TH1F* nuisancehisto  = new TH1F("nuisanceplot",  "", n+1, 0, n+1);

  TIterator* Itr = nuisanceParsList->createIterator();
  for (Int_t i=0; (var = (RooRealVar*)Itr->Next()); ++i) {
    TString varName = var->GetName();
    if(!varName.Contains("alpha") && !varName.Contains("gamma_stat")) continue;
    if(varName.Contains("Lumi") || varName.Contains("nom") || varName.Contains("binWidth") || varName.Contains("err") || varName.Contains("pull")) continue;

    tree->SetBranchAddress(varName.Data(), &varValsNuisance[i]);

    //std::cout << "INFO::Plotting nuisance parameter " << varName << std::endl;
    TH1F* temphisto  = new TH1F("temphisto",  "temphisto", 100, -5, 5);

    for(Int_t j=0; j < tree->GetEntries(); j++){
      tree->GetEntry(j);
      if(status != 0) continue;
      temphisto->Fill(varValsNuisance[i]);
      //nuisancehisto->SetBinContent(i+1, varValsNuisance[i]);
    }
    nuisancehisto->SetBinContent(i+1, temphisto->GetMean());
    delete temphisto;
  }

  delete Itr;

  Int_t ngammas = 0;
  // Fill the error
  TIterator* Itr2 = nuisanceParsList->createIterator();
  for (Int_t i=0; (var = (RooRealVar*)Itr2->Next()); ++i) {
    TString varName = var->GetName() + TString("err");

    if(!varName.Contains("alpha") && !varName.Contains("gamma_stat")) continue;
    if(!varName.Contains("err")) continue;
    if(varName.Contains("Lumi") || varName.Contains("nom") || varName.Contains("binWidth") || varName.Contains("pull")) continue;

    tree->SetBranchAddress(varName.Data(), &varValsNuisance[i]);

    if(varName.Contains("gamma_stat"))
      ngammas++;

    TH1F* temphisto  = new TH1F("temphisto",  "temphisto", 100, -5, 5);

    for(Int_t j=0; j < tree->GetEntries(); j++){
      tree->GetEntry(j);
      if(status != 0) continue;
      temphisto->Fill(varValsNuisance[i]);
      //nuisancehisto->SetBinError(i+1, varValsNuisance[i]);
    }

    nuisancehisto->SetBinError(i+1, temphisto->GetMean());
    delete temphisto;
  }

  delete Itr2;

  nuisancehisto->SetLineColor(1);
  nuisancehisto->SetMarkerStyle(21);
  nuisancehisto->SetMarkerColor(1);
  nuisancehisto->GetYaxis()->SetRangeUser(-2.0,2.0);
  nuisancehisto->GetXaxis()->SetTitle("Systematic ID");
  nuisancehisto->GetYaxis()->SetTitle("Fit Result");
  if(tree->GetEntries() > 1)
    nuisancehisto->GetYaxis()->SetTitle("Mean Fit Result From Toys");
  nuisancehisto->SetTitle("Nuisance Parameters");
  nuisancehisto->SetStats(false);

  TLine *mline = new TLine(0,0.0,n,0.0);
  mline->SetLineColor(kRed);
  TLine *sline = new TLine(0,1.0,n,1.0);
  sline->SetLineColor(kRed);
  TLine *s1line = new TLine(0,-1.0,n,-1.0);
  s1line->SetLineColor(kRed);
  TLine *vline = new TLine(n-ngammas,-2,n-ngammas,2.0);
  vline->SetLineColor(kGreen);

  TCanvas* cNuisanceParameters = new TCanvas("cNuisanceParameters","cNuisanceParamters");
  nuisancehisto->Draw("e");
  mline->Draw("same");
  sline->Draw("same");
  s1line->Draw("same");
  vline->Draw("same");

  return cNuisanceParameters;

}

