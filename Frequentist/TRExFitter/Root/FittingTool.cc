// Class include
#include "TRExFitter/FittingTool.h"

//Framework includes
#include "TRExFitter/Common.h"
#include "TRExFitter/FitUtils.h"
#include "TRExFitter/StatusLogbook.h"
#include "TRExFitter/YamlConverter.h"

//ROOR includes
#include "TCanvas.h"
#include "TH2.h"
#include "TRandom3.h"
#include "TFile.h"

//Roostats includes
#include "Math/MinimizerOptions.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/ModelConfig.h"

//Roofit includes
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooMinimizer.h"
#include "RooFitResult.h"
#include "RooArgSet.h"

//c++ includes
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

//________________________________________________________________________
//
FittingTool::FittingTool():
    m_CPU(1),
    m_useMinos(false),
    m_constPOI(false),
    m_fitResult(nullptr),
    m_noGammas(false),
    m_noSystematics(false),
    m_noNormFactors(false),
    m_noShapeFactors(false),
    m_RangePOI_up(100.),
    m_RangePOI_down(-10.),
    m_randomize(false),
    m_randomNP(0.1),
    m_randSeed(-999),
    m_minNll(999999),
    m_externalConstraints(nullptr),
    m_strategy(-1),
    m_useHesse(true),
    m_hesseBeforeMigrad(false)
{
}

//________________________________________________________________________
//
FittingTool::~FittingTool() {
}

//________________________________________________________________________
//
void FittingTool::SetSubCategories() {
    WriteDebugStatus("FittingTool::SetSubCategories", "finding unique SubCategories");
    // loop over m_subCategoryMap to find all unique SubCategories, save in m_subCategories set
    for(std::map<std::string, std::string>::iterator it = m_subCategoryMap.begin(); it != m_subCategoryMap.end(); ++it) {
        m_subCategories.insert(it->second);
    }
}

//________________________________________________________________________
//
void FittingTool::AddValPOI(const std::string& name, const double value) {
    auto it = std::find_if(m_valPOIs.begin(), m_valPOIs.end(),
    [&name](const std::pair<std::string, double>& element){ return element.first == name;});

    // take only unique ones
    if (it != m_valPOIs.end()) return;

    m_valPOIs.emplace_back(std::make_pair(name, value));
}

//____________________________________________________________________________________
//
void FittingTool::ReplacePOIVal(const std::string& name, const double value) {
    auto it = std::find_if(m_valPOIs.begin(), m_valPOIs.end(), 
        [&](const std::pair<std::string, double>& element){return name == element.first;});

    if (it == m_valPOIs.end()) return;

    it->second = value;
}

//________________________________________________________________________
//
double FittingTool::FitPDF( RooStats::ModelConfig* model, RooAbsPdf* fitpdf, RooAbsData* fitdata, bool fastFit, bool noFit, bool saturatedModel ) {

    WriteDebugStatus("FittingTool::FitPDF", "-> Entering in FitPDF function");

    //
    // Printing the whole model for information
    //
    if(TRExFitter::DEBUGLEVEL >= 2) model->Print();

    //
    // Getting the list of model that can be constrained (parameters of the MC)
    //
    RooArgSet* constrainedParams = fitpdf->getParameters(*fitdata);
    RooStats::RemoveConstantParameters(constrainedParams);
    RooFit::Constrain(*constrainedParams);

    //
    // Get the global observables (nominal values)
    //
    const RooArgSet* glbObs = model->GetGlobalObservables();

    //
    // Create the likelihood based on fitpdf, fitData and the parameters
    //
    std::unique_ptr<RooAbsReal> nll(fitpdf->createNLL(*fitdata,
                                    RooFit::Constrain(*constrainedParams),
                                    RooFit::GlobalObservables(*glbObs),
                                    RooFit::Offset(1),
                                    RooFit::NumCPU(m_CPU,RooFit::Hybrid),
                                    RooFit::Optimize(kTRUE),
                                    RooFit::ExternalConstraints(*m_externalConstraints)
                                    ));

    //
    // Needed for Ranking plot, but also to set random initial values for the NPs
    //
    if(m_randSeed == -999){
        gRandom->SetSeed(time(nullptr));
    }
    else{
        gRandom->SetSeed(m_randSeed);
    }

    //
    // Getting the POIs
    //
    std::vector<RooRealVar*> pois = this->GetVectorPOI(model);

    if (pois.empty()) {
        WriteErrorStatus("FittingTool::FitPDF", "No POIs found!");
        return 0;
    }

    for (auto poi : pois) {
        poi->setConstant(m_constPOI || saturatedModel);

        const std::string name = poi->GetName();
        double value(0);
        auto it = std::find_if(m_valPOIs.begin(), m_valPOIs.end(),
            [&name](const std::pair<std::string, double>& element){ return element.first == name;});

        if (it != m_valPOIs.end()) {
            value = it->second;
        }
        poi->setVal(value);

        // randomize the POI
        if(!m_constPOI && m_randomize){
            poi->setVal( value + m_randomNP*(gRandom->Uniform(2)-1.) );
        }

        WriteInfoStatus("FittingTool::FitPDF", "POI: " + name);
        WriteInfoStatus("FittingTool::FitPDF", "   -> Constant POI : " + std::to_string(poi->isConstant()));
        WriteInfoStatus("FittingTool::FitPDF", "   -> Initial value of POI : " + std::to_string(poi->getVal()));
    }

    if (saturatedModel) {
        FitUtils::FixAllParameters(model, {"saturated_model"});
    } else if (model->GetNuisanceParameters()) {
        for(auto var_tmp : *model->GetNuisanceParameters()) {
            RooRealVar* var = static_cast<RooRealVar*>(var_tmp);
            const std::string np = var->GetName();
            bool found = false;
            //
            // first check if all systs, norm and gammas should be set to constant
            if((np.find("gamma_stat")!=string::npos || np.find("gamma_shape_stat")!=string::npos) && m_noGammas){
                WriteDebugStatus("FittingTool::FitPDF", "setting to constant : " + np + " at value " + std::to_string(var->getVal()));
                var->setConstant( 1 );
                var->setVal( 1 );
                found = true;
            }
            else if((np.find("alpha_")!=string::npos || (np.find("gamma_shape")!=string::npos && np.find("gamma_shape_stat")==string::npos)) && m_noSystematics){
                WriteDebugStatus("FittingTool::FitPDF", "setting to constant : " + np + " at value " + std::to_string(var->getVal()));
                var->setConstant( 1 );
                var->setVal( 0 );
                found = true;
            }
            else if(np.find("alpha_")==string::npos && np.find("gamma_")==string::npos && (m_noNormFactors || saturatedModel)){
                WriteDebugStatus("FittingTool::FitPDF", "setting to constant : " + np + " at value " + std::to_string(var->getVal()));
                var->setConstant( 1 );
                found = true;
            }
            if(found) continue;
            //
            // set to constant the saturatedModel shape factor parameters if saturatedModel is left to false
            if(np.find("saturated_model_sf_")!=std::string::npos && !saturatedModel){
                WriteDebugStatus("FittingTool::FitPDF", "setting to constant : " + np + " at value " + std::to_string(var->getVal()));
                var->setConstant( 1 );
                found = true;
            }
            if(found) continue;
            //
            // loop on the NP specified to have custom starting value
            for( unsigned int i_np = 0; i_np<m_initialNP.size(); i_np++ ){
                if( np == ("alpha_"+m_initialNP[i_np]) || np == m_initialNP[i_np] ){
                    var->setVal(m_initialNPvalue[i_np]);
                    WriteInfoStatus("FittingTool::FitPDF", " ---> Setting " + m_initialNP[i_np] + " to "  +std::to_string(m_initialNPvalue[i_np]));
                    found = true;
                    break;
                }
            }
            //
            // loop on the NP specified to be constant - This should be after setting to initial value
            for( unsigned int i_np = 0; i_np<m_constNP.size(); i_np++ ){
                if( np == ("alpha_"+m_constNP[i_np]) || np == m_constNP[i_np]
                    || np == ("gamma_"+m_constNP[i_np])
                ){
                    WriteInfoStatus("FittingTool::FitPDF", "setting to constant : " + np + " at value " + std::to_string(m_constNPvalue[i_np]));
                    var->setVal(m_constNPvalue[i_np]);
                    var->setConstant(1);
                    found = true;
                    break;
                }
            }
            if(!found){
                if( np.find("alpha_")!=string::npos ){   // for syst NP
                    if(m_randomize) var->setVal( m_randomNP*(gRandom->Uniform(2)-1.) );
                    else            var->setVal(0);
                    var->setConstant(0);
                }
                else {  // for norm factors & gammas
                    if(m_randomize) var->setVal( 1 + m_randomNP*(gRandom->Uniform(2)-1.) );
                    else            var->setVal( 1 );
                }
            }
        }
    }
    
    double nllval = nll->getVal();

    WriteDebugStatus("FittingTool::FitPDF","   -> Initial value of the NLL = " +std::to_string(nllval));
    if(TRExFitter::DEBUGLEVEL >= 2) constrainedParams->Print("v");

    //
    // return here if specified not to perform the fit
    if(noFit) {
        m_minNll = nllval;
        return nllval;
    }
    
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(1);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);

    // Safe fit loop
    int strat = ::ROOT::Math::MinimizerOptions::DefaultStrategy();
    if(m_strategy >= 0) {
        strat = m_strategy;
        WriteInfoStatus("FittingTool::FitPDF", "Manually setting strategy to "+std::to_string(strat));
    }
    
    const TString minimType = ::ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str();
    const TString algorithm = ::ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str();
    const double tol =        ::ROOT::Math::MinimizerOptions::DefaultTolerance(); //AsymptoticCalculator enforces not less than 1 on this

    RooMinimizer minim(*nll);
    minim.optimizeConst(2);
    minim.setStrategy(strat);
    minim.setMinimizerType(minimType);
    minim.setPrintLevel(TRExFitter::DEBUGLEVEL - 1);
    minim.setEps(tol);
    // it doesn't make sense to try more than 2 additional strategies
    const int maxRetries = 3 - strat;

    // fast fit - e.g. for ranking
    if(fastFit){
        minim.setStrategy(0);  // to be the same as ttH comb
        minim.setPrintLevel(0);
    }

    TStopwatch sw;
    sw.Start();

    // always run one fit
    WriteInfoStatus("FittingTool::FitPDF", "");
    WriteInfoStatus("FittingTool::FitPDF", "");
    WriteInfoStatus("FittingTool::FitPDF", "");
    WriteInfoStatus("FittingTool::FitPDF", "Fit try no." + std::to_string(1));
    WriteInfoStatus("FittingTool::FitPDF", "======================");
    WriteInfoStatus("FittingTool::FitPDF", "");

    int status = minim.minimize(minimType.Data(),algorithm.Data());
    if (m_useHesse) {
        if (status == 0 || m_hesseBeforeMigrad) {
            minim.hesse();
        }
    }
    std::unique_ptr<RooFitResult> r(minim.save());
    double edm = r->edm();
    status = r->status();

    // check if the fit converged
    bool fitIsNotGood = (status > 1) || (edm > 0.0001);

    // if not, loop maximum 2 times with increased strategy
    int nrItr = 0;
    while (nrItr<maxRetries && fitIsNotGood){
        WriteWarningStatus("FittingTool::FitPDF", "");
        WriteWarningStatus("FittingTool::FitPDF", "   *******************************");
        WriteWarningStatus("FittingTool::FitPDF", "   * Increasing Minuit strategy (was " + std::to_string(strat) + ")");
        ++strat;
        WriteWarningStatus("FittingTool::FitPDF", "   * Fit failed with : ");
        WriteWarningStatus("FittingTool::FitPDF", "      - minuit status " + std::to_string(status % 100));
        WriteWarningStatus("FittingTool::FitPDF", "      - hess status " + std::to_string(status / 100));
        WriteWarningStatus("FittingTool::FitPDF", "      - Edm = " + std::to_string(edm));
        WriteWarningStatus("FittingTool::FitPDF", "   * Retrying with strategy " + std::to_string(strat));
        WriteWarningStatus("FittingTool::FitPDF", "   ********************************");
        WriteWarningStatus("FittingTool::FitPDF", "");
        PrintMinuitHelp();
        minim.setStrategy(strat);
        status = minim.minimize(minimType.Data(),algorithm.Data());
        if (m_useHesse) {
            if (status <= 1 || m_hesseBeforeMigrad) {
                minim.hesse();
            }
        }
        r.reset(minim.save());
        edm = r->edm();
        status = r->status();

        fitIsNotGood = (status > 1) || (edm > 0.0001);
        nrItr++;
    }

    // if the fit is not good even after retries print an error message
    if (fitIsNotGood) {
        WriteErrorStatus("FittingTool::FitPDF", "");
        WriteErrorStatus("FittingTool::FitPDF", "");
        WriteErrorStatus("FittingTool::FitPDF", "");
        WriteErrorStatus("FittingTool::FitPDF", "");
        WriteErrorStatus("FittingTool::FitPDF", "***********************************************************");
        WriteErrorStatus("FittingTool::FitPDF", "Fit failure unresolved with status " + std::to_string(status));
        WriteErrorStatus("FittingTool::FitPDF", "   Please investigate your workspace");
        WriteErrorStatus("FittingTool::FitPDF", "   Find a wall : you will need it to crash your head on it");
        WriteErrorStatus("FittingTool::FitPDF", "***********************************************************");
        WriteErrorStatus("FittingTool::FitPDF", "");
        WriteErrorStatus("FittingTool::FitPDF", "");
        WriteErrorStatus("FittingTool::FitPDF", "");
        PrintMinuitHelp();
        m_fitResult = nullptr;

        return 0;
    }

    if(m_useMinos && !saturatedModel){
        if (model->GetNuisanceParameters()) {
            std::unique_ptr<RooArgSet> SliceNPs(new RooArgSet( *(model->GetNuisanceParameters()) ));
            SliceNPs->add(*(model->GetParametersOfInterest()));
            WriteDebugStatus("FittingTool::FitPDF", "Size of variables for MINOS: " + std::to_string(m_varMinos.size()));

            if (m_varMinos.at(0)!="all"){
                for(auto var_tmp : *model->GetNuisanceParameters()) {
                    RooRealVar* var = static_cast<RooRealVar*>(var_tmp);
                    const TString vname=var->GetName();
                    bool isthere=false;
                    for (unsigned int m=0;m<m_varMinos.size();++m){
                        if(vname.Contains(m_varMinos.at(m))) {
                            isthere=true;
                            break;
                        }
                    }
                    if (!isthere) SliceNPs->remove(*var, true, true);
                }
                for(auto var_tmp : *model->GetParametersOfInterest()) {
                    RooRealVar* var = static_cast<RooRealVar*>(var_tmp);
                    const TString vname=var->GetName();
                    bool isthere=false;
                    for (unsigned int m=0;m<m_varMinos.size();++m){
                        if(vname.Contains(m_varMinos.at(m))) {
                            isthere=true;
                            break;
                        }
                    }
                    if (!isthere) SliceNPs->remove(*var, true, true);
                }
                minim.minos(*SliceNPs);
            }
            else {
                minim.minos();
            }
        }
    }//end useMinos

    r.reset(minim.save());
    WriteInfoStatus("FittingTool::FitPDF", "");
    WriteInfoStatus("FittingTool::FitPDF", "");
    WriteInfoStatus("FittingTool::FitPDF", "");
    WriteInfoStatus("FittingTool::FitPDF", "***********************************************************");
    WriteInfoStatus("FittingTool::FitPDF", "         FIT FINALIZED SUCCESSFULLY : ");
    WriteInfoStatus("FittingTool::FitPDF", "            - minuit status " + std::to_string(status % 100));
    WriteInfoStatus("FittingTool::FitPDF", "            - hess status " + std::to_string(status / 100));
    WriteInfoStatus("FittingTool::FitPDF", "            - Edm = " + std::to_string(edm));
    if (status == 1) WriteWarningStatus("FittingTool::FitPDF", "            Covariance was made pos-def!");
    WriteInfoStatus("FittingTool::FitPDF", "***********************************************************");
    if (TRExFitter::DEBUGLEVEL >= 2) sw.Print();
    WriteInfoStatus("FittingTool::FitPDF", "");
    WriteInfoStatus("FittingTool::FitPDF", "");
    WriteInfoStatus("FittingTool::FitPDF", "");

    if(r!=nullptr) m_fitResult = std::unique_ptr<RooFitResult>(static_cast<RooFitResult*>(r->Clone()));

    //
    // clean stuff

    nllval = nll->getVal();

    if(TRExFitter::DEBUGLEVEL >= 1) {
        std::streamsize ss = std::cout.precision();
        std::cout << std::fixed << std::setprecision(20);

        WriteInfoStatus("FittingTool::FitPDF", "***********************************************************");
        WriteInfoStatus("FittingTool::FitPDF", "  Final value of the NLL = " + std::to_string(nllval));
        WriteInfoStatus("FittingTool::FitPDF", "***********************************************************");

        std::cout << resetiosflags( ios::fixed | ios::showpoint );
        std::cout << std::setprecision(ss);
    }
    m_minNll = nllval;
    return nllval;
}

//____________________________________________________________________________________
//
void FittingTool::SaveFitResult( const std::string &fileName )
{
    std::unique_ptr<TFile> f(TFile::Open(fileName.c_str(),"RECREATE"));
    if (!f) {
        WriteWarningStatus("FittingTool::SaveFitResult", "Cannot open file: " + fileName + " cannot write the fit result");
        return;
    }
    m_fitResult->Write("",TObject::kOverwrite);
    f->Close();
}

//____________________________________________________________________________________
//
void FittingTool::ExportFitResultInTextFile( const std::string &fileName, const std::vector<std::string>& blinded )
{
    if(!m_fitResult){
        WriteErrorStatus("FittingTool::ExportFitResultInTextFile", "The FitResultObject seems not to be defined.");
        return;
    }
    
    //
    // Also saves fit result in root file with same name as txt file
    //
    TString fName = fileName;
    fName.ReplaceAll(".txt",".root");
    SaveFitResult(fName.Data());

    //
    // Printing the nuisance parameters post-fit values
    //
    ofstream nuisParAndCorr(fileName);
    nuisParAndCorr << "NUISANCE_PARAMETERS\n";

    for (auto var_tmp : m_fitResult->floatParsFinal()) {
        RooRealVar* var = static_cast<RooRealVar*>(var_tmp);

        // Not consider nuisance parameter being not associated to syst (yet)
        std::string varname = var->GetName();
        TString vname=var->GetName();
        vname.ReplaceAll("alpha_","");

        const double pull  = var->getVal(); // GetValue() return value in unit of sigma
        const double errorHi = var->getErrorHi();
        const double errorLo = var->getErrorLo();

        if (blinded.size() == 0){
            FittingTool::CheckUnderconstraint(var);
            nuisParAndCorr << vname << "  " << pull << " +" << fabs(errorHi) << " -" << fabs(errorLo)  << "\n";
        } else {
            std::string vname_s = vname.Data();
            if (std::find(blinded.begin(), blinded.end(), vname_s) == blinded.end()){
                FittingTool::CheckUnderconstraint(var);
                nuisParAndCorr << vname << "  " << pull << " +" << fabs(errorHi) << " -" << fabs(errorLo)  << "\n";
            } else {
                const std::string& hex = Common::DoubleToPseudoHex(pull);
                nuisParAndCorr << vname << "  " << hex << " +" << fabs(errorHi) << " -" << fabs(errorLo)  << "\n";
            }
        }
    }

    //
    // Correlation matrix
    //
    TH2* h2Dcorrelation = m_fitResult -> correlationHist();
    nuisParAndCorr << "\n\nCORRELATION_MATRIX\n";
    nuisParAndCorr << h2Dcorrelation->GetNbinsX() << "   " << h2Dcorrelation->GetNbinsY() << "\n";
    for(int kk=1; kk < h2Dcorrelation->GetNbinsX()+1; kk++) {
        for(int ll=1; ll < h2Dcorrelation->GetNbinsY()+1; ll++) {
            nuisParAndCorr << h2Dcorrelation->GetBinContent(kk,ll) << "   ";
        }
        nuisParAndCorr << "\n";
    }

    //
    // NLL value
    //

    nuisParAndCorr << std::fixed << std::setprecision(5);

    nuisParAndCorr << "\n\nNLL\n";
    nuisParAndCorr << m_minNll << "\n";

    //
    // Closing the output file
    //
    nuisParAndCorr << "\n";
    nuisParAndCorr.close();
}

//____________________________________________________________________________________
//
std::map < std::string, double > FittingTool::ExportFitResultInMap(){
    std::map < std::string, double > result;
    if(!m_fitResult){
        WriteErrorStatus("FittingTool::ExportFitResultInMap", "The FitResultObject seems not to be defined.");
        return result;
    }
    for (auto var_tmp : m_fitResult->floatParsFinal()) {
        RooRealVar* var = static_cast<RooRealVar*>(var_tmp);
        // Not consider nuisance parameter being not associated to syst
        std::string varname = var->GetName();
        const double pull  = var->getVal();
        result.insert( std::pair < std::string, double >(varname, pull) );
    }
    return result;
}

//____________________________________________________________________________________
//
void FittingTool::GetGroupedImpact(RooStats::ModelConfig* model,
                                   RooAbsPdf* fitpdf,
                                   RooAbsData* fitdata,
                                   RooWorkspace* ws,
                                   const std::string& categoryOfInterest,
                                   const std::string& outFileName,
                                   const std::string& fileName,
                                   const std::string& lumiLabel,
                                   const std::string& cmeLabel,
                                   const bool useHEPData) const {
    // obtain constrainedParams and poi like in FitPDF(), but make sure POI is not constant
    RooArgSet* constrainedParams = fitpdf->getParameters(*fitdata);
    RooStats::RemoveConstantParameters(constrainedParams);
    RooFit::Constrain(*constrainedParams);

    std::vector<RooRealVar*> pois = this->GetVectorPOI(model);
    std::vector<std::ofstream> outFiles(pois.size());

    if (pois.empty()) {
        WriteErrorStatus("FittingTool::GetGroupedImpact", "No POIs found!");
        return;
    }

    std::size_t iPOI(0);
    for (auto poi : pois) {
        poi -> setConstant(false);
        double value(0);
        const std::string name = poi->GetName();
        auto it = std::find_if(m_valPOIs.begin(), m_valPOIs.end(),
            [&name](const std::pair<std::string, double>& element){ return element.first == name;});

        if (it != m_valPOIs.end()) {
            value = it->second;
        }
        poi->setVal(value);
        if(!m_constPOI && m_randomize){
            poi->setVal(value + m_randomNP*(gRandom->Uniform(2)-1.) );
        }
        outFiles.at(iPOI).open((outFileName+"_"+name+".txt").c_str());
        ++iPOI;
    }

    // save snapshot of original workspace
    ws->saveSnapshot("snapshot_AfterFit_POI", *(model->GetParametersOfInterest()) );
    if (model->GetNuisanceParameters()) ws->saveSnapshot("snapshot_AfterFit_NP" , *(model->GetNuisanceParameters())   );
    ws->saveSnapshot("snapshot_AfterFit_GO" , *(model->GetGlobalObservables())    );

    std::vector<std::string> associatedParams; // parameters associated to a SubCategory

    // repeat the nominal fit - done so that the initial randomization is the exact same as for the following fit(s)
    // this should help avoid issues with fits ending up in different local minima for groups with very small impact on the POI
    FitExcludingGroup(false, false, fitdata, fitpdf, constrainedParams, model, ws, "Nominal", associatedParams);  // nothing held constant -> "snapshot_AfterFit_POI_Nominal"

    //
    // eventually do it once more
    if(TRExFitter::OPTION["GroupedImpactMoreFit"]>0){
        ws->saveSnapshot("snapshot_AfterFit_POI", *(model->GetParametersOfInterest()) );
        if (model->GetNuisanceParameters()) ws->saveSnapshot("snapshot_AfterFit_NP" , *(model->GetNuisanceParameters())   );
        ws->saveSnapshot("snapshot_AfterFit_GO" , *(model->GetGlobalObservables())    );
        FitExcludingGroup(false, false, fitdata, fitpdf, constrainedParams, model, ws, "Nominal", associatedParams);  // nothing held constant -> "snapshot_AfterFit_POI_Nominal"
    }

    // loop over unique SubCategories
    for (const auto& cat : m_subCategories) {
        if(categoryOfInterest!="all" && cat != categoryOfInterest) continue; // if a category was specified via command line, only process that one

        WriteInfoStatus("FittingTool::GetGroupedImpact","performing grouped systematics impact evaluation for: " + cat);

        // find all associated parameters per SubCategory
        associatedParams.clear();
        for(const auto& itSysts : m_subCategoryMap) {
            if (itSysts.second == cat) {
                associatedParams.push_back(itSysts.first);
            }
        }

        // special case for gammas
        if(cat == "Gammas") {
            FitExcludingGroup(true,  false, fitdata, fitpdf, constrainedParams, model, ws, cat, associatedParams);
        }

        // special case for stat-only fit
        else if(cat == "FullSyst") {
            FitExcludingGroup(true,  true,  fitdata, fitpdf, constrainedParams, model, ws, cat, associatedParams);
        }

        // default: perform a fit where parameters in SubCategory are held constant
        else {
            FitExcludingGroup(false, false, fitdata, fitpdf, constrainedParams, model, ws, cat, associatedParams);
        }
    }

    // load original workspace again
    ws->loadSnapshot("snapshot_AfterFit_GO");
    ws->loadSnapshot("snapshot_AfterFit_POI");
    if (model->GetNuisanceParameters()) {
        ws->loadSnapshot("snapshot_AfterFit_NP");
    }

    WriteInfoStatus("FittingTool::GetGroupedImpact","-----------------------------------------------------");

    for (std::size_t ipoi(0); ipoi < pois.size(); ++ipoi) {
        const std::string& name = pois.at(ipoi)->GetName();
        // report replication of nominal fit
        ws->loadSnapshot("snapshot_AfterFit_POI_Nominal");
        WriteInfoStatus("FittingTool::GetGroupedImpact", "replicated nominal fit");
        WriteInfoStatus("FittingTool::GetGroupedImpact", "POI: " + name);
        WriteInfoStatus("FittingTool::GetGroupedImpact", "     " + std::to_string(pois.at(ipoi)->getVal()) + " +/- " + std::to_string(pois.at(ipoi)->getError()) +
                                                         "    ( +" + std::to_string(pois.at(ipoi)->getErrorHi()) + ", " + std::to_string(pois.at(ipoi)->getErrorLo()) + " )");

        const double NomUp2=(pois.at(ipoi)->getErrorHi()*pois.at(ipoi)->getErrorHi());
        const double NomLo2=(pois.at(ipoi)->getErrorLo()*pois.at(ipoi)->getErrorLo());
        const double Nom2  =(pois.at(ipoi)->getError()*pois.at(ipoi)->getError());

        std::vector<YamlConverter::ImpactContainer> container;

        // report impact calculations, impact is obtained by quadrature subtraction from replicated nominal fit
        for (const auto& cat : m_subCategories) {
            if(categoryOfInterest!="all" && cat != categoryOfInterest) continue; // if a category was specified via command line, only process that one

            ws->loadSnapshot(("snapshot_AfterFit_POI_" + cat).c_str());
            WriteInfoStatus("FittingTool::GetGroupedImpact","-----------------------------------------------------");
            WriteInfoStatus("FittingTool::GetGroupedImpact", "category: " + cat + " (fixed to best-fit values for fit)");
            WriteInfoStatus("FittingTool::GetGroupedImpact", "POI is:   " + std::to_string(pois.at(ipoi)->getVal()) + " +/- " + std::to_string(pois.at(ipoi)->getError()) +
                                                             "    ( +" + std::to_string(pois.at(ipoi)->getErrorHi()) + ", " + std::to_string(pois.at(ipoi)->getErrorLo()) + " )");
            if (cat == "FullSyst") WriteDebugStatus("FittingTool::GetGroupedImpact", "  (corresponds to a stat-only fit)");
            const double impact = -(pois.at(ipoi)->getError()*pois.at(ipoi)->getError()) + Nom2;
            if (impact < 0) {
                WriteWarningStatus("FittingTool::GetGroupedImpact", "NaN impact encountered");
                WriteWarningStatus("FittingTool::GetGroupedImpact", "  - Nominal error squared: " + std::to_string(Nom2));
                WriteWarningStatus("FittingTool::GetGroupedImpact", "  - Category fixed error squared: " + std::to_string(pois.at(ipoi)->getError()*pois.at(ipoi)->getError()));
                WriteWarningStatus("FittingTool::GetGroupedImpact", "  - Nominal error squared minus category error squared: " + std::to_string(impact));
            }
            // cppcheck-suppress invalidFunctionArg
            WriteInfoStatus("FittingTool::GetGroupedImpact", "           --> impact: " + std::to_string(std::sqrt(impact)) +
                                                             "    ( +" + std::to_string(std::sqrt(- (pois.at(ipoi)->getErrorHi()*pois.at(ipoi)->getErrorHi()) + NomUp2)) + ", -" + std::to_string(std::sqrt(-(pois.at(ipoi)->getErrorLo()*pois.at(ipoi)->getErrorLo()) + NomLo2)) + " )" );

            // write results to file
            // cppcheck-suppress invalidFunctionArg
            outFiles.at(ipoi) << cat << "    " << std::sqrt(impact) << "  ( +" << std::sqrt( - (pois.at(ipoi)->getErrorHi()*pois.at(ipoi)->getErrorHi()) + NomUp2 ) << ", -" << std::sqrt( - (pois.at(ipoi)->getErrorLo()*pois.at(ipoi)->getErrorLo()) + NomLo2 ) << " )\n";

            YamlConverter::ImpactContainer imp;
            imp.name    = cat;
            // cppcheck-suppress invalidFunctionArg
            imp.error   = std::sqrt(impact);
            imp.errorHi = std::sqrt(-(pois.at(ipoi)->getErrorLo()*pois.at(ipoi)->getErrorLo()) + NomLo2);
            imp.errorLo = std::sqrt(-(pois.at(ipoi)->getErrorHi()*pois.at(ipoi)->getErrorHi()) + NomUp2);
            container.emplace_back(std::move(imp));
        }

        WriteInfoStatus("FittingTool::GetGroupedImpact", "-----------------------------------------------------");

        YamlConverter converter{};
        converter.WriteImpact(container, outFileName+"_"+name+".yaml");
        if (useHEPData) {
            converter.SetLumi(Common::ReplaceString(lumiLabel, " fb^{-1}", ""));
            converter.SetCME(Common::ReplaceString(cmeLabel, " TeV", "000"));
            converter.WriteImpactHEPData(container, fileName, "_"+name);
        }

        outFiles.at(ipoi).close();
    }

    // load original workspace again
    ws->loadSnapshot("snapshot_AfterFit_GO");
    ws->loadSnapshot("snapshot_AfterFit_POI");
    if (model->GetNuisanceParameters()) {
        ws->loadSnapshot("snapshot_AfterFit_NP");
    }
}

//____________________________________________________________________________________
//
// perform a fit where all parameters in "affectedParams" (usually coming from SubGroup "category") are set to constant, optionally also gammas or all parameters
void FittingTool::FitExcludingGroup(bool excludeGammas, bool statOnly, RooAbsData*& fitdata, RooAbsPdf*& fitpdf, RooArgSet*& constrainedParams,
                                    RooStats::ModelConfig* mc, RooWorkspace* ws, const std::string& category, const std::vector<std::string>& affectedParams) const {

    if (!mc->GetNuisanceParameters()) return;
    // (VD): use this to fix nuisance parameter before the fit
    const RooArgSet* glbObs = mc->GetGlobalObservables();
    ws->loadSnapshot("snapshot_AfterFit_GO");
    ws->loadSnapshot("snapshot_AfterFit_POI");
    ws->loadSnapshot("snapshot_AfterFit_NP");

    WriteInfoStatus("FittingTool::FitExcludingGroup", "-----------------------------------------------------");
    WriteInfoStatus("FittingTool::FitExcludingGroup", "           breakdown for " + category);
    WriteInfoStatus("FittingTool::FitExcludingGroup", "-----------------------------------------------------");

    for (auto var_tmp : *mc->GetNuisanceParameters()) {
        RooRealVar* var = static_cast<RooRealVar*>(var_tmp);
        std::string varname = var->GetName();

        // default: set everything non-constant (but the saturated-model norm-factors!)
        if (varname.find("saturated_model_sf_")==std::string::npos) var->setConstant(0);
        else var->setConstant(1);

        // if excludeGammas==true, set gammas to constant
        if (excludeGammas) {
            if (varname.find("gamma_stat")!=string::npos) var->setConstant(1);
        }

        // set all affectedParams constant
        if (std::find(affectedParams.begin(), affectedParams.end(), varname) != affectedParams.end()) {
            var->setConstant(1);
        }

        // for stat-only fits, set everything constant
        if (statOnly) {
            var->setConstant(1);
        }
    }

    // repeat the fit here ....
    std::unique_ptr<RooAbsReal> nll(fitpdf->createNLL(*fitdata,
                                    RooFit::Constrain(*constrainedParams),
                                    RooFit::GlobalObservables(*glbObs),
                                    RooFit::Offset(1),
                                    NumCPU(m_CPU,RooFit::Hybrid),
                                    RooFit::Optimize(kTRUE),
                                    RooFit::ExternalConstraints(*m_externalConstraints)
                                   ));
    
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(1);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(-1);
    const TString minimType = ::ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str();
    const TString algorithm = ::ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str();
    const double tol =        ::ROOT::Math::MinimizerOptions::DefaultTolerance(); //AsymptoticCalculator enforces not less than 1 on this
    
    RooMinimizer minim2(*nll);
    minim2.optimizeConst(2);
    minim2.setStrategy(1);
    minim2.setMinimizerType(minimType);
    minim2.setPrintLevel(1);
    minim2.setEps(tol);
    const int status = minim2.minimize(minimType.Data(),algorithm.Data());
    const bool HessStatus = minim2.hesse();
    
    std::vector<RooRealVar*> pois = GetVectorPOI(mc);

    bool first(true);
    for (auto poi : pois) {
        RooArgSet minosSet(*poi);
        if(m_useMinos && (m_varMinos.at(0)=="all" || std::find(m_varMinos.begin(), m_varMinos.end(), poi->GetName()) != m_varMinos.end())){
            minim2.minos(minosSet);
        }

        if (status!=0) WriteErrorStatus("FittingTool::FitExcludingGroup", "unable to perform fit correctly! HessStatus: " + std::to_string(HessStatus));

        const double newPOIerr =poi->getError();
        const double newPOIerrU=poi->getErrorHi();
        const double newPOIerrD=poi->getErrorLo();

        const std::string snapshotName = "snapshot_AfterFit_POI_" + category;
        if (first) ws->saveSnapshot(snapshotName.c_str(), *mc->GetParametersOfInterest() );

        ws->loadSnapshot("snapshot_AfterFit_POI_Nominal");
        const double oldPOIerr =poi->getError();
        const double oldPOIerrU=poi->getErrorHi();
        const double oldPOIerrD=poi->getErrorLo();

        // check if uncertainties have increased compared to nominal fit, with 0.5% tolerance
        if ( (std::fabs(newPOIerrU)>std::fabs(oldPOIerrU)*1.005) || (std::fabs(newPOIerrD)>std::fabs(oldPOIerrD)*1.005) ) {
            const std::string& name = poi->GetName();
            WriteWarningStatus("FittingTool::FitExcludingGroup", "POI: " + name);
            WriteWarningStatus("FittingTool::FitExcludingGroup", "uncertainty has increased for " + category + "! please check the fit");
            WriteWarningStatus("FittingTool::FitExcludingGroup", "old: " + std::to_string(oldPOIerr) + " (+" + std::to_string(oldPOIerrU) + ", " + std::to_string(oldPOIerrD) + ")");
            WriteWarningStatus("FittingTool::FitExcludingGroup", "new: " + std::to_string(newPOIerr) + " (+" + std::to_string(newPOIerrU) + ", " + std::to_string(newPOIerrD) + ")");
        }
        first = false;
    }
}

//____________________________________________________________________________________
//
// Check for underconstraints
void FittingTool::CheckUnderconstraint(const RooRealVar* const var) const {
    const std::string name = var->GetName();
    const double errorHi = var->getErrorHi();
    const double errorLo = var->getErrorLo();

    // dont check gamma parameters
    if (name.find("alpha_") == std::string::npos) return;

    if (errorHi > 1.001 || errorLo < -1.001){
        WriteWarningStatus("FittingTool::CheckUnderconstraint","NuisanceParameter: " + name + " is underconstrained! This may indicate fit convergence problems!");
    }
}

//____________________________________________________________________________________
//
std::vector<RooRealVar*> FittingTool::GetVectorPOI(const RooStats::ModelConfig* model) const {
    std::vector<RooRealVar*> result;
    for (auto var_tmp : *model->GetParametersOfInterest()) {
        auto tmp = dynamic_cast<RooRealVar*>(var_tmp);
        if (!tmp) {
            WriteErrorStatus("FittingTool::GetVectorPOI", "Cannot find the parameter of interest !");
            result.clear();
            return result;
        }

        const std::string name = tmp->GetName();
        auto it = std::find_if(result.begin(), result.end(),
            [&name](const RooRealVar* var){return name == var->GetName();});

        if (it == result.end()) {
            result.emplace_back(tmp);
        }
    }
    return result;
}

//__________________________________________________________________________________
//
void FittingTool::PrintMinuitHelp() const {
    if (TRExFitter::DEBUGLEVEL < 2) return;
    WriteInfoStatus("FittingTool::PrintMinuitHelp", "");
    WriteInfoStatus("FittingTool::PrintMinuitHelp", "------------------------------------------");
    WriteInfoStatus("FittingTool::PrintMinuitHelp", "         Status code cheat-sheet          ");
    WriteInfoStatus("FittingTool::PrintMinuitHelp", "------------------------------------------");
    WriteInfoStatus("FittingTool::PrintMinuitHelp", " status = MIGRAD status + 100*Hesse status");
    WriteInfoStatus("FittingTool::PrintMinuitHelp", " MIGRAD status");
    WriteInfoStatus("FittingTool::PrintMinuitHelp", "    - 0 = No problems");
    WriteInfoStatus("FittingTool::PrintMinuitHelp", "    - 1 = Covariance was made pos defined");
    WriteInfoStatus("FittingTool::PrintMinuitHelp", "    - 2 = Hesse is invalid");
    WriteInfoStatus("FittingTool::PrintMinuitHelp", "    - 3 = Edm is above max");
    WriteInfoStatus("FittingTool::PrintMinuitHelp", "    - 4 = Reached call limit");
    WriteInfoStatus("FittingTool::PrintMinuitHelp", "    - 5 = Any other failure");
    WriteInfoStatus("FittingTool::PrintMinuitHelp", " Hesse status");
    WriteInfoStatus("FittingTool::PrintMinuitHelp", "    - 0 = No problems");
    WriteInfoStatus("FittingTool::PrintMinuitHelp", "    - 1 = Hesse failed");
    WriteInfoStatus("FittingTool::PrintMinuitHelp", "    - 2 = Matrix inversion failed");
    WriteInfoStatus("FittingTool::PrintMinuitHelp", "    - 3 = Matrix is not pos defined");
    WriteInfoStatus("FittingTool::PrintMinuitHelp", "");
}
