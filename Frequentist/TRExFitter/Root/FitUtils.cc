#include "TRExFitter/FitUtils.h"

#include "TRExFitter/FittingTool.h"
#include "TRExFitter/NormFactor.h"
#include "TRExFitter/StatusLogbook.h"

#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooMultiVarGaussian.h"
#include "RooRealSumPdf.h"
#include "RooRealVar.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"
#include "RooTFnBinding.h"

#include "TDirectory.h"
#include "TFile.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include "TF3.h"

#include <sstream>

//__________________________________________________________________________________
//
void FitUtils::ApplyExternalConstraints(RooWorkspace* ws,
                                        FittingTool* fitTool,
                                        RooSimultaneous* simPdf,
                                        const std::vector<std::shared_ptr<NormFactor> >& normFactors,
                                        int type
                                       ) {

    // Tikhonov regularization (for unfolding)
    // only add the unique NFs
    std::vector<std::string> names;
    std::vector<double> tauVec;
    std::vector<std::shared_ptr<NormFactor> > uniqueNF;
    for(const auto& nf : normFactors) {
        if (std::find(names.begin(), names.end(), nf->fName) != names.end()) continue;
        names.emplace_back(nf->fName);
        uniqueNF.emplace_back(nf);
    }
    
    // Simple form (individual constraint terms on signal strenghts)
    if(type==0){
        RooArgList l;
        std::vector<double> nomVec;
        for(const auto& nf : uniqueNF) {
            if(nf->fTau == 0) continue;
            
            if(ws->var(nf->fName.c_str())) {
                l.add(*ws->var(nf->fName.c_str()));
            } else if (ws->function(nf->fName.c_str())) {
                l.add(*ws->function(nf->fName.c_str()));
            } else {
                WriteWarningStatus("FitUtils::ApplyExternalConstraints","Cannot apply tau to norm factor " + nf->fName);
                continue;
            }
            
            nomVec.push_back( nf->fNominal );
            tauVec.push_back( nf->fTau );
        }

        if(tauVec.empty()) return;

        TVectorD nominal(nomVec.size());
        TMatrixDSym cov(tauVec.size());
        for(unsigned int i_tau=0;i_tau<tauVec.size();i_tau++){
            nominal(i_tau) = nomVec[i_tau];
            cov(i_tau,i_tau) = (1./tauVec[i_tau]) * (1./tauVec[i_tau]);
        }
        RooMultiVarGaussian r("regularization","regularization",l,nominal,cov);
        ws->import(r);
    }
    // Discretized second derivative
    else if(type==1){
        RooArgList l_prime;
        for(unsigned int i_nf=0;i_nf<uniqueNF.size();i_nf++){
            if(i_nf==0 || i_nf==uniqueNF.size()-1 ) continue;
            if(uniqueNF[i_nf]->fTau == 0) continue;
            tauVec.push_back(uniqueNF[i_nf]->fTau);

            TF3 *f3 = new TF3(Form("f3_bin%d",i_nf+1),"x - 2*y + z", 
                    uniqueNF[i_nf-1]->fMin,uniqueNF[i_nf-1]->fMax, 
                    uniqueNF[i_nf]  ->fMin,uniqueNF[i_nf]  ->fMax, 
                    uniqueNF[i_nf+1]->fMin,uniqueNF[i_nf+1]->fMax
                    );
            RooAbsReal *mu_prime = RooFit::bindFunction( f3, 
                                                         *ws->function(uniqueNF[i_nf-1]->fName.c_str()), 
                                                         *ws->function(uniqueNF[i_nf]  ->fName.c_str()), 
                                                         *ws->function(uniqueNF[i_nf+1]->fName.c_str())
                                                        );
            
            mu_prime->SetName(Form("mu_prime_%d",i_nf+1));
            ws->import(*mu_prime);
            l_prime.add(*ws->function(Form("mu_prime_%d",i_nf+1)));
        }
        
        TVectorD nominal(tauVec.size());
        TMatrixDSym cov(tauVec.size());
        for(unsigned int i_tau=0;i_tau<tauVec.size();i_tau++){
            nominal(i_tau) = 0;
            cov(i_tau,i_tau) = (1./tauVec[i_tau]) * (1./tauVec[i_tau]);
        }
        RooMultiVarGaussian r("regularization","regularization",l_prime,nominal,cov);
        ws->import(r);
    }
    
    ws->defineSet("myConstraints","regularization");
    simPdf->setStringAttribute("externalConstraints","myConstraints");

    if(simPdf->getStringAttribute("externalConstraints")){
        WriteInfoStatus("FitUtils::ApplyExternalConstraints",Form("Building NLL with external constraints %s",simPdf->getStringAttribute("externalConstraints")));
        const RooArgSet* externalConstraints = ws->set(simPdf->getStringAttribute("externalConstraints"));
        fitTool->SetExternalConstraints( externalConstraints );
    }
}

//__________________________________________________________________________________
//
void FitUtils::SetBinnedLikelihoodOptimisation(RooWorkspace* ws) {
    for (auto arg : ws->components()) {
        if (arg->IsA() == RooRealSumPdf::Class()) {
            arg->setAttribute("BinnedLikelihood");
            const std::string temp_string = arg->GetName();
            WriteDebugStatus("FitUtils::SetBinnedLikelihoodOptimisation", "Activating binned likelihood attribute for " + temp_string);
        }
    }
}

//__________________________________________________________________________________
//
void FitUtils::InjectGlobalObservables(RooWorkspace* ws,
                                       const std::map< std::string, double >& npValues) {
    RooStats::ModelConfig* mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    if (!mc) {
        WriteErrorStatus("FitUtils::InjectGlobalObservables", "Cannot read ModelCOnfig");
        exit(EXIT_FAILURE);
    }
    RooArgSet mc_globs = *mc->GetGlobalObservables();

    WriteInfoStatus("FitUtils::InjectGlobalObservables", "Injecting the following NP values to global observables");
    for(const auto& np_value : npValues) {
        const std::string this_name = np_value.first;
        double this_value = np_value.second;
    
        std::ostringstream tmp;
    
        // find the corresponding glob
        const std::string glob_name = "nom_" + this_name;
        const std::string glob_name_alpha = "nom_alpha_" + this_name;
        const std::string glob_name_gamma = "nom_gamma_" + this_name;
        RooRealVar* this_glob = nullptr;
        for ( auto glob_tmp : mc_globs) {
            RooRealVar* glob = static_cast<RooRealVar*>(glob_tmp);
            if(glob->GetName() == glob_name || glob->GetName() == glob_name_alpha || glob->GetName() == glob_name_gamma) {
                this_glob = glob;
                break;
            }
        }

        if(!this_glob) {
            WriteWarningStatus("FitUtils::InjectGlobalObservables", "Could not find global observable "+glob_name);
            continue;
        }

        // set gamma values to gamma*nom_gamma
        // cppcheck-suppress stlIfStrFind
        if(glob_name.find("nom_gamma_") == 0) {
            this_value = this_value * this_glob->getVal();
        }

        tmp << this_glob->GetName() << ": " << this_value;
        WriteInfoStatus("FitUtils::InjectGlobalObservables", tmp.str());
        this_glob->setVal(this_value);
    }
}

//__________________________________________________________________________________
//
void FitUtils::DisableSaturatedModel(RooWorkspace* ws) {
    RooArgSet vars = ws->allVars();
    for(auto var_tmp : vars) {
        RooRealVar* var = static_cast<RooRealVar*>(var_tmp);
        const std::string& name = var->GetName();
        if(name.find("saturated_model_sf_")!=std::string::npos){
            WriteInfoStatus("FitUtils::DisableSaturatedModel","Fixing parameter " + name );
            var->setConstant( 1 );
        }
    }
}

//__________________________________________________________________________________
//
void FitUtils::SetPOIinFile(const std::string& path, const std::string& wsName, const std::string& poi) {
    std::unique_ptr<TFile> f(TFile::Open(path.c_str()));
    if (!f) {
        WriteErrorStatus("FitUtils::SetPOIinFile", "Cannot open input file at: " + path);
        return;
    }

   RooWorkspace* ws = dynamic_cast<RooWorkspace*>(f->Get(wsName.c_str()));
    if (!ws) {
        WriteErrorStatus("FitUtils::SetPOIinFile", "Cannot read workspace");
        return;
    }

    RooStats::ModelConfig* mc = dynamic_cast<RooStats::ModelConfig*>(ws->obj("ModelConfig"));
    if (!mc) {
        WriteErrorStatus("FitUtils::SetPOIinFile", "Cannot read ModelConfig");
        return;
    }

    mc->SetParametersOfInterest(poi.c_str());

    f->Close();
}

//__________________________________________________________________________________
//
std::vector<std::string> FitUtils::GetAllParameters(const RooStats::ModelConfig* mc) {

    if (!mc) {
        WriteErrorStatus("FitUtils::GetAllParameters", "ModelCOnfig is nullptr");
        exit(EXIT_FAILURE);
    }

    std::vector<std::string> params;
    for (const auto& i : *mc->GetParametersOfInterest()) {
        params.emplace_back(i->GetName());
    }
    
    if (mc->GetNuisanceParameters()) {
        for (const auto& i : *mc->GetNuisanceParameters()) {
            params.emplace_back(i->GetName());
        }
    }

    return params;
}


//__________________________________________________________________________________
//
std::size_t FitUtils::NumberOfFreeParameters(const RooStats::ModelConfig* mc) {
    if (!mc) {
        WriteErrorStatus("FitUtils::NumberOfFreeParameters", "ModelCOnfig is nullptr");
        exit(EXIT_FAILURE);
    }

    std::size_t result(0);

    for (const auto& i : *mc->GetParametersOfInterest()) {
        if (!i->isConstant()) ++result;
    }
    
    if (mc->GetNuisanceParameters()) {
        for (const auto& i : *mc->GetNuisanceParameters()) {
            if (!i->isConstant()) ++result;
        }
    }

    return result;
}

//__________________________________________________________________________________
//
void FitUtils::FixAllParameters(RooStats::ModelConfig* mc, const std::vector<std::string>& exceptions) {
    if (!mc) {
        WriteErrorStatus("FitUtils::FixAllParameters", "ModelCOnfig is nullptr");
        exit(EXIT_FAILURE);
    }

    for (const auto& i : *mc->GetParametersOfInterest()) {
        auto* tmp = static_cast<RooRealVar*>(i);
        tmp->setConstant(1);
    }
    
    if (mc->GetNuisanceParameters()) {
        for (const auto& i : *mc->GetNuisanceParameters()) {
            auto* tmp = static_cast<RooRealVar*>(i);
            const std::string name = tmp->GetName();

            auto it = std::find_if(exceptions.begin(), exceptions.end(), 
                [&name](const std::string& element){return name.find(element) != std::string::npos;});

            if (it != exceptions.end()) continue;
            tmp->setConstant(1);
        }
    }
}

//__________________________________________________________________________________
//
void FitUtils::FloatAllParameters(RooStats::ModelConfig* mc, const bool saturatedModel) {
    if (!mc) {
        WriteErrorStatus("FitUtils::FloatAllParameters", "ModelCOnfig is nullptr");
        exit(EXIT_FAILURE);
    }

    for (const auto& i : *mc->GetParametersOfInterest()) {
        auto* tmp = static_cast<RooRealVar*>(i);
        tmp->setConstant(0);
    }
    
    if (mc->GetNuisanceParameters()) {
        for (const auto& i : *mc->GetNuisanceParameters()) {
            auto* tmp = static_cast<RooRealVar*>(i);
            const std::string name = tmp->GetName();
            if (saturatedModel) {
                if (name.find("saturated_model_sf") != std::string::npos) {
                    tmp->setVal(1);
                    tmp->setConstant(1);
                    continue;
                }
            }
            tmp->setConstant(0);
        }
    }
}

//__________________________________________________________________________________
//
std::vector<double> FitUtils::CalculateExpressionRoofit(RooWorkspace* ws, TFile* fitResultFile, const std::string& name) {

    if (!ws || !fitResultFile) {
        WriteErrorStatus("FitUtils::CalculateExpressionRoofit", "Passed nullptr");
        exit(EXIT_FAILURE);
    }

    std::vector<double> result;

    std::unique_ptr<RooFitResult> fr(nullptr);
    for(auto *key : *fitResultFile->GetListOfKeys()){
        const std::string& keyName = key->GetName();
        if(keyName.find("nll_simPdf")!=std::string::npos){
            fr.reset(static_cast<RooFitResult*>(fitResultFile->Get(keyName.c_str())));
            break;
        }
    }

    if (!fr) {
        WriteErrorStatus("FitUtils::CalculateExpressionRoofit", "FitResult is nullptr");
        exit(EXIT_FAILURE);
    }

    std::unique_ptr<RooFormulaVar> muN(nullptr);
    for(auto *arg : ws->allFunctions()) {
        const std::string argName = arg->GetName();
        if(argName == name) {
            muN.reset(static_cast<RooFormulaVar*>(arg));
        }
    }

    if (!muN) {
        WriteWarningStatus("FitUtils::CalculateExpressionRoofit", "Did not find the correct name");
        return result;
    }
    
    for(auto *par : fr->floatParsFinal()){
        ws->var( par->GetName() )->setVal( (static_cast<RooRealVar*>(par))->getVal() );
        ws->var( par->GetName() )->setError( (static_cast<RooRealVar*>(par))->getError() );
    }

    const double mean  = muN->getVal();
    const double error = muN->getPropagatedError(*fr);
    const double up    = mean + error;
    const double down  = mean - error;

    result.emplace_back(mean);
    result.emplace_back(up);
    result.emplace_back(down);
    result.emplace_back(error);

    return result;
}

//__________________________________________________________________________________
//
void FitUtils::DisableSaturatedModelFileAndSetPOI(const std::string& path,
                                                  const std::string& wsName,
                                                  const std::string& wsOutName,
                                                  const std::string& poi) {
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "UPDATE"));
    if (!f) {
        WriteWarningStatus("FitUtils::DisableSaturatedModelFile", "Cammpt open file " + path);
        return;
    }
    RooWorkspace* workspace = dynamic_cast<RooWorkspace*>(f->Get(wsName.c_str()));
    if (!workspace) {
        WriteWarningStatus("FitUtils::DisableSaturatedModelFile", "Cannot read WS " + wsName);
        return;
    }
    FitUtils::DisableSaturatedModel(workspace);
    RooStats::ModelConfig* mc = dynamic_cast<RooStats::ModelConfig*>(workspace->obj("ModelConfig"));
    if (!mc) {
        WriteWarningStatus("FitUtils::DisableSaturatedModelFileAndSetPOI", "Cannot read ModelConfig");
        return;
    }

    mc->SetParametersOfInterest(poi.c_str());
    workspace->Write(wsOutName.c_str());
    f->Close();
}

//__________________________________________________________________________________
//
RooArgSet FitUtils::GetPOIsWithout(const RooStats::ModelConfig* mc, const std::vector<std::string>& names) {
    RooArgSet result;
    if (mc->GetNuisanceParameters()) {
        for (const auto i : *mc->GetNuisanceParameters()) {
            auto* tmp = static_cast<RooRealVar*>(i);
            const std::string name = tmp->GetName();
            auto it = std::find(names.begin(), names.end(), name);
            if (it != names.end()) continue;

            result.add(*i);
        }
     }

    return result;
}

