// Class include
#include "TRExFitter/EFTProcessor.h"

// Framework includes
#include "TRExFitter/Common.h"
#include "TRExFitter/HistoTools.h"
#include "TRExFitter/NormFactor.h"
#include "TRExFitter/Region.h"
#include "TRExFitter/Sample.h"
#include "TRExFitter/StatusLogbook.h"

// ROOT includes
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TStyle.h"
#include "TSystem.h"

// c++ includes
#include <memory>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>


using namespace std;


// -------------------------------------------------------------------------------------------------
// SampleHist

//_____________________________________________________________________________
//
EFTProcessor::EFTProcessor() :
    fName(""),
    fTitle(""),
    fRefMap{},
    fExtrapMap{} {
}

EFTProcessor::EFTProcessor(std::string param, std::string title) :
    fName(param),
    fTitle(title),
    fRefMap{},
    fExtrapMap{} {
}

//_____________________________________________________________________________
//
EFTProcessor::~EFTProcessor(){
}

//_____________________________________________________________________________
//
EFTProcessor::EFTExtrap::EFTExtrap() :
    nBins(0),
    p0s{},
    p1s{},
    p2s{} {
}


//_____________________________________________________________________________
//
void EFTProcessor::Print() const {
    for (const auto& imap : fRefMap) {
        WriteDebugStatus("EFTProcessor::Print", "  - "+ imap.first);
        for(const auto& eft_s : imap.second){
	        WriteDebugStatus("EFTProcessor::Print", "    - "+eft_s);
        }
    }
}

//_____________________________________________________________________________
// this draws the control plots for each EFT param-sample combination
void EFTProcessor::DrawEFTInputs(std::vector < Region* > Regions) const {

    for(const auto& ireg : Regions) {
        const std::string rname = ireg->fName;
        WriteDebugStatus("EFTProcessor::DrawEFTInputs", " Region: " + rname);

        for (const auto& imap : fRefMap) {
            WriteDebugStatus("EFTProcessor::DrawEFTInputs", "    SM Ref: " + imap.first);
            ////////////////////////
            // Get and draw SMRef and EFT histograms

            // Get SM Reference sample
            std::shared_ptr<SampleHist> SM_sh = ireg->GetSampleHist(imap.first);
            if( !SM_sh ){
                WriteDebugStatus("EFTProcessor::DrawEFTInputs", "     - Skipping region");
                continue;
            }

            // Make output dir
            gSystem->mkdir((SM_sh->fFitName+"/EFT").c_str());

            std::unique_ptr<TH1> SM_hist(static_cast<TH1*>(SM_sh->fHist->Clone(TString::Format("tmp_%s",imap.first.c_str()))));
            SM_hist->SetDirectory(nullptr);
            const std::string SMRef = imap.first;

            std::vector<Int_t> cols;
            cols.emplace_back(kRed+1);
            cols.emplace_back(kRed-4);
            cols.emplace_back(kGreen+1);
            cols.emplace_back(kGreen-4);
            cols.emplace_back(kBlue+1);
            cols.emplace_back(kBlue-4);
            cols.emplace_back(kMagenta+1);
            cols.emplace_back(kMagenta-4);
            //@TODO be more clever about this...


            // Plotting
            TCanvas c1{};

            // Plot SM hist
            SM_hist->Draw("HIST");
            SM_hist->SetMaximum(SM_hist->GetMaximum()*2);
            SM_hist->SetLineWidth(3);
            SM_hist->SetFillColor(0);
            SM_hist->SetLineColor(1);
            SM_hist->GetXaxis()->SetTitle(SM_sh->fVariableTitle.c_str());
            SM_hist->GetYaxis()->SetTitle("Number of events");

            TLatex tex{};
            tex.SetNDC();
            tex.DrawLatex(0.2,0.79,TString::Format("%s, %s",fTitle.c_str(),SMRef.c_str()));
            tex.DrawLatex(0.2,0.72,rname.c_str());


            std::unique_ptr<TLegend> leg = std::make_unique<TLegend>(0.6,0.7,0.9,0.87);
            leg->SetBorderSize(0);
            leg->SetNColumns(2);

            Int_t n=0;
            double ymax=-999;

            // Cycle through all EFT variations and plot - each corresponds to a differrent coefficient value
            for(const auto& EFT_name : imap.second){     //@TODO: make order coherent for legend
                WriteDebugStatus("EFTProcessor::DrawEFTInputs", "    Working on: " + EFT_name);

                std::shared_ptr<SampleHist> EFT_sh = ireg->GetSampleHist(EFT_name);

                std::unique_ptr<TH1> EFT_hist(static_cast<TH1*>(EFT_sh->fHist->Clone()));
                EFT_hist->SetDirectory(nullptr);

                EFT_hist->SetLineColor(cols[n]);
                EFT_hist->SetLineWidth(3);

                if(EFT_sh->GetSample()->fEFTValue>0)EFT_hist->SetLineStyle(kDashed);
                else EFT_hist->SetLineStyle(kDotted);

                if(EFT_hist->GetMaximum() > ymax)ymax=EFT_hist->GetMaximum();

                EFT_hist->DrawCopy("HISTSAME");
                leg->AddEntry(EFT_hist->Clone(TString::Format("tmp_%s",EFT_name.c_str())),TString::Format("%s = %.1f",fTitle.c_str(),EFT_sh->GetSample()->fEFTValue),"l");
                n++;
            }

            SM_hist->SetMaximum(ymax*1.5);

            leg->Draw();

            for(const auto& format : TRExFitter::IMAGEFORMAT) {
            	c1.SaveAs(TString::Format("%s/EFT/%s_%s_%s.%s",SM_sh->fFitName.c_str(),rname.c_str(),fName.c_str(),SMRef.c_str(),format.c_str()));
            }
        }
    }
}

//_____________________________________________________________________________
// this fits the EFT quadratic bin-by-bin
void EFTProcessor::FitEFTInputs(std::vector < Region* > Regions, const std::string& fileName) {

    // Open file for quadratic fit result output
    std::ofstream NFTXTFile(fileName.c_str());
    NFTXTFile << "NORMFACTORS\n";

    for(const auto& ireg : Regions) {
        const std::string rname = ireg->fName;
        WriteDebugStatus("EFTProcessor::FitEFTInputs", " Region: " + rname);

        for (const auto& imap : fRefMap) {
            WriteDebugStatus("EFTProcessor::FitEFTInputs", "    SM Ref: " + imap.first);
            ////////////////////////
            // Get and draw SMRef and EFT histograms

            // Get SM Reference sample
            std::shared_ptr<SampleHist> SM_sh = ireg->GetSampleHist(imap.first);
            if( !SM_sh ){
                WriteDebugStatus("EFTProcessor::FitEFTInputs", "     - Skipping region");
                continue;
            }

            // Make output dir
            gSystem->mkdir((SM_sh->fFitName+"/EFT").c_str());

            std::unique_ptr<TH1> SM_hist(static_cast<TH1*>(SM_sh->fHist->Clone(TString::Format("tmp_%s",imap.first.c_str()))));
            const std::string SMRef = imap.first;

            double max_rel_eft=-999;
            double min_c=999;
            double max_c=-999;

            ////////////////////////
            // Fitting quadratic for each bin

            double canv_x = std::min((SM_hist->GetNbinsX()+1)*200,2000);
            TCanvas b1("b1","b1",canv_x,400);
            b1.Divide(SM_hist->GetNbinsX()+1,1,0,0);

            b1.cd(1);
            // Global EFT param, Region and Sample labels
            TLatex btex{};
            btex.SetNDC();
            btex.SetTextSize(20);
            btex.DrawLatex(0.17,0.9,fTitle.c_str());
            btex.DrawLatex(0.17,0.79,SMRef.c_str());
            btex.DrawLatex(0.17,0.72,rname.c_str());

            // Draw legend (need some dummy objects)
            TH1D hleg("", "", 1, 0, 1);
            hleg.SetLineColor(kRed);
            hleg.SetLineWidth(2);
            TGraph gleg{};
            gleg.SetTitle("");
            gleg.SetMarkerStyle(8);
            gleg.SetMarkerSize(1);

            TLegend bleg (0.1,0.2,0.9,0.5);
            bleg.SetBorderSize(0);
            bleg.SetTextSize(18);
            bleg.AddEntry(&gleg,TString::Format("#frac{#sigma(%s)}{#sigma(SM)}",fTitle.c_str()),"p");
            bleg.AddEntry(&hleg,"Fit","l");
            bleg.Draw();

            fExtrapMap[SMRef]->nBins=SM_hist->GetNbinsX();
            double bottom_margin=0.3;


            // Build TGraph per Region bin and fit it
            for(int b=0; b<SM_hist->GetNbinsX(); b++){
                b1.cd(b+2);

                std::vector<double> x;
                std::vector<double> x_err;
                std::vector<double> y;
                std::vector<double> y_err;

                const std::string NFname = Form("Expression_muEFT_%s_%s_bin%d_%s",imap.first.c_str(),rname.c_str(),b,fName.c_str());

                // Skip empty bins
                if(SM_hist->GetBinContent(b+1) <= 1e-6){
                    fExtrapMap[SMRef]->p0s.emplace_back(1.);
                    fExtrapMap[SMRef]->p1s.emplace_back(0);
                    fExtrapMap[SMRef]->p2s.emplace_back(0);
                    continue;
                }

                // Set SM point
                x.emplace_back(0.);
                x_err.emplace_back(0.);
                y.emplace_back(1.);
                if (SM_hist->GetBinContent(b+1) != 0) {
                    y_err.emplace_back(SM_hist->GetBinError(b+1)/SM_hist->GetBinContent(b+1));
                } else {
                    y_err.emplace_back(9999999);
                }

                // Cycle through all EFT variations to get yields and extract relative changes
                for(const auto& EFT_name : imap.second){

                    std::shared_ptr<SampleHist> EFT_sh = ireg->GetSampleHist(EFT_name);

                    std::unique_ptr<TH1> EFT_hist(static_cast<TH1*>(EFT_sh->fHist->Clone()));
                    EFT_hist->SetDirectory(nullptr);

                    if(EFT_sh->GetSample()->fEFTValue>0)EFT_hist->SetLineStyle(kDashed);
                    else EFT_hist->SetLineStyle(kDotted);

                    double rel_eft=EFT_hist->GetBinContent(b+1)/SM_hist->GetBinContent(b+1);
                    double rel_eft_err=rel_eft*std::sqrt(((EFT_hist->GetBinError(b+1)/EFT_hist->GetBinContent(b+1))*(EFT_hist->GetBinError(b+1)/EFT_hist->GetBinContent(b+1))) +
                      				 ((SM_hist->GetBinError(b+1)/SM_hist->GetBinContent(b+1))*(SM_hist->GetBinError(b+1)/SM_hist->GetBinContent(b+1))) );

                    if(rel_eft > max_rel_eft) max_rel_eft = rel_eft;


                    if(EFT_sh->GetSample()->fEFTValue > max_c) max_c=EFT_sh->GetSample()->fEFTValue;
                    if(EFT_sh->GetSample()->fEFTValue < min_c) min_c=EFT_sh->GetSample()->fEFTValue;

                    x.emplace_back(EFT_sh->GetSample()->fEFTValue);
                    x_err.emplace_back(0.);
                    y.emplace_back(rel_eft);
                    y_err.emplace_back(rel_eft_err);
                }

                //Plot graph
                TGraphErrors g(imap.second.size()+1, &(x[0]), &(y[0]), &(x_err[0]), &(y_err[0]));
                g.SetTitle("");
                g.Fit("pol2");
                TF1* fit = g.GetFunction("pol2");
                if (!fit) {
                    WriteErrorStatus("EFTProcessor::FitEFTInputs", "Cannot retrieve function!");
                    continue;
                }
                Double_t p0 = fit->GetParameter(0);
                Double_t p1 = fit->GetParameter(1);
                Double_t p2 = fit->GetParameter(2);

                // Write results to file
                NFTXTFile << NFname << "  " << p0 << " " << p1 << " " << p2 << "\n";

                // Save fitter quadratic coeffs to be added to NF later
                fExtrapMap[SMRef]->p0s.push_back(p0);
                fExtrapMap[SMRef]->p1s.push_back(p1);
                fExtrapMap[SMRef]->p2s.push_back(p2);

                fit->SetLineColor(kRed);
                fit->SetLineWidth(2);

                g.SetMarkerStyle(8);
                g.SetMarkerSize(1);

                // Plot the fit for this bin
                TH1 *frame = gPad->DrawFrame(min_c-(0.2*std::abs(min_c)),0.8,max_c+(0.2*std::abs(max_c)),1.1*max_rel_eft,TString::Format("frame%d",b+2));
                frame->GetXaxis()->SetTitle(fTitle.c_str());
                bottom_margin=gPad->GetBottomMargin();

                g.Draw("P");

                // Also draw the fitted equation
                TLatex gtex{};
                gtex.SetNDC();
                gtex.SetTextSize(14);
                gtex.DrawLatex(0.17,0.90,TString::Format(ireg->fVariableTitle.c_str()));
                gtex.DrawLatex(0.17,0.85,TString::Format("Bin%d",b));
                gtex.SetTextSize(13);
                gtex.DrawLatex(0.17,0.80,TString::Format("y=%.2fx^{2} +",p2));
                gtex.DrawLatex(0.17,0.76,TString::Format("    %.2fx +",p1));
                gtex.DrawLatex(0.17,0.72,TString::Format("      %.2f",p0));


            }

            // Need this hack to get y-axis shown in a nice way
            b1.cd(1);
            gPad->Range(-10,-1,10,1);
            TGaxis yaxis(10,(-1+2*bottom_margin),10,1.,0.8,1.1*max_rel_eft,512,"");
            yaxis.SetTitle(TString::Format("#sigma(%s)/#sigma(SM)",fTitle.c_str()));
            yaxis.SetTitleSize(0.08);
            yaxis.SetLabelSize(0.07);
            yaxis.Draw();

            // Go back and fix all y-axis ranges to same value
            for(int ip=2;ip<=SM_hist->GetNbinsX()+1;ip++){
                b1.cd(ip);
                TH1* tmp = static_cast<TH1*>(gPad->GetPrimitive("hframe"));
                tmp->SetMaximum(1.1*max_rel_eft);
                gPad->RedrawAxis();
            }

            for(const auto& format : TRExFitter::IMAGEFORMAT) {
                b1.SaveAs(TString::Format("%s/EFT/QuadFits_%s_%s_%s.%s",SM_sh->fFitName.c_str(),rname.c_str(),fName.c_str(),SMRef.c_str(),format.c_str()));
            }
        }// SM ref loop
    }// region loop

    NFTXTFile << "\n";
    NFTXTFile.close();
}


//_____________________________________________________________________________
// Apply mu NormFactors to each bin based on quadratic fit output.
void EFTProcessor::ApplyMuFactExpressions(std::vector < Region* > Regions,
                                          std::vector < std::shared_ptr<NormFactor> > &NormFactors){

    for(auto& ireg : Regions) {
        const std::string rname=ireg->fName;

        for (const auto& imap : fRefMap) {
            // Get SM Reference sample
            std::shared_ptr<SampleHist> SM_sh = ireg->GetSampleHist(imap.first);
            if( !SM_sh ) continue;

            for(int ib=0; ib < ireg->fNbins; ib++){

                // Form unique name for this region-bin-sample combination
                const std::string NFname = Form("Expression_muEFT_%s_%s_bin%d_%s",imap.first.c_str(),rname.c_str(),ib,fName.c_str());

                // Add the expression based on the result of the quadratic fit
                const std::string quad_eqn = TString::Format("((%f*%s*%s)+(%f*%s)+%f)",fExtrapMap[imap.first]->p2s[ib],fName.c_str(),fName.c_str(),fExtrapMap[imap.first]->p1s[ib],fName.c_str(),fExtrapMap[imap.first]->p0s[ib]).Data();
                const std::string dependancies = TString::Format("%s[%f,%f,%f]",fName.c_str(),0.0,-100.0,100.0).Data();

                WriteDebugStatus("EFTProcessor::ApplyMuFactExpressions", "    Defining EFT mu NormFactor Expression for: " + NFname);
                WriteDebugStatus("EFTProcessor::ApplyMuFactExpressions", "             " + quad_eqn + ":" + dependancies);

                // Add to matching NF
                for(auto norm : NormFactors){
                    if(norm->fName!=NFname) continue;
                    WriteDebugStatus("EFTProcessor::ApplyMuFactExpressions", "      --> Updating NF " + NFname);
                    norm->fExpression = std::make_pair(quad_eqn,dependancies);
                    //Need to name like this because formula is used later on from fTitle...
                    const std::string newname = Form("Expression_%s",quad_eqn.c_str());
                    norm->fTitle = newname;
                    TRExFitter::SYSTMAP[norm->fName] = newname;
                    // nuis-par will contain the nuis-par of the norm factor the expression depends on FIXME
                    const std::string newnpname = Form("Expression_%s",dependancies.c_str());
                    norm->fNuisanceParameter = newnpname;
                    TRExFitter::NPMAP[norm->fName] = newnpname;

                    WriteDebugStatus("EFTProcessor::ApplyMuFactExpressions", "         --> Now NF " + norm->fTitle);
                    break;
                }
            }
        }
    }
}


//_____________________________________________________________________________
// Read quadratic fit results from file and pass them to relevant EFTProcessor objects
void EFTProcessor::ReadEFTFitResults(std::vector < Region* > Regions, const std::string& fileName) {

    WriteDebugStatus("EFTProcessor::ReadEFTFitResults"," Opening file \"" + fileName + "\"");

    std::ifstream NFTXTFile;
    NFTXTFile.open(fileName.c_str());

    if (!NFTXTFile.is_open()) {
        WriteErrorStatus("EFTProcessor::ReadEFTFitResults","Could not open the file \"" + fileName + "\"");
        return;
    }


    std::string input;
    std::string line;
    bool readingNF = false;
    //
    std::string name;
    double p0,p1,p2;

    //
    // read file line by line
    while(std::getline(NFTXTFile, line)){
        if(line=="") continue;
        if(line=="NORMFACTORS"){
            WriteDebugStatus("EFTProcessor::ReadEFTFitResults", "--------------------");
            WriteDebugStatus("EFTProcessor::ReadEFTFitResults", "Reading Norm Factors...");
            WriteDebugStatus("EFTProcessor::ReadEFTFitResults", "--------------------");
            readingNF = true;
            continue;
        }
        std::istringstream iss(line);
        if(readingNF){
            iss >> input;
            if(input==""){ //leaving NF block
                readingNF = false;
            }
            while(input.find("\\")!=std::string::npos) input = input.replace(input.find("\\"),1,"");
            name = input;
	        iss >> p0 >> p1 >> p2;

            WriteVerboseStatus("EFTProcessor::ReadEFTFitResults", name + ": p0=" + std::to_string(p0) + " p1=" + std::to_string(p1) + " p2=" + std::to_string(p2));

            //Find appropriate SM Ref
            bool foundEntry=false;
            for(const auto& ireg : Regions) {
                const std::string rname=ireg->fName;

                for (const auto& imap : fRefMap) {

                    // Get SM Reference sample
                    std::shared_ptr<SampleHist> SM_sh = ireg->GetSampleHist(imap.first);
                    if( !SM_sh )continue;

                    const std::string SMRef = imap.first;

                    for(int b=0; b<ireg->fNbins; b++){
                        const std::string NFname=Form("Expression_muEFT_%s_%s_bin%d_%s",imap.first.c_str(),rname.c_str(),b,fName.c_str());
                        if ( name != NFname ) continue;
                        WriteVerboseStatus("EFTProcessor::ReadEFTFitResults", "  -> Found matching ExtrapMap");

                        fExtrapMap[SMRef]->p0s.push_back(p0);
                        fExtrapMap[SMRef]->p1s.push_back(p1);
                        fExtrapMap[SMRef]->p2s.push_back(p2);

                        foundEntry=true;
                        break;
                    }
    	            if(foundEntry) break;
                }
                if(foundEntry) break;
            }
        }
    }

    NFTXTFile.close();
}
