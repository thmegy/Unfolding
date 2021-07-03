#include "SmoothSystematics/SmoothHist.h"
#include "SmoothSystematics/SmoothingTool.h"
#include "TStyle.h"
#include "TColor.h"
#include "TH1F.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TFile.h"

#include "TLegend.h"
#include "TLegendEntry.h"

#include <iostream>

TH1* getPull(TH1* hNOM, TH1* hSYS, bool fillErr=false){
  string sName = hSYS->GetTitle();
  sName += "hPull";
  TH1F* hPull =(TH1F*)hNOM->Clone(sName.data()); //clone empty hist.
  hPull->Reset("ICESM"); 
  int nBins = hNOM->GetNbinsX();
  for (int iBin=0; iBin<=nBins+1; iBin++){
    float bConNom = hNOM->GetBinContent(iBin);
    float bConSys = hSYS->GetBinContent(iBin);
    if(bConNom == 0 ){
      hPull->SetBinContent(iBin,0.);
    }else{
      hPull->SetBinContent(iBin,(bConSys-bConNom)*100./bConNom);
    }
    if(fillErr){
      float bErrNom = hNOM->GetBinError(iBin);
      if(bConNom == 0 ){
        hPull->SetBinError(iBin,0.);
      }else{
        hPull->SetBinError(iBin,bErrNom*100./bConNom);
      }
    }else{
      hPull->SetBinError(iBin,0);
    }
  }
  return hPull;
}


                                                      
int main(){
  gStyle->SetOptStat(0); //Turn off stat box.
  
  //Open input file.
  TFile* inFile = TFile::Open("13TeV_TwoLepton_2tag2pjet_150ptv_SR_mvaTEST.root");

  //Read nominal histogram to memory.
  TH1F* hnom = (TH1F*)inFile->Get("Zbb");
  hnom->SetLineColor(kBlue); 
 
  /*Get systematics histogram to mamory. 
   * Sys hists are in a TDirectory.*/
  TDirectory* tdSys = (TDirectory*) inFile->Get("Systematics");
  TH1F* hsys = (TH1F*)tdSys->Get("Zbb_SysJET_21NP_JET_EffectiveNP_2__1up");
  hsys->SetLineColor(kRed);
  
  int RebinFactor =20; //= hnom->GetNbinsX() / 10;
  hnom->Rebin(RebinFactor);
  hsys->Rebin(RebinFactor);
  

  //Creat Smoothing Alg
  SmoothHist smoothTool; 

  /************* This method is used in RF*********************
   * This method will merge bins based on stat>0.05.
   * Better to clone sys hist. otherwise original hist might be modified.
   * Below is the example of configuring smoothTool for this method.*/

  /*Set max number of maximum extrema.
   * Here it set it to 1. This means total number of
   * extrema will be 3 (1 max extrema + 2 min extrema).
   * */
  //smoothTool.setNmax(1); //set maximum number of extrema
  
  /*Set threshold for stat error.
   * if stat error of bins are larger then this threshold
   * merge bins between two extrema untill number of extrema reach
   * the value set by setNmax(N).
   * If total stat error of nominal hist larger then the threshold,
   * skip smoothing will appy.
   * */
  smoothTool.setStatErrThreshold(0.05);
  
  //Apply smoothing
  TH1* hsmooth = smoothTool.Smooth(hnom, (TH1*)hsys->Clone("hsmooth"), "smoothRebinParabolic");
  hsmooth->SetTitle("hsmoothRebinParabolic"); 
  //Set styles for smoothed hist.
  hsmooth->SetLineColor(30);
  hsmooth->SetLineStyle(10);

  /************** This method is used in WSM *******************
   * Alternative method: use kernel fuction to smooth histograms.
   * Better to clone sys hist. otherwise original hist might be modified.
   * Below is the example configuring smoothTool for this method.
   * */
   
   /* There are two kernel fuction available: "box" and "normal".
   * "box"-> take average of several bins.
   * "normal"-> Gaussian fuction.
   * */
  
  /*There are two ways to apply smoothing:
   * 1) smoothe the ratio:  Sys/nominal.
   * 2) Smooth the delta: Sys-Nominal
   * */

  /*This is a list of radii for kernel fuction.
   * If you dont know what to set then dont call this fuction.
   * The default list will be used.*/
  smoothTool.setSpans({0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8,
                                 1.0, 1.2, 1.4, 1.6, 1.8, 2.0
                                 });
  
  /*Rebin histograms befor smoothing.
   * Not sure if this is needed.
   * -1 means dont rebin.*/
  //smoothTool.setRebinFactor(20);
  
  //Apply kernel based smoothing.
  TH1* hsmoothKernel = smoothTool.Smooth(hnom, (TH1*)hsys->Clone("hsmoothRatioUniformKernel"),"smoothRatioUniformKernel");
  
  hsmoothKernel->SetLineColor(kGreen+2);
  hsmoothKernel->SetLineStyle(10);
  hsmoothKernel->SetTitle("hsmoothRatioUniformKernel");
  
  //Appy TRExFitter Smoothing
  //smoothTool.setIndependentVar(false);
  TH1* hsmoothTtres = smoothTool.Smooth(hnom, (TH1*)hsys->Clone("hsmoothTtresDependent"),"smoothTtresDependent");
  hsmoothTtres->SetLineColor(kMagenta+2);
  hsmoothTtres->SetLineStyle(10);
  hsmoothTtres->SetTitle("hsmoothTtresDependent");
  
  //TREx default smoothing.
  //smoothTool.setTRExTolerance(0.08);
  //smoothTool.setTRExNbins(1);
  TH1* hsmoothTRExDefault = smoothTool.Smooth(hnom, (TH1*)hsys->Clone("hsmoothTRExDefault"), "smoothTRExDefault");
  hsmoothTRExDefault->SetLineColor(kYellow+2);
  hsmoothTRExDefault->SetLineStyle(10);
  hsmoothTRExDefault->SetTitle("hsmoothTRExDefault");
  
  /*******************
   * Plot histograms *
   * *************** *
   */
  TCanvas* canv = new TCanvas("canvTest","Nom/Sys",800,800);
  TPad *pad1 = new TPad("pad1Test", "pad1", 0, 0.5, 1, 1.0);
  pad1->SetBottomMargin(0);
  pad1->SetTicky(2);
  pad1->Draw();
  pad1->cd();
  
  //TLegend* leg = new TLegend(0.6,0.4,0.9,0.9);
  hnom->Draw("Hist");
  hsys->Draw("Hist same");
  hsmooth->Draw("hist same");
  hsmoothKernel->Draw("hist same");
  hsmoothTtres->Draw("hist same");
  hsmoothTRExDefault->Draw("hist same");
  //leg->AddEntry(hnom, "nominal","l");
  //leg->AddEntry(hsys,"systematic", "l");
  //leg->AddEntry(hsmooth, "user defined smoothing", "l");
  //leg->AddEntry(hsmoothKernel, "kernel based smoothing", "l");
  //leg->Draw();
  pad1->BuildLegend(0.6,0.4,0.9,0.9);

  canv->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.498);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); 
  pad2->SetGridy();
  pad2->SetTicky(2);
  pad2->SetTickx(2);
  pad2->Draw();
  pad2->cd();
  TH1* hnomR = getPull(hnom,hnom, true);
  hnomR->SetFillStyle(3001);
  hnomR->SetFillColor(kGray);
  hnomR->GetXaxis()->SetTitle(hnom->GetXaxis()->GetTitle());
  TH1 *hsysR = getPull(hnom,hsys); //(TH1 *)hsys->Clone("Sys");
  float minR = 1.3*abs(hsysR->GetMinimum());
  float maxR = 1.3*abs(hsysR->GetMaximum());
  maxR = std::max(minR, maxR);
  if (maxR>10.){
    maxR = 10.0;
  }
  if(maxR>1.){
    hnomR->SetMinimum(-maxR);  // Define Y ..
    hnomR->SetMaximum(maxR); // .. range
  }
  hnomR->SetBit(TH1::kNoTitle, true);
  hnomR->SetLineColor(kBlack);
  hnomR->GetYaxis()->SetLabelSize(15);
  hnomR->GetYaxis()->SetLabelFont(43);
  hnomR->GetYaxis()->SetTitleFont(43);
  hnomR->GetYaxis()->SetTitleSize(20);
  hnomR->GetYaxis()->SetTitleOffset(1.55);
  hnomR->GetYaxis()->SetTitle("#frac{Sys-Nom}{Nom} [%]");
  hnomR->GetXaxis()->SetLabelSize(15);
  hnomR->GetXaxis()->SetLabelFont(43);
  hnomR->GetXaxis()->SetTitleFont(43);
  hnomR->GetXaxis()->SetTitleSize(25);
  hnomR->GetXaxis()->SetTitleOffset(4.0);
  hnomR->Draw("hist e2");

  hsysR->SetLineColor(kRed);
  hsysR->SetLineStyle(2);
  hsysR->SetLineWidth(3);
  hsysR->Draw("hist same");
  TH1 *hsmoothR = getPull(hnom,hsmooth);
  hsmoothR->SetLineColor(30);
  hsmoothR->SetLineStyle(7);
  hsmoothR->SetLineWidth(3);
  hsmoothR->Draw("hist same");

  TH1 *hsmoothKernelR = getPull(hnom,hsmoothKernel);
  hsmoothKernelR->SetLineColor(kGreen+2);
  hsmoothKernelR->SetLineStyle(7);
  hsmoothKernelR->SetLineWidth(3);
  hsmoothKernelR->Draw("hist same");
  
  TH1 *hsmoothTtresR = getPull(hnom,hsmoothTtres);
  hsmoothTtresR->SetLineColor(kMagenta+2);
  hsmoothTtresR->SetLineStyle(7);
  hsmoothTtresR->SetLineWidth(3);
  hsmoothTtresR->Draw("hist same");
  
  TH1* hsmoothTRExDefaultR = getPull(hnom, hsmoothTRExDefault);
  hsmoothTRExDefaultR->SetLineColor(kYellow+2);
  hsmoothTRExDefaultR->SetLineStyle(7);
  hsmoothTRExDefaultR->SetLineWidth(3);
  hsmoothTRExDefaultR->Draw("hist same");


  canv->Print("testPlot.pdf");

  delete canv;
  //delete linear;
}

