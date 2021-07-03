import sys, os, argparse
import numpy as np
import ROOT
from ROOT import vector
ENV_WORK_DIR=os.environ["WorkDir_DIR"]
ROOT.gSystem.Load("%s/lib/libCommonSystSmoothingToolLib.so"%ENV_WORK_DIR)
ROOT.gStyle.SetOptStat(0)

parse = argparse.ArgumentParser(description='Smooth systematics hists from root file.', prog='testInputFile.py', add_help=True)

parse.add_argument('--inputFile', default='13TeV_TwoLepton_2tag2pjet_150ptv_SR_mvaTEST.root',type=str, help='Input root file containing nominal and systematic histograms')
parse.add_argument('--outputFile', default='13TeV_TwoLepton_2tag2pjet_150ptv_SR_mvaSmoothed.root',type=str, help='Target file: nominal, systematic and smoothed histograms will be saved.')

parse.add_argument('--sysTdir', default='Systematics', type=str, help='TDirectory name of systematics')
parse.add_argument('--nomTdir', default='', type=str, help='TDirectory name of nominal hists.')
parse.add_argument('--smoothOption', default='smoothRebinParabolic', type=str, choices=['smoothRebinParabolic', 'smoothRebinMonotonic', 'smoothRatioUniformKernel','smoothDeltaUniformKernel','smoothRatioGaussKernel','smoothDeltaGaussKernel','smoothTtresDependent','smoothTtresIndependent','smoothTRExDefault'], help='Preconfigured methods')

parse.add_argument('--rebin', default=20, type=int, help='rebin factor')
parse.add_argument('--generatelist',default=10, type=int, help='number of histograms to be smoothed, give -1 to smooth all hists.')
parse.add_argument('--doPlot',action='store_true', default=False)

parser=parse.parse_args()

ST = ROOT.SmoothingTool()
ST.addInputFileNom(parser.inputFile, parser.nomTdir)
ST.addInputFileSys(parser.inputFile, parser.sysTdir)

ST.addOutputFileSys(parser.outputFile)

ST.setSmoothingOption(parser.smoothOption)
ST.setRebinFactor(parser.rebin)

skipVector=ROOT.vector('string')()

skiplist=["data",   "ggZll",     "stop",     "Wcc",   "Wl",  "Wcl","Zl",     "Zcc",       "Zcl",      "ttbar", "ttH", "Wbc","Wb",     "Wbb",       "WlvH125",  "WZ",    "Zbc", "Zbl","MadZee", "qqZllH125", "MadZmumu", "ZZ","SysEL_EFF_Reco_TOTAyL_1NPCOR_PLUS_UNCOR", "SysEL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR","SysEL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR", "SysEL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR","SysMUON_EFF_SYS","SysMUON_EFF_STAT", "SysMUON_ISO_SYS","SysMUON_ISO_STAT","SysMUON_TTVA_SYS","SysMUON_TTVA_STAT","SysFT_EFF_extrapolation_from_charm","SysFT_EFF_extrapolation","SysFT_EFF_Eigen_Light_0","SysFT_EFF_Eigen_Light_1","SysFT_EFF_Eigen_Light_2","SysFT_EFF_Eigen_Light_3","SysFT_EFF_Eigen_Light_4","SysFT_EFF_Eigen_Light_5","SysFT_EFF_Eigen_Light_6","SysFT_EFF_Eigen_Light_7","SysFT_EFF_Eigen_Light_8","SysFT_EFF_Eigen_Light_9","SysFT_EFF_Eigen_Light_10","SysFT_EFF_Eigen_Light_11","SysFT_EFF_Eigen_Light_12","SysFT_EFF_Eigen_Light_13","SysFT_EFF_Eigen_C_0","SysFT_EFF_Eigen_C_1","SysFT_EFF_Eigen_C_2","SysFT_EFF_Eigen_C_3","SysFT_EFF_Eigen_B_0","SysFT_EFF_Eigen_B_1","SysFT_EFF_Eigen_B_2","SysFT_EFF_Eigen_B_3","SysFT_EFF_Eigen_B_4","SysFT_EFF_Eigen_B_5","SysVVMbbME","SysVVPTVME","SysVHQCDscalePTV","SysVHQCDscaleMbb","SysVHPDFPTV","SysVHQCDscalePTV_ggZH","SysVHQCDscaleMbb_ggZH","SysVHPDFPTV_ggZH","SysVHUEPSPTV","SysVHUEPSMbb","SysVHNLOEWK","SysStoptPTV","SysStoptMBB","SysStopWtPTV","SysStopWtMBB"]

for iSkip in skiplist:
  skipVector.push_back(iSkip)



ST.skipHists(skipVector)



ST.generateList(parser.generatelist)
ST.makePlot(parser.doPlot)
ST.setRebinFactorPlot(-1)
ST.runSmoothing()

