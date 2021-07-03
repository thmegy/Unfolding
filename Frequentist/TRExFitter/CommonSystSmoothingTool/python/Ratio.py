import ROOT
import sys


def openRootFile(fName, option="read"):
  oFile = ROOT.TFile.Open(fName, option)
  if oFile == 0:
    if option == "read":
      print("ERROR: Cannot open file: %s"%fName)
    else: 
      print("ERROR: Cannot create file: %s"%fName)
    sys.exit(0)
  return oFile


def getHist(source, hName):
  source.cd()
  tmpH = source.Get(hName)
  if (tmpH.InheritsFrom(ROOT.TH1.Class())):
    return tmpH
  else:
    print("ERROR: %s is not a TH1."%hName)
    sys.exit(0)


def plotRatio(hbase, htarget):
  ROOT.gStyle.SetOptStat(0)
  c = ROOT.TCanvas("c","",800,200)
  htarget.Divide(hbase);
  htarget.SetMarkerStyle(20)
  htarget.SetMaximum(1.5)
  htarget.SetMinimum(0.5)
  htarget.Draw("p")
  c.SetGridy(1);
  return c


def main():

  f1Name = "13TeV_TwoLepton_2tag2pjet_150ptv_SR_mvaSmoothed.root"
  f2Name = "13TeV_TwoLepton_2tag2pjet_150ptv_SR_mvaTEST_smoothedRatio.root"

  file1 = openRootFile(f1Name)
  tdir1 = file1.Get("Systematics")
  file2 = openRootFile(f2Name)
  tdir2 = file2.Get("Systematics")

  h1 = getHist(tdir1, "Zbb_SysJET_21NP_JET_EffectiveNP_2__1upSmooth")
  h2 = getHist(tdir2, "Zbb_SysJET_21NP_JET_EffectiveNP_2__1up")

  canv = plotRatio(h1, h2)
  canv.Print("tmpRatio.pdf")

if __name__== "__main__":
  main()

