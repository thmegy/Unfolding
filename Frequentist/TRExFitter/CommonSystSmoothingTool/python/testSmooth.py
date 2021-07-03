
import sys,os,argparse
import ROOT

ENV_WORK_DIR=os.environ["WorkDir_DIR"]

ROOT.gSystem.Load("%s/lib/libCommonSystSmoothingToolLib.so"%ENV_WORK_DIR)

ROOT.gStyle.SetOptStat(0)

parse = argparse.ArgumentParser(description='Test Smoothing methods.',prog='testSmooth.py', add_help=True)
parse.add_argument('--inputFile', dest='inputFile', default='13TeV_TwoLepton_2tag2pjet_150ptv_SR_mvaTEST.root',type=str, help='Input root file.')
parse.add_argument('--hnom', dest='hnom', default='Zbb',type=str, help='Name of nominal histogram')
parse.add_argument('--hsys', dest='hsys',default='Zbb_SysJET_21NP_JET_EffectiveNP_2__1up', type=str, help='Name of systematic histogram')
parse.add_argument('--methods', dest='methods', type=str, nargs='+', default=['smoothRebinMonotonic', 'smoothRebinParabolic', 'smoothDeltaUniformKernel','smoothTtresDependent','smoothTRExDefault'])

parse.add_argument('--output', dest='output', type=str, default='testPlot.pdf', help="Name of plot <name.pdf>")


parser=parse.parse_args()
inFile = ROOT.TFile.Open(parser.inputFile)
hnom = inFile.Get(parser.hnom)
hnom.SetLineColor(ROOT.kBlue)

tdSys = inFile.Get("Systematics")
hsys = tdSys.Get(parser.hsys)
hsys.SetLineColor(ROOT.kRed)

RebinFactor=20
hnom.Rebin(RebinFactor)
hsys.Rebin(RebinFactor)

smoothTool=ROOT.SmoothHist()
plotTool = ROOT.PlotHist()

smoothTool.setStatErrThreshold(0.05)
#smoothTool.setSpans([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8,
#                                 1.0, 1.2, 1.4, 1.6, 1.8, 2.0])

smoothTool.setTRExTolerance(0.08)
smoothTool.setTRExNbins(1)

class Method:
  def __init__(self, Name, LineStyle, LineColor, RatioLineStyle):
    self.Name=Name
    self.LineStyle=LineStyle
    self.LineColor=LineColor
    self.RatioLineStyle=RatioLineStyle

  def hStyle(self,hist):
    hist.SetLineColor(self.LineColor)
    hist.SetLineStyle(self.LineStyle)
    hist.SetTitle("h%s"%self.Name)
  def hRatioStyle(self,hist):
    hist.SetLineStyle(self.RatioLineStyle)
    hist.SetLineColor(self.LineColor)
    hist.SetLineWidth(3)


MethodCollection={"smoothRebinParabolic": Method("smoothRebinParabolic",10,30,7),
    "smoothRebinMonotonic":Method("smoothRebinMonotonic",10,31,7),
    "smoothRatioUniformKernel":Method("smoothRebinMonotoni",10,ROOT.kGreen+2,7),
    "smoothDeltaUniformKernel":Method("smoothDeltaUniformKernel",10,ROOT.kGreen+4,7),
    "smoothRatioGaussKernel":Method("smoothRatioGaussKernel",10,ROOT.kGreen+6,7),
    "smoothDeltaGaussKernel": Method("smoothDeltaGaussKernel", 10,ROOT.kGreen,7),
    "smoothTtresDependent": Method("smoothTtresDependent",10,ROOT.kMagenta+2,7),
    "smoothTtresIndependent": Method("smoothTtresIndependent",10,ROOT.kMagenta+4,7),
    "smoothTRExDefault": Method("smoothTRExDefault",10,ROOT.kYellow+2,7)}

#plot
canv = ROOT.TCanvas("canvTest","Nom/Sys",800,800)
pad1 = ROOT.TPad("pad1Test", "pad1", 0, 0.5, 1, 1.0)
pad1.SetBottomMargin(0)
pad1.SetTicky(2)
pad1.Draw()
pad1.cd()

hnom.Draw("Hist")
hsys.Draw("Hist same")

leg = ROOT.TLegend(0.6,0.3,0.9,0.9)
leg.AddEntry(hnom, "nominal","l")
leg.AddEntry(hsys,hsys.GetName(), "l")

hRlist=[]
for iMethod in parser.methods:
  tmpHist=smoothTool.Smooth(hnom, hsys.Clone("h%s"%iMethod),iMethod)
  MethodCollection[iMethod].hStyle(tmpHist)
  tmpHist.Draw("Hist same")
  leg.AddEntry(tmpHist, iMethod,"l")

  tmpHistR=plotTool.getPull(hnom, tmpHist)
  MethodCollection[iMethod].hRatioStyle(tmpHistR)
  hRlist.append(tmpHistR)

leg.Draw()

canv.cd()
pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.498)
pad2.SetTopMargin(0)
pad2.SetBottomMargin(0.2)
pad2.SetGridx() 
pad2.SetGridy()
pad2.SetTicky(2)
pad2.SetTickx(2)
pad2.Draw()
pad2.cd()

hnomR = plotTool.getPull(hnom,hnom, True)
hnomR.SetFillStyle(3001)
hnomR.SetFillColor(ROOT.kGray)
hnomR.GetXaxis().SetTitle(hnom.GetXaxis().GetTitle())
hsysR = plotTool.getPull(hnom,hsys)
minR = 1.3*abs(hsysR.GetMinimum())
maxR = 1.3*abs(hsysR.GetMaximum())
maxR = max(minR, maxR)
if maxR>10.: maxR=10
if maxR>1.: 
  hnomR.SetMinimum(-maxR)
  hnomR.SetMaximum(maxR)
hnomR.SetBit(ROOT.TH1.kNoTitle, True)
hnomR.SetLineColor(ROOT.kBlack)
hnomR.GetYaxis().SetLabelSize(15)
hnomR.GetYaxis().SetLabelFont(43)
hnomR.GetYaxis().SetTitleFont(43)
hnomR.GetYaxis().SetTitleSize(20)
hnomR.GetYaxis().SetTitleOffset(1.55)
hnomR.GetYaxis().SetTitle("#frac{Sys-Nom}{Nom} [%]")
hnomR.GetXaxis().SetLabelSize(15)
hnomR.GetXaxis().SetLabelFont(43)
hnomR.GetXaxis().SetTitleFont(43)
hnomR.GetXaxis().SetTitleSize(25)
hnomR.GetXaxis().SetTitleOffset(4.0)
hnomR.Draw("hist e2")

hsysR.SetLineColor(ROOT.kRed)
hsysR.SetLineStyle(2)
hsysR.SetLineWidth(3)
hsysR.Draw("hist same")


for ihR in hRlist:
  ihR.Draw("hist same")

canv.Print("%s"%parser.output)




