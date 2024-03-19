import ROOT
from functions import getRootFiles

#input reading
input = open("test/comparison-input.txt","r")
inputLines = input.readlines()
inputLines = [line.rstrip() for line in inputLines]

# file calling
normalDfile = ROOT.TFile(inputLines[1])
norecoDfile = ROOT.TFile(inputLines[3])
histlist = getRootFiles(normalDfile)

# looping over individual plots
for hist in histlist:
    # skipping plots for other stuff
    if " vs " in hist:
        continue

    normalHist = normalDfile.Get(hist)
    normalHist.SetName(inputLines[0])
    normalHist.SetLineColor(ROOT.kRed)

    norecoHist = norecoDfile.Get(hist)
    norecoHist.SetName(inputLines[2])

    title = hist+";"+normalHist.GetXaxis().GetTitle()+";"+normalHist.GetYaxis().GetTitle()
    stack = ROOT.THStack(hist,title)
    stack.Add(norecoHist)
    stack.Add(normalHist)

    legend = ROOT.TLegend(0.8 ,0.8, 0.9, 0.9)
    legend.AddEntry(normalHist, inputLines[0])
    legend.AddEntry(norecoHist, inputLines[2])

    c1 = ROOT.TCanvas()
    stack.Draw("nostack")
    legend.Draw("same")

    if "Spectrum" in hist:
        c1.SetLogy()
        c1.SetLogx()
    c1.Print(inputLines[4]+hist+".png")