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
    if "Spectrum" not in hist:
        continue

    normalHist = normalDfile.Get(hist)
    normalHist.SetName(inputLines[0])
    normalHist.SetLineColor(ROOT.kRed)
    normalHist.SetStats(0)

    norecoHist = norecoDfile.Get(hist)
    norecoHist.SetName(inputLines[2])
    norecoHist.SetStats(0)

    title = hist+";"+normalHist.GetXaxis().GetTitle()+";"+normalHist.GetYaxis().GetTitle()
    #stack = ROOT.THStack(hist,title)
    #stack.Add(norecoHist)
    #stack.Add(normalHist)

    legend = ROOT.TLegend(0.8,0.8,0.9,0.9)
    legend.AddEntry(normalHist, inputLines[0])
    legend.AddEntry(norecoHist, inputLines[2])

    c1 = ROOT.TCanvas()
    if "Spectrum" in hist:
        c1.SetLogy()
        c1.SetLogx()
    #stack.Draw("nostack")
    rp = ROOT.TRatioPlot(norecoHist,normalHist)
    rp.Draw()
    legend.Draw("same")
    c1.Print(inputLines[4]+hist+".png")