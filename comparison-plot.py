import ROOT
from functions import getRootFiles

# file calling
normalDfile = ROOT.TFile("test/Dtotals.root")
norecoDfile = ROOT.TFile("runs/Dnoreco/points/0/root/totals.root")
histlist = getRootFiles(normalDfile)

# looping over individual plots
for hist in histlist:
    # skipping plots for other stuff
    if "GeV" not in hist:
        continue

    normalHist = normalDfile.Get(hist)
    normalHist.SetName("Normal")
    normalHist.SetLineColor(ROOT.kRed)

    norecoHist = norecoDfile.Get(hist)
    norecoHist.SetName("No reco")

    stack = ROOT.THStack(hist,hist+";\Delta\phi;1/N_{D}dN_{pairs}/d\Delta\phi")
    stack.Add(normalHist)
    stack.Add(norecoHist)

    legend = ROOT.TLegend(0.8 ,0.8, 0.9, 0.9)
    legend.AddEntry(normalHist, "Total" )
    legend.AddEntry(norecoHist, "No Reco" )

    c1 = ROOT.TCanvas()
    stack.Draw("nostack")
    legend.Draw("same")
    c1.Print("test/Drecocomparison/"+hist+".png")