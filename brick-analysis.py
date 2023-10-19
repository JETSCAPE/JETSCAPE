# makes set of brick systems for comparing photon absoroption
from multiprocessing import connection
import os
from datetime import date
from functions import *
import pandas as pd
from operator import itemgetter
import xml.etree.ElementTree as ET
import ROOT
from functions import parton

# setting directory for analysis
analysisDir = sys.argv[1]
directories = getDirs(analysisDir)
outdir = analysisDir+"analysis/"
outfile = ROOT.TFile.Open(outdir+'analysis.root', "RECREATE")

# initializing histograms
#pTs = [1,3,5,7,9]
pTs = [0.5,1.5]
lengths = [1,3,5,7,9,11,13]
totalhist = ROOT.TH2D('Photon Rate','Photon Rate; Brick Length (fm); pT (GeV); Photons Out (%)',
                      len(lengths)-1,lengths[0],lengths[len(lengths)-1],
                      len(pTs)-1,pTs[0],pTs[len(pTs)-1])
eventcounts = []

# directory loop
for dir in directories:
    xml = ET.parse(analysisDir+"points/"+dir+"/xml/brick_0.xml")
    root = xml.getroot()
    events = 0

    # reading length
    length = 0
    for temp in root.iter('brick_length'):
        length = float(temp.text)

    # reading partons
    datfile = open(analysisDir+"points/"+dir+"/dat/Brick_Bin.dat","r")
    datlines = datfile.readlines()
    for line in datlines:
        if 'Event' in line:
            events = events + 1

        if '#' in line:
            continue

        part = parton(line)
        if part.PID == 22 and part.E > 0.95:
            totalhist.Fill(length,part.pT)

    eventcounts.append(events)


# fits
totalhist.Scale(100.0/eventcounts[0])
for i in range(len(pTs)-1):
    name = str(pTs[i])+"to"+str(pTs[i+1])
    hist = totalhist.ProjectionX(name,i+1,i+1)

    hist.GetYaxis().SetTitle("Photons Out (%)") 
    hist.Fit("expo")
    hist.SetLineWidth(3)
    rp = ROOT.TRatioPlot(hist)
    fit = hist.GetFunction("expo")

    # output
    print("Mean free path: " + str(-1.0/fit.GetParameter(1)) + " +/- " + str(-1*fit.GetParError(1)/fit.GetParameter(1)) + " fm")
    c1 = ROOT.TCanvas("c1","c1",1400,1200)
    rp.Draw()
    c1.Print(outdir+name+".png")
    c1.Close()
    hist.Write()

# end behavior
totalhist.Write()
outfile.Close()