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

# initializing histograms
hist = ROOT.TH1D('Photon Rate','Photon Rate; Brick Length (fm); Photons Out',4,1.0,9.0)
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
        if part.PID == 22:
            hist.Fill(length)

    eventcounts.append(events)

# fits
hist.Scale(100.0/eventcounts[0])
hist.Fit("expo")
hist.SetLineWidth(3)
rp = ROOT.TRatioPlot(hist)
fit = hist.GetFunction("expo")

# output
print("Mean free path: " + str(-1.0/fit.GetParameter(1)) + " +/- " + str(-1*fit.GetParError(1)/fit.GetParameter(1)) + " fm")
c1 = ROOT.TCanvas("c1","c1",1400,1200)
rp.Draw()
c1.Print(analysisDir+"fit.png")
c1.Close()

# end behavior
outfile = ROOT.TFile.Open(analysisDir+'analysis.root', "RECREATE")
hist.Write()
outfile.Close()