"""
  utils for histogramming / plotting
  
  Author: James Mulligan (james.mulligan@berkeley.edu)
"""

# General
import os
import sys
import yaml
import argparse

# Data analysis and plotting
import ROOT
import ctypes
import numpy as np

# Base class
sys.path.append('.')
from jetscape_analysis.base import common_base

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

################################################################
class PlotUtils(common_base.CommonBase):

    # ---------------------------------------------------------------
    # Constructor
    # ---------------------------------------------------------------
    def __init__(self, **kwargs):
        super(PlotUtils, self).__init__(**kwargs)

    # ---------------------------------------------------------------
    # Get bin array specified in config block
    # ---------------------------------------------------------------
    def bins_from_config(self, block, sqrts, observable_type, observable, suffix=''):
    
        if 'bins' in block:
            print(f'  Histogram with custom binning found for {observable} {suffix}')
            return np.array(block['bins'])
        elif f'hepdata' in block:
            print(f'  Histogram with hepdata binning found for {observable} {suffix}')
            return self.bins_from_hepdata(block, sqrts, observable_type, observable, suffix)
        else:
            print(f'  Warning: No binning found for {observable} {suffix}')
            return np.array([])
            
    # ---------------------------------------------------------------
    # Get bin array from hepdata file specified in config block
    # ---------------------------------------------------------------
    def bins_from_hepdata(self, block, sqrts, observable_type, observable, suffix=''):

        # Open the HEPData file
        hepdata_dir = f'plot/data/{sqrts}/{observable_type}/{observable}'
        hepdata_filename = os.path.join(hepdata_dir, block['hepdata'])
        f = ROOT.TFile(hepdata_filename, 'READ')
        
        # Find the relevant directory:
        # - The list of dir/hist names may also contain a suffix,
        #   which specifies e.g. the pt bin, jetR, or other parameters
        
        # First, check for dir/hist names in config
        if f'hepdata_AA_dir{suffix}' in block:
            dir_key = f'hepdata_AA_dir{suffix}'
        elif f'hepdata_AA_dir' in block:
            dir_key = f'hepdata_AA_dir'
        else:
            print(f'hepdata_AA_dir{suffix} not found!')
            
        if f'hepdata_AA_hname{suffix}' in block:
            h_key = f'hepdata_AA_hname{suffix}'
        elif f'hepdata_AA_hname' in block:
            h_key = f'hepdata_AA_hname'
        else:
            print(f'hepdata_AA_hname{suffix} not found!')
            
        # Get the appropriate entry in the dir/hist list
        if type(block[dir_key]) is list:
            dir_name = block[dir_key][0]
            h_name = block[h_key]
        else:
            dir_name = block[dir_key]
            h_name = block[h_key][0]
        
        # Get the histogram, and return the bins
        dir = f.Get(dir_name)
        h = dir.Get(h_name)
        bins = np.array(h.GetXaxis().GetXbins())

        # For certain Soft Drop observables, we need to exclude the "untagged" bin so that it will become underflow
        if observable == 'zg_alice' or observable == 'tg_alice':
            bins = bins[1:]

        f.Close()
        
        return bins
        
    # ---------------------------------------------------------------
    # Get tgraph from hepdata file specified in config block
    # ---------------------------------------------------------------
    def tgraph_from_hepdata(self, block, sqrts, observable_type, observable, suffix='', pt_suffix=''):

        # Open the HEPData file
        hepdata_dir = f'plot/data/{sqrts}/{observable_type}/{observable}'
        hepdata_filename = os.path.join(hepdata_dir, block['hepdata'])
        f = ROOT.TFile(hepdata_filename, 'READ')
        
        # Find the relevant directory:
        # - The list of dir/hist names may contain a suffix,
        #   which specifies e.g. the pt bin, jetR, or other parameters

        # First, check for dir names in config
        if f'hepdata_pp_dir{suffix}' in block:
            dir_key = f'hepdata_pp_dir{suffix}'
        elif f'hepdata_pp_dir{suffix}{pt_suffix}' in block:
            dir_key = f'hepdata_pp_dir{suffix}{pt_suffix}'
        elif f'hepdata_pp_dir' in block:
            dir_key = f'hepdata_pp_dir'
        else:
            #print(f'hepdata_pp_dir{suffix} not found!')
            return None
            
        # Check for hist names in config
        if f'hepdata_pp_gname{suffix}' in block:
            g_key = f'hepdata_pp_gname{suffix}'
        elif f'hepdata_pp_gname{suffix}{pt_suffix}' in block:
            g_key = f'hepdata_pp_gname{suffix}{pt_suffix}'
        elif f'hepdata_pp_gname' in block:
            g_key = f'hepdata_pp_gname'
        else:
            #print(f'hepdata_pp_gname{suffix} not found!')
            return None
            
        # Get the appropriate centrality entry in the dir/hist list
        if type(block[dir_key]) is list:
            dir_name = block[dir_key][0]
            g_name = block[g_key]    
        elif type(block[g_key]) is list:
            dir_name = block[dir_key]
            g_name = block[g_key][0]
        else:
            dir_name = block[dir_key]
            g_name = block[g_key]
        
        # Get the tgraph, and return the bins
        dir = f.Get(dir_name)
        g = dir.Get(g_name)
        f.Close()
        
        return g
        
    #---------------------------------------------------------------
    # Divide a histogram by a tgraph, point-by-point
    #---------------------------------------------------------------
    def divide_histogram_by_tgraph(self, h, g):

        # Clone tgraph, in order to return a new one
        g_new = g.Clone(f'{h.GetName()}_divided')
    
        nBins = h.GetNbinsX()
        for bin in range(1, nBins+1):

            # Get histogram (x,y)
            h_x = h.GetBinCenter(bin)
            h_y = h.GetBinContent(bin)
            h_error = h.GetBinError(bin)

            # Get TGraph (x,y) and errors
            g_x = ctypes.c_double(0)
            g_y = ctypes.c_double(0)
            g.GetPoint(bin-1, g_x, g_y)
            yErrLow = g.GetErrorYlow(bin-1)
            yErrUp  = g.GetErrorYhigh(bin-1)
            
            gx = g_x.value
            gy = g_y.value

            # Skip if tgraph starts below hist (since hist has min cut)
            if gx < h_x:
                continue
                
            if not np.isclose(h_x, gx):
                print(f'ERROR: hist x: {h_x}, graph x: {gx} -- will not plot ratio')
                return None
          
            new_content = h_y / gy
            
            # Combine tgraph and histogram relative uncertainties in quadrature
            if gy > 0. and h_y > 0.:
                new_error_low = np.sqrt( pow(yErrLow/gy,2) + pow(h_error/h_y,2) ) * new_content
                new_error_up = np.sqrt( pow(yErrUp/gy,2) + pow(h_error/h_y,2) ) * new_content
            else:
                new_error_low = (yErrLow/gy) * new_content
                new_error_up = (yErrUp/gy) * new_content

            g_new.SetPoint(bin-1, h_x, new_content)
            g_new.SetPointError(bin-1, 0, 0, new_error_low, new_error_up)
        return g_new

    #---------------------------------------------------------------
    # Divide a tgraph by a tgraph, point-by-point: g1/g2
    # NOTE: Ignore uncertainties on denominator
    #---------------------------------------------------------------
    def divide_tgraph_by_tgraph(self, g1, g2):
    
        # Clone tgraph, in order to return a new one
        g_new = g1.Clone(f'{g1.GetName()}_divided')
        
        if g1.GetN() != g2.GetN():
            sys.exit(f'ERROR: TGraph {g1.GetName()} has {g1.GetN()} points, but {g2.GetName()} has {g2.GetN()} points')
        
        for i in range(0, g1.GetN()):

            # Get TGraph (x,y) and errors
            g1_x = ctypes.c_double(0)
            g1_y = ctypes.c_double(0)
            g1.GetPoint(i, g1_x, g1_y)
            y1ErrLow = g1.GetErrorYlow(i)
            y1ErrUp  = g1.GetErrorYhigh(i)
            g1x = g1_x.value
            g1y = g1_y.value
            
            g2_x = ctypes.c_double(0)
            g2_y = ctypes.c_double(0)
            g2.GetPoint(i, g2_x, g2_y)
            g2x = g2_x.value
            g2y = g2_y.value
            
            if not np.isclose(g1x, g2x):
                sys.exit(f'ERROR: TGraph {g1.GetName()} point {i} at {g1x}, but {g2.GetName()} at {g2x}')

            new_content = g1y / g2y
            new_error_low = y1ErrLow/g1y * new_content
            new_error_up = y1ErrUp/g1y * new_content
            
            g_new.SetPoint(i, g1x, new_content)
            g_new.SetPointError(i, 0, 0, new_error_low, new_error_up)
        return g_new

    #-------------------------------------------------------------------------------------------
    # Set legend parameters
    #-------------------------------------------------------------------------------------------
    def setup_legend(self, leg, textSize, sep):

        leg.SetTextFont(42)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetFillColor(0)
        leg.SetMargin(0.25)
        leg.SetTextSize(textSize)
        leg.SetEntrySeparation(sep)

    #-------------------------------------------------------------------------------------------
    def setOptions(self):

        font = 42

        ROOT.gStyle.SetFrameBorderMode(0)
        ROOT.gStyle.SetFrameFillColor(0)
        ROOT.gStyle.SetCanvasBorderMode(0)
        ROOT.gStyle.SetPadBorderMode(0)
        ROOT.gStyle.SetPadColor(10)
        ROOT.gStyle.SetCanvasColor(10)
        ROOT.gStyle.SetTitleFillColor(10)
        ROOT.gStyle.SetTitleBorderSize(1)
        ROOT.gStyle.SetStatColor(10)
        ROOT.gStyle.SetStatBorderSize(1)
        ROOT.gStyle.SetLegendBorderSize(1)

        ROOT.gStyle.SetDrawBorder(0)
        ROOT.gStyle.SetTextFont(font)
        ROOT.gStyle.SetStatFont(font)
        ROOT.gStyle.SetStatFontSize(0.05)
        ROOT.gStyle.SetStatX(0.97)
        ROOT.gStyle.SetStatY(0.98)
        ROOT.gStyle.SetStatH(0.03)
        ROOT.gStyle.SetStatW(0.3)
        ROOT.gStyle.SetTickLength(0.02,"y")
        ROOT.gStyle.SetEndErrorSize(3)
        ROOT.gStyle.SetLabelSize(0.05,"xyz")
        ROOT.gStyle.SetLabelFont(font,"xyz")
        ROOT.gStyle.SetLabelOffset(0.01,"xyz")
        ROOT.gStyle.SetTitleFont(font,"xyz")
        ROOT.gStyle.SetTitleOffset(1.2,"xyz")
        ROOT.gStyle.SetTitleSize(0.045,"xyz")
        ROOT.gStyle.SetMarkerSize(1)
        ROOT.gStyle.SetPalette(1)

        ROOT.gStyle.SetOptTitle(0)
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptFit(0)
