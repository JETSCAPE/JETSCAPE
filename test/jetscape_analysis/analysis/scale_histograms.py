# Macro to scale histograms of all Pt-hard bins, using xsec from JetscapeAnalysis task.
# This script expects files X/AnalysisResults.root, and will output scaled histograms
# to the same file, in a new output list with suffix "Scaled". The script will automatically loop over
# all output lists, subject to some simple criteria that covers basic use cases (can be adapted as needed).
#
# There is an option "bRemoveOutliers" to remove outliers from certain histograms. The features are
# currently hard-coded below so you will need to modify the code as needed. This feature is adapted from code of Raymond Ehlers.
#
# Author: James Mulligan (james.mulligan@yale.edu)
#

import ctypes
import os

import ROOT

# Prevent ROOT from stealing focus when plotting
ROOT.gROOT.SetBatch(True)

###################################################################################
# Main function
def scale_histograms(outputDirBin, bin, bRemoveOutliers=False):

    # Option to remove outliers from specified histograms
    # If the average bin content stays below the "outlierLimit" for "outlierNBinsThreshold" bins, it is removed
    outlierLimit = 2
    outlierNBinsThreshold = 4

    # Option to print out detailed info about scaling and outlier removal
    verbose = False

    # Read the cross-section, and scale histograms
    print(f"ooo Scaling Pt-hard bin {bin+1}")
    filename = os.path.join(outputDirBin, 'AnalysisResults.root')
    f = ROOT.TFile(filename, "UPDATE")
    n_events = f.Get("hNevents").GetBinContent(bin + 1)

    # Note that we should be careful of getting cross-section from file
    # since if files are hadd'ed they will contain the wrong value
    cross_section = f.Get("hCrossSection").GetBinContent(bin + 1)

    scaleFactor = cross_section / n_events
    print(f"ooo scaleFactor: {scaleFactor}  (cross_section={cross_section}, n_events={n_events})")

    # Now, scale all the histograms
    keys = f.GetListOfKeys()
    for key in keys:
        name = key.GetName()
        if "Scaled" in name:
            continue
        obj = f.Get(name)
        if obj:
            scale_all_histograms(
                obj,
                scaleFactor,
                f,
                verbose,
                outputDirBin,
                bRemoveOutliers,
                outlierLimit,
                outlierNBinsThreshold,
                bin,
                name
            )
        else:
            print("obj not found!")

        obj.Write("%sScaled" % obj.GetName())

    f.Close()


###################################################################################
# Function to iterate recursively through an object to scale all TH1/TH2/THnSparse
def scale_all_histograms(
    obj,
    scaleFactor,
    f,
    verbose,
    outputDirBin,
    bRemoveOutliers=False,
    limit=2,
    nBinsThreshold=4,
    pTHardBin=0,
    taskName=""
):

    # Set Sumw2 if not already done
    if obj.GetSumw2N() == 0:
        obj.Sumw2()
        if verbose:
            print("Set Sumw2 on {}".format(obj.GetName()))

    if obj.InheritsFrom(ROOT.TProfile.Class()):
        if verbose:
            print("TProfile %s not scaled..." % obj.GetName())
    elif obj.InheritsFrom(ROOT.TH2.Class()):
        obj.Scale(scaleFactor)
        if verbose:
            print("TH2 %s was scaled..." % obj.GetName())
    elif obj.InheritsFrom(ROOT.TH1.Class()):
        if bRemoveOutliers:
            name = obj.GetName()
            # only perform outlier removal on these couple histograms
            if "Pt" in name:
                remove_outliers(
                    pTHardBin, obj, verbose, outputDirBin, limit, nBinsThreshold, 1, taskName,
                )
        obj.Scale(scaleFactor)
        if verbose:
            print("TH1 %s was scaled..." % obj.GetName())
    elif obj.InheritsFrom(ROOT.THnSparse.Class()):
        obj.Scale(scaleFactor)
        if verbose:
            print("THnSparse %s was scaled..." % obj.GetName())
    else:
        if verbose:
            print("Not a histogram!")
            print(obj.GetName())
        for subobj in obj:
            ScaleAllHistograms(
                subobj,
                scaleFactor,
                f,
                verbose,
                outputDirBin,
                bRemoveOutliers,
                limit,
                nBinsThreshold,
                pTHardBin,
                taskName,
            )


###################################################################################
# Function to remove outliers from a TH3 (i.e. truncate the spectrum), based on projecting to the y-axis
# It truncates the 3D histogram based on when the 1D projection 4-bin moving average has been above
# "limit" for "nBinsThreshold" bins.
def remove_outliers(
    pTHardBin, hist, verbose, outputDirBin, limit=2, nBinsThreshold=4, dimension=3, taskName="",
):

    # Project to the pT Truth axis
    if dimension == 3:
        histToCheck = hist.ProjectionY("{}_projBefore".format(hist.GetName()))
    if dimension == 2:
        histToCheck = hist.ProjectionX("{}_projBefore".format(hist.GetName()))
    if dimension == 1:
        histToCheck = hist

    # Check with moving average
    foundAboveLimit = False
    cutLimitReached = False
    # The cut index is where we decided cut on that row
    cutIndex = -1
    nBinsBelowLimitAfterLimit = 0
    # nBinsThreshold= n bins that are below threshold before all bins are cut

    if verbose:
        (preMean, preMedian) = mean_and_median(histToCheck)

    for index in range(0, histToCheck.GetNcells()):
        if verbose:
            print("---------")
        avg = moving_average(histToCheck, index=index, numberOfCountsBelowIndex=2, numberOfCountsAboveIndex=2,)
        if verbose:
            print(
                "Index: {0}, Avg: {1}, BinContent: {5}, foundAboveLimit: {2}, cutIndex: {3}, cutLimitReached: {4}".format(
                    index, avg, foundAboveLimit, cutIndex, cutLimitReached, histToCheck.GetBinContent(index),
                )
            )
        if avg > limit:
            foundAboveLimit = True

        if not cutLimitReached:
            if foundAboveLimit and avg <= limit:
                if cutIndex == -1:
                    cutIndex = index
                nBinsBelowLimitAfterLimit += 1

            if nBinsBelowLimitAfterLimit != 0 and avg > limit:
                # Reset
                cutIndex = -1
                nBinsBelowLimitAfterLimit = 0

            if nBinsBelowLimitAfterLimit > nBinsThreshold:
                cutLimitReached = True
                break  # no need to continue the loop - we found our cut index

    # Do not perform removal here because then we miss values between the avg going below
    # the limit and crossing the nBinsThreshold
    if verbose:
        print("Hist checked: {0}, cut index: {1}".format(histToCheck.GetName(), cutIndex))

    # Use on both TH1 and TH2 since we don't start removing immediately, but instead only after the limit
    if cutLimitReached:
        if verbose:
            print("--> --> --> Removing outliers")
        # Check for values above which they should be removed by translating the global index
        x = ctypes.c_int(0)
        y = ctypes.c_int(0)
        z = ctypes.c_int(0)
        for index in range(0, hist.GetNcells()):
            # Get the bin x, y, z from the global bin
            hist.GetBinXYZ(index, x, y, z)
            if dimension == 3:
                if y.value >= cutIndex:
                    if hist.GetBinContent(index) > 1e-3:
                        if verbose:
                            print("Cutting for index {}. y bin {}. Cut index: {}".format(index, y, cutIndex))
                        hist.SetBinContent(index, 0)
                        hist.SetBinError(index, 0)
            if dimension == 2:
                # for the response matrix the pT Truth is on the y-Axis
                if hist.GetName() == "hResponseMatrixEMCal":
                    x.value = y.value
                if x.value >= cutIndex:
                    if hist.GetBinContent(index) > 1e-3:
                        if verbose:
                            print("Cutting for index {}. x bin {}. Cut index: {}".format(index, x, cutIndex))
                        hist.SetBinContent(index, 0)
                        hist.SetBinError(index, 0)
            if dimension == 1:
                if x.value >= cutIndex:
                    if hist.GetBinContent(index) > 1e-3:
                        if verbose:
                            print("Cutting for index {}. x bin {}. Cut index: {}".format(index, x, cutIndex))
                        hist.SetBinContent(index, 0)
                        hist.SetBinError(index, 0)

    else:
        if verbose:
            print("Hist {} did not have any outliers to cut".format(hist.GetName()))

    # Check the mean and median
    # Use another temporary hist
    if dimension == 3:
        histToCheckAfter = hist.ProjectionY()
    if dimension == 2:
        if hist.GetName() == "hResponseMatrixEMCal":
            histToCheckAfter = hist.ProjectionY()
        else:
            histToCheckAfter = hist.ProjectionX()
    if dimension == 1:
        histToCheckAfter = hist

    if verbose:
        (postMean, postMedian) = get_hist_mean_and_median(histToCheckAfter)
        print("Pre  outliers removal mean: {}, median: {}".format(preMean, preMedian))
        print("Post outliers removal mean: {}, median: {}".format(postMean, postMedian))
    outlierFilename = "{}OutlierRemoval_{}.pdf".format(outputDirBin, hist.GetName())
    if "Pt" in hist.GetName():
        plot_outlier_PDF(
            histToCheck, histToCheckAfter, pTHardBin, outlierFilename, verbose, "hist E", True,
        )


########################################################################################################
def get_hist_mean_and_median(hist):
    # Median
    # See: https://root-forum.cern.ch/t/median-of-histogram/7626/5
    x = ctypes.c_double(0)
    q = ctypes.c_double(0.5)
    # Apparently needed to be safe(?)
    hist.ComputeIntegral()
    hist.GetQuantiles(1, x, q)

    mean = hist.GetMean()
    return (mean, x.value)


########################################################################################################
def moving_average(hist, index, numberOfCountsBelowIndex=0, numberOfCountsAboveIndex=2):
    """
  # [-2, 2] includes -2, -1, 0, 1, 2
  """
    # Check inputs
    if numberOfCountsBelowIndex < 0 or numberOfCountsAboveIndex < 0:
        print("Moving average number of counts above or below must be >= 0. Please check the values!")

    count = 0.0
    average = 0.0
    for i in range(index - numberOfCountsBelowIndex, index + numberOfCountsAboveIndex + 1):
        # Avoid going over histogram limits
        if i < 0 or i >= hist.GetNcells():
            continue
        # print("Adding {}".format(hist.GetBinContent(i)))
        average += hist.GetBinContent(i)
        count += 1

    # if count != (numberOfCountsBelowIndex + numberOfCountsAboveIndex + 1):
    #    print("Count: {}, summed: {}".format(count, (numberOfCountsBelowIndex + numberOfCountsAboveIndex + 1)))
    # exit(0)

    return average / count


########################################################################################################
# Plot basic histogram    ##############################################################################
########################################################################################################
def plot_outlier_PDF(
    h, hAfter, pTHardBin, outputFilename, verbose, drawOptions="", setLogy=False,
):

    c = ROOT.TCanvas("c", "c: hist", 600, 450)
    c.cd()
    if setLogy:
        c.SetLogy()
    h.GetXaxis().SetRangeUser(0, 250)
    h.Draw("hist")
    h.SetLineColor(616)
    hAfter.SetLineColor(820)
    hAfter.Draw("same hist")

    leg1 = ROOT.TLegend(0.17, 0.7, 0.83, 0.85, "outlier removal of Bin {}".format(pTHardBin))
    leg1.SetFillColor(10)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetTextSize(0.04)
    leg1.AddEntry(h, "before", "l")
    leg1.AddEntry(hAfter, "after", "l")
    leg1.Draw("same")

    c.Print("{}".format(outputFilename))
    c.Close()

# ---------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    print("Executing scale_histograms.py...")
    print("")
