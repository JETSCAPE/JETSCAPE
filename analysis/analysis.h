#ifndef ANALYSIS_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define ANALYSIS_H

//C++ header
#include "string"
#include <iostream>
#include <fstream>
#include "iomanip"
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>
#include <sstream>

#include <memory>
#include <chrono>
#include <thread>

#include "gzstream.h"
#include "PartonShower.h"
#include "JetScapeLogger.h"
#include "JetScapeReader.h"
#include "fjcore.hh"
#include "Pythia8/Pythia.h"

#include <GTL/dfs.h>
//ROOT headers
#include <TH1.h>
#include <TFile.h>
#include <TVector.h>
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TRatioPlot.h"

using namespace Jetscape;

vector<string> splitString(string input, string divider);

TGraphErrors histToGraph(TH1D* hist);

void ratioPlot(TH1D* dataHist, TH1D* predictionHist, string title, bool xlog = false);

void ratioPlot(TGraphErrors* dataHist, TH1D* predictionHist, string title, bool xlog = false, bool ylog = false, string xname = "pT (GeV)");

vector<double> getThrustSphericity(vector<shared_ptr<Hadron>> hadrons);

int histToCSV(TH1D* hist, string title);

void makeDatFile(vector<vector<double>> data, string title, string header);

void scaleBins(TH1D* hist, double scale);

std::vector<std::string> get_directories(const std::string& s);

std::vector<std::vector<string>> get_dat_bounds(const std::string& s);

std::vector<std::string> doubleSort(std::vector<std::string> input);

void smoothBins(TH1D* hist);

void smoothBinsOld(TH1D* hist);

//reading first line of param file to get center of mass energy
string getEcm(string baseDir);

//Getting list of xsecs for a run based of the ECM in the param file and pthat bounds fed in
vector<vector<double>> getXsecs(vector<string> pTHatMin, vector<string> pTHatMax);

//integral for root histograms
double integral(TH1D* hist);

//assigning weights to hadrons to blend hard and soft QCD
double blending(int binID, double pT, double center, double width, double softcut);

#endif