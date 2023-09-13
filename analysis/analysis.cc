#include "analysis.h"
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
#include "TLine.h"
#include "TProfile.h"

using namespace Jetscape;

//method to split line inputs by a divider into a vector of individual strings
vector<string> splitString(string input, string divider){
	vector<string> output;
	int pos = 0;
	
	while( pos < input.length() ){
		int end = input.find(divider, pos);
		output.push_back(input.substr(pos,end-pos));
		pos = end + 1;
        if(end == -1) break;
    }
	
	return output;
}

//converts histogram to a tgrapherrors for plotting
//TGraphErrors* exampleGraph = (TGraphErrors*)histToGraph(exampleHist).Clone();
TGraphErrors histToGraph(TH1D* hist){
    int n = hist->GetNbinsX();
    double xcoords[n], ycoords[n], xerrors[n], yerrors[n];

    for(int i = 1; i <= n; i++){
        xcoords[i-1] = hist->GetBinCenter(i);
        ycoords[i-1] = hist->GetBinContent(i);
        xerrors[i-1] = hist->GetBinWidth(i)/2;
        yerrors[i-1] = hist->GetBinError(i);
    }

    TGraphErrors* graph = new TGraphErrors(n,xcoords,ycoords,xerrors,yerrors);
    return *graph;
}

//custom ratio plot generator since existing one is broken
void ratioPlot(TH1D* dataHist, TH1D* predictionHist, string title, bool xlog){
    //values for plots
    int bins = dataHist->GetNbinsX();
    double data[bins], prediction[bins], xcoords[bins], dataErrors[bins], binWidths[bins], asym[bins], asymErrors[bins];

    //cycling through bins, note first one is an overflow bin so it is skipped
    for(int i = 0; i < bins; i++){
        //reading values from histograms
        data[i] = dataHist->GetBinContent(i+1);
        dataErrors[i] = dataHist->GetBinError(i+1);
        prediction[i] = predictionHist->GetBinContent(i+1);
        xcoords[i] = predictionHist->GetBinCenter(i+1);
        binWidths[i] = predictionHist->GetBinWidth(i+1) / 2;
        asym[i] = (prediction[i]/data[i]);
        asymErrors[i] = dataErrors[i]/prediction[i];

        //debugging
        //cout << xcoords[i] << " " << dataHist->GetBinContent(i+1)*1000000 << " " << dataHist->GetBinError(i+1)*1000000 << " " << asym[i] << " " << binWidths[i] << endl;
        
        //correcting for negatives when doing log
        //if(data[i] < 0) data[i] = 0;
        //if(dataErrors[i] > data[i]) dataErrors[i] = data[i]*0.99;
    }

    //Data graph
    TGraphErrors* dataPlot = new TGraphErrors(bins,xcoords,data,binWidths,dataErrors);
    dataPlot->SetMarkerStyle(kFullDotLarge);
    dataPlot->SetMarkerColor(kRed);
    dataPlot->SetLineColor(kRed);

    //prediction graph
    TGraph* predictionPlot = new TGraph(bins,xcoords,prediction);
    predictionPlot->SetLineColor(kBlue);
    predictionPlot->SetLineWidth(3);

    //layered graph
    TMultiGraph* total = new TMultiGraph();
    total->Add(predictionPlot,"l");
    total->Add(dataPlot,"AP");
    total->SetTitle("Jetscape Comparison");
    total->GetYaxis()->SetTitle(title.c_str());
    total->GetXaxis()->SetRangeUser(xcoords[0]-binWidths[0],xcoords[bins-1]+binWidths[bins-1]);
    //total->GetYaxis()->SetTitleSize(0.05);
    total->GetYaxis()->SetLabelSize(0.04);

    //ratio plot
    TGraphErrors* comparisonPlot = new TGraphErrors(bins, xcoords, asym, binWidths, asymErrors);
    comparisonPlot->SetTitle("");
    comparisonPlot->GetXaxis()->SetRangeUser(xcoords[0]-binWidths[0],xcoords[bins-1]+binWidths[bins-1]);
    comparisonPlot->GetXaxis()->SetTitle("pT (GeV)");
    comparisonPlot->GetXaxis()->SetTitleSize(0.06);
    comparisonPlot->GetXaxis()->SetLabelSize(0.06);
    comparisonPlot->GetYaxis()->SetTitle("Asymmetry");
    comparisonPlot->GetYaxis()->SetTitleSize(0.06);
    comparisonPlot->GetYaxis()->SetLabelSize(0.06);
    comparisonPlot->SetMarkerStyle(kFullDotLarge);

    //Drawing main plot
    TCanvas* c = new TCanvas("c1","c1",1000,800);
    TPad* upper = new TPad("plot","plot",0,0.4,1,1);
    upper->SetBottomMargin(0);
    upper->Draw();
    upper->cd();
    upper->SetLogy();
    if(xlog) upper->SetLogx();
	total->Draw("apl");

    //Legend
    TLegend leg(.7,.7,.9,.9,"Sources");
    leg.AddEntry(predictionPlot,"Jetscape","l");
    leg.AddEntry(dataPlot,"CMS Data","ep");
    leg.DrawClone("Same");

    //drawing ratio plot
    c->cd();
    TPad* lower = new TPad("plot","plot",0,0,1,0.4);
    lower->SetTopMargin(0);
    lower->SetBottomMargin(1);
    lower->Draw();
    lower->cd();
    if(xlog) lower->SetLogx();
	comparisonPlot->Draw("AP");
    TLine* line = new TLine(xcoords[0],1,xcoords[bins-1]+binWidths[bins-1],1);
    line->SetLineStyle(9);
    line->Draw();

    string filename  = splitString(title," (")[0]+".png"; //trims off the units in the file name
    c->Print(filename.c_str());
    c->Close();
}

//overload for graph input
void ratioPlot(TGraphErrors* dataHist, TH1D* predictionHist, string title, bool xlog, bool ylog, string xname){
    //values for plots
    int bins = predictionHist->GetNbinsX();
    double data[bins], prediction[bins], xcoords[bins], dataErrors[bins], binWidths[bins], asym[bins], asymErrors[bins];

    //cycling through bins, note first one is an overflow bin so it is skipped
    for(int i = 0; i < bins; i++){
        //reading values from histograms
        prediction[i] = predictionHist->GetBinContent(i+1);
        xcoords[i] = predictionHist->GetBinCenter(i+1);
        binWidths[i] = predictionHist->GetBinWidth(i+1) / 2;
        data[i] = dataHist->Eval(xcoords[i]);
        dataErrors[i] = dataHist->GetErrorY(i+1);
        asym[i] = (prediction[i]/data[i]);
        asymErrors[i] = dataErrors[i]/prediction[i];

        //debugging
        //cout << xcoords[i] << " " << dataHist->GetBinContent(i+1)*1000000 << " " << dataHist->GetBinError(i+1)*1000000 << " " << asym[i] << " " << binWidths[i] << endl;
        
        //correcting for negatives when doing log
        //if(data[i] < 0) data[i] = 0;
        //if(dataErrors[i] > data[i]) dataErrors[i] = data[i]*0.99;
    }

    //Data graph
    TGraphErrors* dataPlot = (TGraphErrors*) dataHist->Clone();
    dataPlot->SetMarkerStyle(kFullDotLarge);
    dataPlot->SetMarkerColor(kRed);
    dataPlot->SetLineColor(kRed);

    //prediction graph
    TGraph* predictionPlot = new TGraph(bins,xcoords,prediction);
    predictionPlot->SetLineColor(kBlue);
    predictionPlot->SetLineWidth(3);

    //layered graph
    TMultiGraph* total = new TMultiGraph();
    total->Add(predictionPlot,"l");
    total->Add(dataPlot,"AP");
    total->SetTitle("Jetscape Comparison");
    total->GetYaxis()->SetTitle(title.c_str());
    total->GetXaxis()->SetRangeUser(xcoords[0]-binWidths[0],xcoords[bins-1]+binWidths[bins-1]);
    //total->GetYaxis()->SetTitleSize(0.05);
    total->GetYaxis()->SetLabelSize(0.04);

    //ratio plot
    TGraphErrors* comparisonPlot = new TGraphErrors(bins, xcoords, asym, binWidths, asymErrors);
    comparisonPlot->SetTitle("");
    comparisonPlot->GetXaxis()->SetRangeUser(xcoords[0]-binWidths[0],xcoords[bins-1]+binWidths[bins-1]);
    comparisonPlot->GetYaxis()->SetRangeUser(0,2);
    comparisonPlot->GetXaxis()->SetTitle(xname.c_str());
    comparisonPlot->GetXaxis()->SetTitleSize(0.06);
    comparisonPlot->GetXaxis()->SetLabelSize(0.06);
    comparisonPlot->GetYaxis()->SetTitle("Asymmetry");
    comparisonPlot->GetYaxis()->SetTitleSize(0.06);
    comparisonPlot->GetYaxis()->SetLabelSize(0.06);
    comparisonPlot->SetMarkerStyle(kFullDotLarge);

    //Drawing main plot
    TCanvas* c = new TCanvas("c1","c1",1000,800);
    TPad* upper = new TPad("plot","plot",0,0.4,1,1);
    upper->SetBottomMargin(0);
    upper->Draw();
    upper->cd();
    if(ylog) upper->SetLogy();
    if(xlog) upper->SetLogx();
	total->Draw("apl");

    //Legend
    TLegend leg(.7,.7,.9,.9,"Sources");
    leg.AddEntry(predictionPlot,"Jetscape","l");
    leg.AddEntry(dataPlot,"Data","ep");
    leg.DrawClone("Same");

    //drawing ratio plot
    c->cd();
    TPad* lower = new TPad("plot","plot",0,0,1,0.4);
    lower->SetTopMargin(0);
    lower->SetBottomMargin(1);
    lower->Draw();
    lower->cd();
    if(xlog) lower->SetLogx();
	comparisonPlot->Draw("AP");
    TLine* line = new TLine(xcoords[0],1,xcoords[bins-1]+binWidths[bins-1],1);
    line->SetLineStyle(9);
    line->Draw();

    string filename  = "plots/"+splitString(title," (")[0]+".png"; //trims off the units in the file name
    c->Print(filename.c_str());
    c->Close();
}

//thrust calculation classes and methods
class Vec4 {
    public:

    // Constructors.
    Vec4(double xIn = 0., double yIn = 0., double zIn = 0., double tIn = 0.) : xx(xIn), yy(yIn), zz(zIn), tt(tIn) { }
    Vec4(const Vec4& v) : xx(v.xx), yy(v.yy), zz(v.zz), tt(v.tt) { }
    Vec4& operator=(const Vec4& v) { if (this != &v) { xx = v.xx; yy = v.yy; zz = v.zz; tt = v.tt; } return *this; }
    Vec4& operator=(double value) { xx = value; yy = value; zz = value; tt = value; return *this; }

    // Member functions for input.
    void reset() {xx = 0.; yy = 0.; zz = 0.; tt = 0.;}
    void p(double xIn, double yIn, double zIn, double tIn) {xx = xIn; yy = yIn; zz = zIn; tt = tIn;}
    void p(Vec4 pIn) {xx = pIn.xx; yy = pIn.yy; zz = pIn.zz; tt = pIn.tt;}
    void px(double xIn) {xx = xIn;}
    void py(double yIn) {yy = yIn;}
    void pz(double zIn) {zz = zIn;}
    void e(double tIn) {tt = tIn;}

    // Member functions for output.
    double px() const {return xx;}
    double py() const {return yy;}
    double pz() const {return zz;}
    double e() const {return tt;}
    double& operator[](int i) {
        if      (i == 1) return xx;
        else if (i == 2) return yy;
        else if (i == 3) return zz;
        else             return tt;
    }
    double pAbs() const {return sqrt(xx*xx + yy*yy + zz*zz);}
    double pAbs2() const {return xx*xx + yy*yy + zz*zz;}

    // Scalar and cross product of 3-vector parts.
    friend double dot3(const Vec4& v1, const Vec4& v2) {return v1.xx*v2.xx + v1.yy*v2.yy + v1.zz*v2.zz;}
    friend Vec4 cross3(const Vec4& v1, const Vec4& v2) {
        Vec4 v;
        v.xx = v1.yy * v2.zz - v1.zz * v2.yy;
        v.yy = v1.zz * v2.xx - v1.xx * v2.zz;
        v.zz = v1.xx * v2.yy - v1.yy * v2.xx;
        return v;
    }

    // Operator overloading with member functions
    inline Vec4 operator-() const {Vec4 tmp; tmp.xx = -xx; tmp.yy = -yy; tmp.zz = -zz; tmp.tt = -tt; return tmp;}
    inline Vec4& operator+=(const Vec4& v) {xx += v.xx; yy += v.yy; zz += v.zz; tt += v.tt; return *this;}
    inline Vec4& operator-=(const Vec4& v) {xx -= v.xx; yy -= v.yy; zz -= v.zz; tt -= v.tt; return *this;}
    inline Vec4& operator*=(double f) {xx *= f; yy *= f; zz *= f; tt *= f; return *this;}
    inline Vec4& operator/=(double f) {xx /= f; yy /= f; zz /= f; tt /= f; return *this;}
    inline Vec4 operator+(const Vec4& v) const {Vec4 tmp = *this; return tmp += v;}
    inline Vec4 operator-(const Vec4& v) const {Vec4 tmp = *this; return tmp -= v;}
    inline Vec4 operator*(double f) const {Vec4 tmp = *this; return tmp *= f;}
    inline Vec4 operator/(double f) const {Vec4 tmp = *this; return tmp /= f;}
    inline double operator*(const Vec4& v) const {return tt*v.tt - xx*v.xx - yy*v.yy - zz*v.zz;}


    private:

    // The four-vector data members.
    double xx, yy, zz, tt;
  
};

inline Vec4 operator*(double f, const Vec4& v1){Vec4 v = v1; return v *= f;}/*
// Scalar and cross product of 3-vector parts.
friend double dot3(const Vec4& v1, const Vec4& v2){return v1.xx*v2.xx + v1.yy*v2.yy + v1.zz*v2.zz;}
friend Vec4 cross3(const Vec4& v1, const Vec4& v2){
	Vec4 v;
	v.xx = v1.yy * v2.zz - v1.zz * v2.yy;
	v.yy = v1.zz * v2.xx - v1.xx * v2.zz;
	v.zz = v1.xx * v2.yy - v1.yy * v2.xx;
	return v;
}
*/
double pow2(double x){return x*x;}
double pow3(double x){return x*x*x;}

//determines if particle is charged or not
int getCharge(int p_id){
	int id[3] = {0,0,0}; int pos_id = std::abs(p_id);
	
	id[0] = (pos_id/1000)%10; id[1] = (pos_id/100)%10; id[2] = (pos_id/10)%10;
	
	if(pos_id > 5999){return false;}
	if(pos_id < 111){return ((pos_id == 11) || (pos_id == 13) || (pos_id == 15) || (pos_id == 17)) ? true : false;}
	
	//finding net charge:
	int charge = 0; //charge /3 is the actual charge, but we just want != 0
	bool meson=false; if(id[0]==0){meson=true;}
	for(int i=0;i<3;++i){
		int factor = 1; if(meson && i==2){factor = -1;}
		if(id[i] == 0){continue;}
		else if(id[i]%2 == 0){charge+=2*factor;}
		else{charge-=factor;}
	}
	
	return charge/3; //quark charges are scaled by 3
}

//returns vector of {thrust, sphericity} for a list of hadrons
vector<double> getThrustSphericity(vector<shared_ptr<Hadron>> hadrons){
    //error vector to return if cuts not met
    vector<double> error = {-1,-100};
    double phiDeg = -100;

    // Initializing vars used in thrust calculation
    // Initial values and counters zero.
    // Thr:
    double Thr_eVal1 = 0.;
    Vec4 Thr_eVec1 = 0.;
    int nStudy = 0;
    std::vector<Vec4> pOrder;
    Vec4 pSum, nRef, pPart, pFull, pMax;
    
    int nHad = 0; int nChg = 0;
    int runmore = 1; int curptn_event = 1; int prevptn_event = 1; int nevents = 0;
    double chgE = 0; double thr;
    double prevevt_thrust = 0.;

    // Sph: (default 'power' is 2 => default 'powerInt' is 2, 'powerMod' is 0)
    int powerInt = 2; double powerMod = 0.;

    double denom = 0.;
    double Sph_eVal1 = 0.; double Sph_eVal2 = 0.; double Sph_eVal3 = 0.;
    Vec4 Sph_eVec1 = 0.; Vec4 Sph_eVec2 = 0.; Vec4 Sph_eVec3 = 0.;

    double tt[4][4] = {0.0}; 

    for(unsigned int i=0; i<hadrons.size(); i++){
        //vars to store current particle in
        int SN = i;
        double pid = hadrons[i].get()->pid();
        double e  = hadrons[i].get()->e();
        double px = hadrons[i].get()->px();
        double py = hadrons[i].get()->py();
        double pz = hadrons[i].get()->pz();
        double Phi = hadrons[i].get()->phi();
        double pStat = hadrons[i].get()->pstat();
        double chg = getCharge(pid);
        double pT, p, eta;
        
        //if this particle isn't one that we want to include in the calculation, skip it
        //if(!ischarged(pid) || don'twantparticle || !isfinal){continue;}
        //if(chg == 0){continue;}
        
        // current particle is in *this* event - add it to current event listings and then move to next particle
        pT = sqrt(px*px + py*py);
        p = sqrt(px*px + py*py + pz*pz);
        eta = atanh(pz/p);

//CUTS (1)

        //Neutral Hadron Selection
        if (chg == 0 /*&& e > 0.4 && std::abs(eta) < 2.29*/) {
            nHad++;
        }
        else if (chg == 0) {continue;}

        //Charged Hadron Selection
        if (chg != 0 /*&& pT > 0.2 && std::abs(eta) < 1.74*/) {
            nChg++;
            nHad++;
            chgE += e;
        }
        else if (chg !=0) {continue;}

        // Thr event fill:
        Vec4 pNow(px,py,pz,e);
        pNow.e(pNow.pAbs());
        pSum += pNow;
        pOrder.push_back(pNow);
        prevevt_thrust = thr;
        
        //Sph event fill:
        double p2Now = px*px + py*py + pz*pz;
        //double pWeight = 1.;
        //if (powerInt == 1) pWeight = 1. / sqrt(max(P2MIN, p2Now));
        //else if (powerInt == 0) pWeight = pow( max(P2MIN, p2Now), powerMod);
        double pWeight = (powerInt == 1) ? (1./sqrt(std::max(1.0e-20,p2Now))) : 1.;
        pWeight = (powerInt == 0) ? pow(std::max(1.0e-20, p2Now),powerMod) : 1.;
        //for (int j=1;j<4;++j){for(int k=j;k<4;++k){tt[j][k] += pWeight * pNow[j] * pNow[k];} [1,1]px*px, [1,2]px*py, [1,3]px*pz, [2,2]py*py, [2,3]py*pz, [3,3]pz*pz
        tt[1][1]+=pWeight*px*px; 
        tt[1][2]+=pWeight*px*py; 
        tt[1][3]+=pWeight*px*pz; 
        tt[2][2]+=pWeight*py*py; 
        tt[2][3]+=pWeight*py*pz; 
        tt[3][3]+=pWeight*pz*pz;
        denom += pWeight * p2Now;
        
        ++nStudy; 
        
        if(i != hadrons.size() - 1) continue; //skips the rest unless this is the last particle
        
        ///////////////////////////////////////////////////////////////////////////////////////////////
        // only happens after youve reached the end of the list
        
        // if this event doesn't have enough particles to calculate thrust/sphericity, skip it (can also skip otherwise unwanted events here too)
        if(nStudy < 2 || nChg < 5 || nHad < 13 || chgE < 15.0){
            //std::cout << "Event " << prevptn_event << " had too few particles to find thrust...\n\n";
            //fileout << "0\n";
            return error;
        }
        
        //Thrust and sphericity calculations proper
        // Thr: Try all combinations of reference vector orthogonal to two particles.
        for (int i1 = 0; i1 < nStudy - 1; ++i1) {
            for (int i2 = i1 + 1; i2 < nStudy; ++i2) {
                nRef = cross3( pOrder[i1], pOrder[i2]);
                nRef /= nRef.pAbs();
                pPart = 0.;

                // Add all momenta with sign; two choices for each reference particle.
                for (int i = 0; i < nStudy; ++i) if (i != i1 && i != i2) {
                if (dot3(pOrder[i], nRef) > 0.) pPart += pOrder[i];
                else                            pPart -= pOrder[i];
                }
                for (int j = 0; j < 4; ++j) {
                if      (j == 0) pFull = pPart + pOrder[i1] + pOrder[i2];
                else if (j == 1) pFull = pPart + pOrder[i1] - pOrder[i2];
                else if (j == 2) pFull = pPart - pOrder[i1] + pOrder[i2];
                else             pFull = pPart - pOrder[i1] - pOrder[i2];
                pFull.e(pFull.pAbs());
                if (pFull.e() > pMax.e()) pMax = pFull;
                }
            }
        }
        
        //cout << pMax.e() << " " << pSum.e() << endl;
        // Thr: Maximum gives thrust axis and value.
        Thr_eVal1 = pMax.e() / pSum.e();
        Thr_eVec1 = pMax / pMax.e();
        Thr_eVec1.e(0.);

        
        // Sph: Normalize tensor to trace = 1.
        for (int j=1;j<4;++j) {
            for(int k=j;k<4;++k) {
                tt[j][k]/=denom;
            }
        }
        
        // Sph: Find eigenvalues to matrix (third degree equation).
        double qCoef = ( tt[1][1] * tt[2][2] + tt[1][1] * tt[3][3] + tt[2][2] * tt[3][3] - pow2(tt[1][2]) - pow2(tt[1][3]) - pow2(tt[2][3]) ) / 3. - 1./9.;
        double qCoefRt = sqrt( -qCoef);
        double rCoef = -0.5 * ( qCoef + 1./9. + tt[1][1] * pow2(tt[2][3]) + tt[2][2] * pow2(tt[1][3]) + tt[3][3] * pow2(tt[1][2]) - tt[1][1] * tt[2][2] * tt[3][3] ) + tt[1][2] * tt[1][3] * tt[2][3] + 1./27.;
        double pTemp = std::max( std::min( rCoef / pow3(qCoefRt), 1.), -1.);
        double pCoef = cos( acos(pTemp) / 3.);
        double pCoefRt = sqrt( 3. * (1. - pow2(pCoef)) );
        Sph_eVal1 = 1./3. + qCoefRt * std::max( 2. * pCoef,  pCoefRt - pCoef);
        Sph_eVal3 = 1./3. + qCoefRt * std::min( 2. * pCoef, -pCoefRt - pCoef);
        Sph_eVal2 = 1. - Sph_eVal1 - Sph_eVal3;
        
        // Sph: Begin find first and last eigenvector.
        for (int iVal=0;iVal<2;++iVal) {
            double eVal = (iVal == 0) ? Sph_eVal1 : Sph_eVal3;
            
            // If all particles are back-to-back then simpleminded third axis.
            //if (iVal > 0 && Sph_eVal2 < EIGENVALUEMIN) {
            if (iVal > 0 && Sph_eVal2 < 1.0e-10) {
                if ( std::abs(Sph_eVec1.pz()) > 0.5) {
                    Sph_eVec3 = Vec4( 1., 0., 0., 0.);
                }
                else{                           
                    Sph_eVec3 = Vec4( 0., 0., 1., 0.);
                }
                Sph_eVec3 -= dot3( Sph_eVec1, Sph_eVec3) * Sph_eVec1;
                Sph_eVec3 /= Sph_eVec3.pAbs();
                Sph_eVec2  = cross3( Sph_eVec1, Sph_eVec3);
                continue;
            }
            
            // Set up matrix to diagonalize.
            double dd[4][4];
            for (int j=1;j<4;++j) { 
                dd[j][j]=tt[j][j]-eVal;
                for(int k=j+1; k<4; ++k) {
                    dd[j][k]=tt[j][k]; dd[k][j]=tt[j][k];
                } 
            }
            
            // Find largest = pivotal element in matrix.
            int jMax = 0; int kMax = 0; double ddMax = 0.;
            for (int j = 1; j < 4; ++j) {
                for (int k = 1; k < 4; ++k) {
                    if (std::abs(dd[j][k]) > ddMax) {
                        jMax=j; 
                        kMax=k; 
                        ddMax=std::abs(dd[j][k]);
                    }
                }
            }
            
            // Subtract one row from the other two; find new largest element.
            int jMax2 = 0; ddMax = 0.;
            for(int j=1;j<4;++j){
                if (j != jMax) {
                    double pivot = dd[j][kMax] / dd[jMax][kMax];
                    for(int k=1;k<4;++k){
                        dd[j][k] -= pivot * dd[jMax][k];
                        if (std::abs(dd[j][k]) > ddMax){ jMax2=j; ddMax=std::abs(dd[j][k]); }
                    }
                }
            }
            
            // Construct eigenvector. Normalize to unit length; sign irrelevant.
            int k1 = kMax + 1; if (k1 > 3){k1 -= 3;}
            int k2 = kMax + 2; if (k2 > 3){k2 -= 3;}
            double eVec[4];
            eVec[k1]   = -dd[jMax2][k2];
            eVec[k2]   =  dd[jMax2][k1];
            eVec[kMax] = (dd[jMax][k1] * dd[jMax2][k2] - dd[jMax][k2] * dd[jMax2][k1]) / dd[jMax][kMax];

        //    std::cout << dd[1][1] << " " << dd[1][2] << " " << dd[1][3] << " " << dd[2][2] << " " << dd[2][3] << dd[3][3] << "\n";
        //  std::cout << k1 << " " << k2 << " " << kMax << " " << jMax << " " << jMax2 << "\n\n";


            double length = sqrt( pow2(eVec[1]) + pow2(eVec[2]) + pow2(eVec[3]) );
            
        //std::cout << eVec[1] << " : " << eVec[2] << " : " << eVec[3] << " length: " << length  <<  "\n";

            // Store eigenvectors.
            if (iVal == 0) Sph_eVec1 = Vec4( eVec[1] / length, eVec[2] / length, eVec[3] / length, 0.);
            else Sph_eVec3 = Vec4( eVec[1] / length, eVec[2] / length, eVec[3] / length, 0.);
        }
        
        // Sph: Middle eigenvector is orthogonal to the other two; sign irrelevant.
        Sph_eVec2 = cross3( Sph_eVec1, Sph_eVec3);

        
        //std::cout << Sph_eVec1.px() << " : " << Sph_eVec1.py() << " : " << Sph_eVec1.pz() <<  "\n";
        
    /*	  // Sph: Return info on results of analysis.
        double sphericity()      const {return 1.5 * (eVal2 + eVal3);}
        double aplanarity()      const {return 1.5 * eVal3;}
        double eigenValue(int i) const {return (i < 2) ? eVal1 : ( (i < 3) ? eVal2 : eVal3 ) ;}
        Vec4 eventAxis(int i)    const {return (i < 2) ? eVec1 : ( (i < 3) ? eVec2 : eVec3 ) ;}
    *//*	
        
        // can cut on sphericity here - if you don't want the event, maybe just write out a line with just a '0' below, or skip entirely...
        // if the file write is skipped, may want to skip it back on line 196 as well ( "fileout << "0\n";" )
        
        // Writing thrust to file (Thr_eVal1) - thrust normal vector is stored in eVec1
        // WARNING, thrust normal vector may contain +/-inf and/or +/-nan (when sphericity == 0.)!!!

        //..Sphericity Cut Here:
        //..Want sph axis (given by Sph_eVec1) to have a polar angle above 35 degrees.
        */
        double rho = sqrt(Sph_eVec1.px()*Sph_eVec1.px() + Sph_eVec1.py()*Sph_eVec1.py());
        double phi = atan(rho / Sph_eVec1.pz());
        phiDeg = phi * (180.0/3.1415926535897);
        
        if (std::abs(phiDeg) <= 35.0) {
            return error;
        }
    }

    //double sphericity = 1.5 * (Sph_eVal2 + Sph_eVal3);
    if(Thr_eVal1 == 0) Thr_eVal1 = -1;
    vector<double> output = {Thr_eVal1, phiDeg};
    return output;
}

//histogram to csv method
int histToCSV(TH1D* hist, string title){
    int n = hist->GetNbinsX();
    string xTitle = hist->GetXaxis()->GetTitle();
    string yTitle = hist->GetYaxis()->GetTitle();
    string filename = "QVir_Analysis/csv/" + title;
    
    //headers
    ofstream csvfile(filename.c_str());
    csvfile << xTitle << "," << yTitle << "," << "width" << endl;

    for(int i = 1; i <= n; i++){
        csvfile << hist->GetBinCenter(i) << "," << hist->GetBinContent(i) << "," << hist->GetBinWidth(i) << endl;
    }

    csvfile.close();
    return 0;
}

//compiles data for input to Bayesian code
void makeDatFile(vector<vector<double>> data, string title, string header){
    string filename = "QVir_Analysis/" + title + ".dat";
    cout << "Making " << filename << endl;
    ofstream datfile(filename);
    datfile << header << endl;

    int imax = data.at(0).size();
    int jmax = data.size();
    for(int i  = 0; i < imax; i++){
        for(int j = 0; j < jmax; j++) datfile << data[j][i] << " ";
        datfile << endl;
    }

    datfile.close();
}

//scales histogram and adjusts bin content for bin widths and centers since ROOT can only do width
void scaleBins(TH1D* hist, TProfile* prof, double scale = 1){
    double sum = 0;

    for(int i = 1; i <= hist->GetNbinsX(); i++){
        double raw = hist->GetBinContent(i);
        double center = prof->GetBinContent(i); if(center == 0 || isnan(center)) continue;
        double scaled = raw*scale/(hist->GetBinWidth(i)*center);
        //cout << raw << " " << scaled << " " << hist->GetBinWidth(i) << " " << scale << endl;
        hist->SetBinContent(i, scaled);
        sum += scaled*hist->GetBinWidth(i);
    }

    //cout << sum << endl;
}
//scales histogram and adjusts bin content for bin widths and centers since ROOT can only do width
void scaleBins(TH1D* hist, double scale = 1.0){
    double sum = 0;

    for(int i = 1; i <= hist->GetNbinsX(); i++){
        double raw = hist->GetBinContent(i);
        double center = hist->GetBinCenter(i); if(center==0 || isnan(center) || raw==0) continue;
        double scaled = raw*scale/(hist->GetBinWidth(i)*center);
        //cout << raw << " " << scaled << " " << hist->GetBinWidth(i) << " " << scale << endl;
        hist->SetBinContent(i, scaled);
        sum += scaled*hist->GetBinWidth(i);
    }

    //cout << sum << endl;
}

//to iterate over directories for q vir analysis
std::vector<std::string> get_directories(const std::string& s){
    std::vector<std::string> r;
    for(auto& p : boost::make_iterator_range(boost::filesystem::directory_iterator(s)))
        if (boost::filesystem::is_directory(p)){
            string dir = p.path().string();
            r.push_back(dir);
        }
    return r;
}

//getting directories for comparisons
std::vector<std::string> getComparisonDirs(int argc, char* argv[]){
    std::vector<std::string> directories;

    for(int i = 1; i < argc; i++){
        chdir(argv[i]);
        vector<string> tempdirs = get_directories(".");

        //removing the results dir
        int rmindex = 0;
        for(int j = 0; j < tempdirs.size(); j++){
            tempdirs[j].erase(0,2);
            if(tempdirs[j].find("QVir_Analysis") != string::npos) rmindex = j;
        }
        tempdirs.erase(tempdirs.begin() + rmindex); 

        vector<string> sorteddirs = doubleSort(tempdirs); //sorting for consistency
        for(int j = 0; j < sorteddirs.size(); j++){
            string temp = argv[i] + sorteddirs[j];
            sorteddirs[j] = temp;
        }
        directories.insert(directories.end(), sorteddirs.begin(), sorteddirs.end()); //inserting into the end of total vector
    }
    chdir("/scratch/user/cameron.parker/newJETSCAPE/JETSCAPE/build");

    return directories;
}

//gets ptHat bounds associated with dat files in a directories
std::vector<std::vector<string>> getDatBounds(const std::string& s){
    //reading the list of files
    std::vector<std::string> r;
    for(auto& p : boost::make_iterator_range(boost::filesystem::directory_iterator(s))){
        r.push_back(p.path().string());
        //cout << r.back() << endl; //debugging line
    }

    sort(r.begin(),r.end());

    //parsing the files for the seperate bounds and adding them to lists
    std::vector<string> lows;
    std::vector<string> highs;
    cout << "Dat files:" << endl;
    for(int i = 0; i < r.size(); i++){
        cout << "\t" << r[i] << endl;
        int first = r[i].find("Bin") + 3;
        int second = r[i].find("_",first) + 1;
        int end = r[i].find(".dat");

        string lower = r[i].substr(first,second-first-1);
        string upper = r[i].substr(second, end-second);

        lows.push_back(lower);
        highs.push_back(upper);
    }

    std::vector<string> sortedlows = doubleSort(lows);
    std::vector<string> sortedhighs = doubleSort(highs);
    sortedlows.erase(unique(sortedlows.begin(), sortedlows.end()), sortedlows.end());
    sortedhighs.erase(unique(sortedhighs.begin(), sortedhighs.end()), sortedhighs.end());

    //output
    std::vector<std::vector<string>> output; output.push_back(sortedlows); output.push_back(sortedhighs);
    return output;
}

//sorting string lists by their double values
std::vector<std::string> doubleSort(std::vector<std::string> input){
    //creating and filling vector of doubles
    std::vector<double> doublevec;
    for(int i  = 0; i < input.size(); i++) doublevec.push_back(stod(input[i]));
    sort(doublevec.begin(), doublevec.end());

    //matching sorted vec to strings
    std::vector<std::string> output;
    for(int i  = 0; i < doublevec.size(); i++){
        for(int j = 0; j < input.size(); j++){
            if(doublevec[i] == stod(input[j])){
                output.push_back(input[j]);
                //cout << doublevec[i] << " " << input[j] <<endl; //debugging line
                continue;
            }
        }
    }

    //returning sorted vector
    return output;
}

//smoothing spikes in identified hadrons
void smoothBinsOld(TH1D* hist){
    for(int i = 5; i < hist->GetNbinsX(); i++){
        if(hist->GetBinContent(i) > hist->GetBinContent(i-1) && hist->GetBinCenter(i) > 1){ //checking to see if its needed
            //cout << hist->GetBinCenter(i) << " " << hist->GetBinContent(i) << " " << hist->GetBinContent(i-1) << endl;
            //setting initial values for calculations
            double lower = hist->GetBinContent(i-1)*hist->GetBinWidth(i-1);
            double center = hist->GetBinContent(i)*hist->GetBinWidth(i);
            double higher = hist->GetBinContent(i+1)*hist->GetBinWidth(i+1);
            double del = hist->GetBinCenter(i+1) - hist->GetBinCenter(i-1);
            double strength = exp(-1); //how much is shifted out of the initial peak bin
            double shift = 0;
            if(lower > higher) shift = (lower-higher)/(higher+lower); //making sure we dont end up with negative values

            //apply slanted gaussian spread to bins higher than their predecessor
            double newcenter = center*(1-strength);
            double newlower = lower + center*strength*(1+shift)/2;
            double newhigher = higher + center*strength*(1-shift)/2;

            //setting new values for smoothed bins
            hist->SetBinContent(i-1, newlower/hist->GetBinWidth(i-1));
            hist->SetBinContent(i, newcenter/hist->GetBinWidth(i));
            hist->SetBinContent(i+1, newhigher/hist->GetBinWidth(i+1));
        }
    }
}

void smoothBins(TH1D* hist){
    for(int i = 5; i < hist->GetNbinsX(); i++){
        if(hist->GetBinCenter(i) > 1){ //checking to see if its needed
            //cout << hist->GetBinCenter(i) << " " << hist->GetBinContent(i) << " " << hist->GetBinContent(i-1) << endl;
            //setting initial values for calculations
            double lower = hist->GetBinContent(i-1)*hist->GetBinWidth(i-1);
            double center = hist->GetBinContent(i)*hist->GetBinWidth(i);
            double higher = hist->GetBinContent(i+1)*hist->GetBinWidth(i+1);
            double del = hist->GetBinCenter(i+1) - hist->GetBinCenter(i-1);
            double mid = hist->GetBinCenter(i+1) - hist->GetBinCenter(i);

            //strength of what we move out
            double higherweighted = hist->GetBinContent(i+1);
            double lowerweighted = hist->GetBinContent(i-1)-higherweighted;
            double centerweighted = hist->GetBinContent(i)-higherweighted;
            double strength = (centerweighted - (lowerweighted*mid)/(del))/(1.5*hist->GetBinContent(i)); //how much is shifted out of the initial peak bin
            if(lowerweighted < 0) strength = exp(-1); //default for when 2 bins in a row are higher than the lower one
            //strength = strength + strength*strength; //testing higher order corrections

            //apply slanted spread to bins higher than their predecessor
            double shift = 0;
            if(lower >= higher) shift = (lower-higher)/(higher+lower); //making sure we dont end up with negative values
            double newcenter = center*(1-strength);
            double newlower = 0, newhigher = 0;
            
            //cases to handle there being 2 bins in a row increasing
            /*if(lower > higher*1.1){
                newlower = lower + center*strength*(1+shift)/2;
                newhigher = higher + center*strength*(1-shift)/2;
            }else{
                newlower = lower;
                newhigher = higher + center*strength;    
            }*/
            newlower = lower + center*strength*(1+shift)/2;
            newhigher = higher + center*strength*(1-shift)/2;

            //setting new values for smoothed bins with filter for only events that need it
            if(strength > 0){
                hist->SetBinContent(i-1, newlower/hist->GetBinWidth(i-1));
                hist->SetBinContent(i, newcenter/hist->GetBinWidth(i));
                hist->SetBinContent(i+1, newhigher/hist->GetBinWidth(i+1));
            }
        }
    }
}

//Reading first line of the param file to get ECM
string getEcm(){
    std::ifstream infile("../QVir_Analysis/parameters.txt");
    string a,b;
    infile >> a >> b;
    infile.close();

    return b;
}

//Getting list of xsecs for a run based of the ECM in the param file and pthat bounds fed in
//Requires that the xsec file have the correct name of xsec[ECM].dat
vector<vector<double>> getXsecs(vector<string> pTHatMin, vector<string> pTHatMax){
    string Ecm = getEcm();
	vector<double> xsecList(pTHatMin.size(), 0);
    vector<double> xsecErrorList(pTHatMin.size(), 0);
    double xsectotal = 0;
    cout << "Center of Mass energy is: " << Ecm << endl;

    //reading file
    string line;
    ifstream myfile("/scratch/user/cameron.parker/pygen/xsec"+Ecm+".dat");
	if (myfile.is_open()){
        //go line by line and see if it matches any of the pT bins
		while ( getline (myfile,line) ){
			vector<string> result = splitString(line,", ");
			for(int i = 0; i < pTHatMin.size(); i++){
                //assiging values if there is a match
				if(stod(result[0]) == stod(pTHatMin[i]) && stod(result[1]) == stod(pTHatMax[i])){ 
					xsecList[i] = stod(result[2]);
                    xsecErrorList[i] = stod(result[3]);
					//cout << pTHatMin[i] << " " << pTHatMax[i] << " " << xsecList[i]*100000 << "\n"; //debugging line
                    break;
				}
			}
		}
		myfile.close();
	}

    for(int k = 0; k<xsecList.size(); k++) xsectotal += xsecList[k];

    //correcting total xsec to appropriate energy
    double actualXsec = 0;
    switch(stoi(Ecm)){
        case 2760:
            //https://arxiv.org/pdf/1208.4968.pdf
            actualXsec = 62.8;
            break;
        
        case 200:
            //https://arxiv.org/pdf/2005.00776.pdf
            actualXsec = 42.07;
            break;

        default:
            actualXsec = xsectotal;
            break;
    }
    cout << "Total Cross-section Corrected to: " << actualXsec << endl;
    for(int k = 0; k<xsecList.size(); k++) xsecList[k] *= actualXsec/xsectotal;

    //output
    vector<vector<double>> xsecout;
    xsecout.push_back(xsecList); xsecout.push_back(xsecErrorList);
    return xsecout;
}

//getting single xsec for a file
vector<double> getXsec(char HadronFile[300]){
    ifstream myfile(HadronFile);
    string line;
    double xsec = 0;
    double xsecerr = 0;

	if (myfile.is_open()){
        //go line by line and see if it matches the line for xsecs
		while( getline(myfile,line) ){
            std::size_t found = line.find("# JetScape writersigmaGen ");
            if(found!=std::string::npos){
                xsec = stod(line.substr(26));

                //getting err
                getline(myfile,line);
                xsecerr = stod(line.substr(26));
            }
        }
    }

    //output
    vector<double> output; output.push_back(xsec); output.push_back(xsecerr);
    return output;
}

//scrapes all dat files for the xsec and returns them
vector<vector<double>> getAllXsecs(vector<string> pTHatMin, vector<string> pTHatMax){
    string Ecm = getEcm();
	vector<double> xsecList;
    vector<double> xsecErrorList;
    double xsectotal = 0;
    char HadronFile[300];

    for(int i = 0; i < pTHatMin.size(); i++){
        sprintf(HadronFile,"dat/PP_Bin%s_%s.dat", pTHatMin[i].c_str(), pTHatMax[i].c_str());
        vector<double> temp = getXsec(HadronFile);
        xsecList.push_back(temp[0]);
        xsecErrorList.push_back(temp[1]);
    }

    //output
    vector<vector<double>> xsecout;
    xsecout.push_back(xsecList); xsecout.push_back(xsecErrorList);
    return xsecout;
}

//integral for root histograms
double integral(TH1D* hist){
    double total = 0;
    for(int i = 1; i <= hist->GetNbinsX(); i++) total += hist->GetBinContent(i);
    return total;
}

//blending for soft to hard QCD transition
double blending(int binID, double pT, double center, double width, double softcut){
    if(binID == 0) width *= 4; //increasing width for soft bin

    //defining blending region and value
    double lowlimit = center - width/2;
    double highlimit = center + width/2;
    double strength = pow((pT - lowlimit)/width,1);

    //Formatted in bin, then pt range, with blending region left for last
    if(pT < lowlimit) return 1;
    if(pT > highlimit) return 0;
    return 1 - strength;

    //catch all returning 0 in case of error
    return 0;
}

//subtracting floors out of azi correlation hists
TH1D* getZeroedHist(TH1D* hist){
    //cloning hist so we dont overwrite old files
    TH1D* histout = (TH1D*)hist->Clone();

    //subtracting out floor
    double floor = histout->GetMinimum();
    for(int i = 1; i<=histout->GetNbinsX(); i++) histout->Fill(histout->GetBinCenter(i),-1.0*floor);

    //returning output
    return histout;
}
