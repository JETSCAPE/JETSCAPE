
#include "ThermPtnSampler.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <unordered_map>
#include <vector>

#include "JetScapeConstants.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"

using namespace Jetscape;

ThermalPartonSampler::ThermalPartonSampler(unsigned int ran_seed) {
  rng_engine.seed(ran_seed);  // set the random seed from the framework

  // Gaussian weights for integration
  GWeight[0] = 0.0621766166553473;
  GWeight[1] = 0.0619360674206832;
  GWeight[2] = 0.0614558995903167;
  GWeight[3] = 0.0607379708417702;
  GWeight[4] = 0.0597850587042655;
  GWeight[5] = 0.0586008498132224;
  GWeight[6] = 0.0571899256477284;
  GWeight[7] = 0.0555577448062125;
  GWeight[8] = 0.0537106218889962;
  GWeight[9] = 0.0516557030695811;
  GWeight[10] = 0.0494009384494663;
  GWeight[11] = 0.0469550513039484;
  GWeight[12] = 0.0443275043388033;
  GWeight[13] = 0.0415284630901477;
  GWeight[14] = 0.0385687566125877;
  GWeight[15] = 0.0354598356151462;
  GWeight[16] = 0.0322137282235780;
  GWeight[17] = 0.0288429935805352;
  GWeight[18] = 0.0253606735700124;
  GWeight[19] = 0.0217802431701248;
  GWeight[20] = 0.0181155607134894;
  GWeight[21] = 0.0143808227614856;
  GWeight[22] = 0.0105905483836510;
  GWeight[23] = 0.0067597991957454;
  GWeight[24] = 0.0029086225531551;
  GWeight[25] = 0.0621766166553473;
  GWeight[26] = 0.0619360674206832;
  GWeight[27] = 0.0614558995903167;
  GWeight[28] = 0.0607379708417702;
  GWeight[29] = 0.0597850587042655;
  GWeight[30] = 0.0586008498132224;
  GWeight[31] = 0.0571899256477284;
  GWeight[32] = 0.0555577448062125;
  GWeight[33] = 0.0537106218889962;
  GWeight[34] = 0.0516557030695811;
  GWeight[35] = 0.0494009384494663;
  GWeight[36] = 0.0469550513039484;
  GWeight[37] = 0.0443275043388033;
  GWeight[38] = 0.0415284630901477;
  GWeight[39] = 0.0385687566125877;
  GWeight[40] = 0.0354598356151462;
  GWeight[41] = 0.0322137282235780;
  GWeight[42] = 0.0288429935805352;
  GWeight[43] = 0.0253606735700124;
  GWeight[44] = 0.0217802431701248;
  GWeight[45] = 0.0181155607134894;
  GWeight[46] = 0.0143808227614856;
  GWeight[47] = 0.0105905483836510;
  GWeight[48] = 0.0067597991957454;
  GWeight[49] = 0.0029086225531551;

  // Gaussian abscissas for integration
  GAbs[0] = 0.0310983383271889;
  GAbs[1] = 0.0931747015600861;
  GAbs[2] = 0.1548905899981459;
  GAbs[3] = 0.2160072368760418;
  GAbs[4] = 0.2762881937795320;
  GAbs[5] = 0.3355002454194373;
  GAbs[6] = 0.3934143118975651;
  GAbs[7] = 0.4498063349740388;
  GAbs[8] = 0.5044581449074642;
  GAbs[9] = 0.5571583045146501;
  GAbs[10] = 0.6077029271849502;
  GAbs[11] = 0.6558964656854394;
  GAbs[12] = 0.7015524687068222;
  GAbs[13] = 0.7444943022260685;
  GAbs[14] = 0.7845558329003993;
  GAbs[15] = 0.8215820708593360;
  GAbs[16] = 0.8554297694299461;
  GAbs[17] = 0.8859679795236131;
  GAbs[18] = 0.9130785566557919;
  GAbs[19] = 0.9366566189448780;
  GAbs[20] = 0.9566109552428079;
  GAbs[21] = 0.9728643851066920;
  GAbs[22] = 0.9853540840480058;
  GAbs[23] = 0.9853540840480058;
  GAbs[24] = 0.9988664044200710;
  GAbs[25] = -0.0310983383271889;
  GAbs[26] = -0.0931747015600861;
  GAbs[27] = -0.1548905899981459;
  GAbs[28] = -0.2160072368760418;
  GAbs[29] = -0.2762881937795320;
  GAbs[30] = -0.3355002454194373;
  GAbs[31] = -0.3934143118975651;
  GAbs[32] = -0.4498063349740388;
  GAbs[33] = -0.5044581449074642;
  GAbs[34] = -0.5571583045146501;
  GAbs[35] = -0.6077029271849502;
  GAbs[36] = -0.6558964656854394;
  GAbs[37] = -0.7015524687068222;
  GAbs[38] = -0.7444943022260685;
  GAbs[39] = -0.7845558329003993;
  GAbs[40] = -0.8215820708593360;
  GAbs[41] = -0.8554297694299461;
  GAbs[42] = -0.8859679795236131;
  GAbs[43] = -0.9130785566557919;
  GAbs[44] = -0.9366566189448780;
  GAbs[45] = -0.9566109552428079;
  GAbs[46] = -0.9728643851066920;
  GAbs[47] = -0.9853540840480058;
  GAbs[48] = -0.9853540840480058;
  GAbs[49] = -0.9988664044200710;

  // Adjustable parameters
  xmq = 0.33 / GEVFM;  // use pythia values here
  xms = 0.5 / GEVFM;   // use pythia values here
  T = 0.165 / GEVFM;   // 165 MeV into fm^-1 // is set from outside with
                       // brick_temperature() or from hypersurface
  NUMSTEP = 1048577;   // 2^20+1, for steps of CDF Table, changes coarseness of
                       // momentum sampling

  // Adjustable params for 3+1d
  CellDX = 0.2;
  CellDY = 0.2;
  CellDZ = 0.2;
  CellDT = 0.1;  // these are the values used in SurfaceFinder.cc

  // Flags
  SetNum = false;         // Set 'true' to set number of particles by hand-
                          // !!!Statistics use above temperature!!!
  SetNumLight = 1000000;  // If SetNum == true, this many UD quarks are
                          // generated
  SetNumStrange = 0;   // If SetNum == true, this many S quarks are generated
  ShuffleList = true;  // Should list of particles be shuffled at the end

  // Brick Info
  L = 4.0;     // Thickness from box edge
  W = 4.0;     // x/y width of box
  Time = 2.0;  // Time of the brick partons

  Vx = 0.;  // 'Uniform' no flow in x-dir
  Vy = 0.;  // 'Uniform' no flow in y-dir
  Vz = 0.;  // 'Uniform' no flow in z-dir

  surface.clear();

  num_ud = 0;
  num_s = 0;

  for (int iCDF = 0; iCDF < NUMSTEP; ++iCDF) {
    std::vector<double> temp = {0., 0.};
    CDFTabLight.push_back(temp);
    CDFTabStrange.push_back(temp);
  }
}

double ThermalPartonSampler::FermiPDF(double P, double M, double T, double mu) {
  return 1. / (exp((sqrt(M * M + P * P) - mu) / T) + 1.);
}

bool ThermalPartonSampler::SplitSort(double goal, int floor, int ceiling,
                                     int quark) {
  std::vector<std::vector<double>>& CDFTab =
      (quark == 1) ? CDFTabLight : CDFTabStrange;

  int TargetPoint = ((floor + ceiling) / 2);
  double TargetVal = CDFTab[TargetPoint][1];

  return (goal > TargetVal);
}

void ThermalPartonSampler::CDFGenerator(double Temp, double M, int quark) {
  double PStep;
  double PMax = 10. * Temp;      // CutOff for Integration
  PStep = PMax / (NUMSTEP - 1);  // Stepsize in P

  std::vector<std::vector<double>> CDFTab;
  if (quark == 1) {
    CDFTab = CDFTabLight;
  } else if (quark == 2) {
    CDFTab = CDFTabStrange;
  } else {
    JSWARN << "This is not a valid quark input for CDFGenerator()";
    return;
  }

  // Initialize tabulated results for CDF
  CDFTab[0][0] = 0;  // For zero momentum or less...
  CDFTab[0][1] = 0;  // There is zero chance

  // Precompute constant values used in the loop
  double muPi0 = 0.0;  // You need to provide the correct value of muPi0

  // Tabulate CDF(x) = int(0->x) PDF
  for (int i = 1; i < NUMSTEP; i++) {
    CDFTab[i][0] =
        CDFTab[i - 1][0] + PStep;  // Calculate Momentum of the next step
    double Fermi0 = FermiPDF(CDFTab[i - 1][0], M, Temp, muPi0);  // PDF
    double Fermi1 = FermiPDF(CDFTab[i][0], M, Temp, muPi0);
    CDFTab[i][1] = CDFTab[i - 1][1] +
                   (PStep / 2) * (Fermi0 * CDFTab[i - 1][0] * CDFTab[i - 1][0] +
                                  Fermi1 * CDFTab[i][0] * CDFTab[i][0]);
  }

  // Normalize the CDF
  for (int i = 0; i < NUMSTEP; i++) {
    CDFTab[i][1] = CDFTab[i][1] / CDFTab[NUMSTEP - 1][1];
  }

  if (quark == 1) {
    CDFTabLight = CDFTab;
  } else if (quark == 2) {
    CDFTabStrange = CDFTab;
  }
}

void ThermalPartonSampler::MCSampler(double Temp,
                                     int quark) {  // input temperature in fm^-1
  double PMag, CosT, Phi, PStep;
  double PMax = 10. * Temp;      // CutOff for Integration
  PStep = PMax / (NUMSTEP - 1);  // Stepsize in P

  std::vector<std::vector<double>>& CDFTab =
      (quark == 1) ? CDFTabLight : CDFTabStrange;

  int floor = 0;
  int ceiling = NUMSTEP - 1;
  double PRoll;
  double denominator;

  bool sample = true;
  int test = 0;
  while (sample) {
    PRoll = ran();

    for (int i = 0; i < 25; i++) {  // Use 25 iterations for both quark types
      if (SplitSort(PRoll, floor, ceiling, quark)) {
        floor = ((floor + ceiling) / 2);
      } else {
        ceiling = ((floor + ceiling) / 2);
      }
    }

    denominator = CDFTab[ceiling][1] - CDFTab[floor][1];
    if (std::fabs(denominator) > 1e-16) {
      PMag =
          PStep * (PRoll - CDFTab[floor][1]) / denominator + CDFTab[floor][0];
      sample = false;
    }
  }

  CosT = (ran() - 0.5) * 2.0;
  Phi = ran() * 2 * PI;

  NewX = PMag * sqrt(1 - CosT * CosT) * cos(Phi);
  NewY = PMag * sqrt(1 - CosT * CosT) * sin(Phi);
  NewZ = PMag * CosT;
  NewP = PMag;
}

void ThermalPartonSampler::samplebrick() {
  // preliminary parameter checks
  if (L < 0.) {
    L = -L;
    JSWARN << "Negative brick length - setting to positive " << L << " fm.";
  }
  if (W < 0.) {
    W = -W;
    JSWARN << "Negative brick width - setting to positive " << W << " fm.";
  }

  // Input read from cells
  double LFSigma[4];  // LabFrame hypersurface (tau/t,x,y,eta/z=0)
  double CMSigma[4];  // Center of mass hypersurface (tau/t,x,y,eta/z=0)
  double Vel[4];      // Gamma & 3-velocity of cell (gamma, Vx, Vy, Vz) NOT FOUR
                      // VELOCITY

  // Calculated global quantities
  int PartCount;            // Total Count of Particles over ALL cells
  double NumLight;          // Number DENSITY of light quarks at set T
  double NumStrange;        // Number DENSITY of s quarks at set T
  double UDDeg = 4. * 6.;   // Degeneracy of UD quarks
  double OddDeg = 2. * 6.;  // Degeneracy of S quarks

  // counter for total number of light and strange quarks
  int nL_tot = 0;
  int nS_tot = 0;

  // Calculated quantities in cells
  double LorBoost[4][4];    // Lorentz boost defined as used - form is always
                            // Lambda_u^v
  int GeneratedParticles;   // Number of particles to be generated this cell
  double new_quark_energy;  // store new quark energy for boost

  // End definition of static variables

  // Define hypersurface for brick here
  // Default t=const hypersurface
  // Vector must be covariant (negative signs on spatial components)
  LFSigma[0] = 1.;
  LFSigma[1] = 0.;
  LFSigma[2] = 0.;
  LFSigma[3] = 0.;

  Vel[1] = Vx;
  Vel[2] = Vy;
  Vel[3] = Vz;

  double vsquare = Vel[1] * Vel[1] + Vel[2] * Vel[2] + Vel[3] * Vel[3];

  if (vsquare > 1.) {
    JSWARN << "v^2 = " << vsquare;
    JSWARN << "Unphysical velocity (brick flow)! Set to \"No Flow\" case";
    Vel[1] = 0.;
    Vel[2] = 0.;
    Vel[3] = 0.;
  }

  Vel[0] = 1. / std::sqrt(1. - vsquare);  // gamma - Vel is not four velocity

  // Lambda_u ^v from (lab frame to rest frame with flow velocity)
  if (vsquare == 0) {
    LorBoost[0][0] = Vel[0];
    LorBoost[0][1] = Vel[0] * Vel[1];
    LorBoost[0][2] = Vel[0] * Vel[2];
    LorBoost[0][3] = Vel[0] * Vel[3];
    LorBoost[1][0] = Vel[0] * Vel[1];
    LorBoost[1][1] = 1.;
    LorBoost[1][2] = 0.;
    LorBoost[1][3] = 0.;
    LorBoost[2][0] = Vel[0] * Vel[2];
    LorBoost[2][1] = 0.;
    LorBoost[2][2] = 1.;
    LorBoost[2][3] = 0.;
    LorBoost[3][0] = Vel[0] * Vel[3];
    LorBoost[3][1] = 0.;
    LorBoost[3][2] = 0.;
    LorBoost[3][3] = 1.;
  } else {
    LorBoost[0][0] = Vel[0];
    LorBoost[0][1] = Vel[0] * Vel[1];
    LorBoost[0][2] = Vel[0] * Vel[2];
    LorBoost[0][3] = Vel[0] * Vel[3];
    LorBoost[1][0] = Vel[0] * Vel[1];
    LorBoost[1][1] = (Vel[0] - 1.) * Vel[1] * Vel[1] / vsquare + 1.;
    LorBoost[1][2] = (Vel[0] - 1.) * Vel[1] * Vel[2] / vsquare;
    LorBoost[1][3] = (Vel[0] - 1.) * Vel[1] * Vel[3] / vsquare;
    LorBoost[2][0] = Vel[0] * Vel[2];
    LorBoost[2][1] = (Vel[0] - 1.) * Vel[1] * Vel[2] / vsquare;
    LorBoost[2][2] = (Vel[0] - 1.) * Vel[2] * Vel[2] / vsquare + 1.;
    LorBoost[2][3] = (Vel[0] - 1.) * Vel[2] * Vel[3] / vsquare;
    LorBoost[3][0] = Vel[0] * Vel[3];
    LorBoost[3][1] = (Vel[0] - 1.) * Vel[1] * Vel[3] / vsquare;
    LorBoost[3][2] = (Vel[0] - 1.) * Vel[2] * Vel[3] / vsquare;
    LorBoost[3][3] = (Vel[0] - 1.) * Vel[3] * Vel[3] / vsquare + 1.;
  }

  // Code parity with hypersurface case
  // Lambda_u^v Sigma_v = CMSigma_u
  // Caution: CMSigma is a covariant vector
  CMSigma[0] = (LorBoost[0][0] * LFSigma[0] + LorBoost[0][1] * LFSigma[1] +
                LorBoost[0][2] * LFSigma[2] + LorBoost[0][3] * LFSigma[3]);
  CMSigma[1] = (LorBoost[1][0] * LFSigma[0] + LorBoost[1][1] * LFSigma[1] +
                LorBoost[1][2] * LFSigma[2] + LorBoost[1][3] * LFSigma[3]);
  CMSigma[2] = (LorBoost[2][0] * LFSigma[0] + LorBoost[2][1] * LFSigma[1] +
                LorBoost[2][2] * LFSigma[2] + LorBoost[2][3] * LFSigma[3]);
  CMSigma[3] = (LorBoost[3][0] * LFSigma[0] + LorBoost[3][1] * LFSigma[1] +
                LorBoost[3][2] * LFSigma[2] + LorBoost[3][3] * LFSigma[3]);

  /* Define Parton Densities */

  double cut = 10. * T;   // Each coordinate of P is integrated this far
  double E_light;         // Store energy of light quarks
  double E_strange;       // Store energy of strange quarks
  double GWeightProd;     // Needed for integration
  double pSpatialdSigma;  // Needed for integration
  PartCount = 0;
  NumLight = 0;  // Initialize density for light and strange quarks
  NumStrange = 0;

  // Define the lambda function to compute the energy
  auto computeEnergy = [](double mass, double cut, double GAbsL, double GAbsM,
                          double GAbsK) {
    return sqrt(mass * mass + (cut * GAbsL) * (cut * GAbsL) +
                (cut * GAbsM) * (cut * GAbsM) + (cut * GAbsK) * (cut * GAbsK));
  };

  // Define the lambda function for the Fermi distribution
  auto FermiDistribution = [](double energy, double temperature) {
    return 1. / (exp(energy / temperature) + 1.);
  };

  // GAUSSIAN INTEGRALS <n> = int f(p)d3p
  for (int l = 0; l < GPoints; l++) {
    for (int m = 0; m < GPoints; m++) {
      for (int k = 0; k < GPoints; k++) {
        GWeightProd = GWeight[l] * GWeight[m] * GWeight[k];
        pSpatialdSigma = (cut * GAbs[l]) * CMSigma[1] +
                         (cut * GAbs[m]) * CMSigma[2] +
                         (cut * GAbs[k]) * CMSigma[3];

        // For UD quarks
        E_light = computeEnergy(xmq, cut, GAbs[l], GAbs[m], GAbs[k]);
        NumLight += GWeightProd * FermiDistribution(E_light, T) *
                    (E_light * CMSigma[0] + pSpatialdSigma) / E_light;

        // For S quarks
        E_strange = computeEnergy(xms, cut, GAbs[l], GAbs[m], GAbs[k]);
        NumStrange += GWeightProd * FermiDistribution(E_strange, T) *
                      (E_strange * CMSigma[0] + pSpatialdSigma) / E_strange;
      }
    }
  }

  // Normalization factors: degeneracy, gaussian integration, and 1/(2Pi)^3
  NumLight *= UDDeg * cut * cut * cut / (8. * PI * PI * PI);
  NumStrange *= OddDeg * cut * cut * cut / (8. * PI * PI * PI);

  // Define rest-frame cumulative functions
  CDFGenerator(T, xmq, 1);  // for light quarks
  CDFGenerator(T, xms, 2);  // for strange quarks

  // U, D, UBAR, DBAR QUARKS
  // <N> = V <n>
  double NumHere = NumLight * L * W * W;

  // S, SBAR QUARKS
  // <N> = V <n>
  double NumOddHere = NumStrange * L * W * W;

  std::poisson_distribution<int> poisson_ud(NumHere);
  std::poisson_distribution<int> poisson_s(NumOddHere);

  // In case of overwriting
  if (SetNum) {
    NumHere = SetNumLight;
    NumOddHere = SetNumStrange;
  }

  // Generating light quarks
  GeneratedParticles = poisson_ud(getRandomGenerator());

  // List of particles ( pos(x,y,z,t), mom(px,py,pz,E), species)
  for (int partic = 0; partic < GeneratedParticles; partic++) {
    // adding space to PList for output quarks
    std::vector<double> temp = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
    Plist.push_back(temp);

    // Species - U,D,UBar,Dbar are equally likely
    double SpecRoll = ran();  // Probability of species die roll
    if (SpecRoll <= 0.25) {   // UBar
      Plist[PartCount][1] = -2;
    } else if (SpecRoll <= 0.50) {  // DBar
      Plist[PartCount][1] = -1;
    } else if (SpecRoll <= 0.75) {  // D
      Plist[PartCount][1] = 1;
    } else {  // U
      Plist[PartCount][1] = 2;
    }

    // Position
    // Located at x,y pos of area element
    double XRoll = ran() - 0.5;  // center at x=0
    double YRoll = ran() - 0.5;  // center at y=0
    double ZRoll = ran() - 0.5;  // center at z=0

    Plist[PartCount][7] = XRoll * L;
    Plist[PartCount][8] = YRoll * W;
    Plist[PartCount][9] = ZRoll * W;

    // Time
    Plist[PartCount][10] = Time;  // Tau = L/2.: assume jet at light speed

    // Momentum
    // Sample rest frame momentum given T and mass of light quark
    MCSampler(T, 1);  // NewP=P, NewX=Px, ...

    // PLab^u = g^u^t Lambda_t ^v Pres^w g_w _v  = Lambda ^u _w Pres_w (with
    // velocity -v) (Lambda _u ^t with velocity v) == (Lambda ^u _t with
    // velocity -v) This is boost as if particle sampled in rest frame
    new_quark_energy = sqrt(xmq * xmq + NewP * NewP);
    Plist[PartCount][6] =
        (LorBoost[0][0] * new_quark_energy + LorBoost[0][1] * NewX +
         LorBoost[0][2] * NewY + LorBoost[0][3] * NewZ) *
        GEVFM;
    Plist[PartCount][3] =
        (LorBoost[1][0] * new_quark_energy + LorBoost[1][1] * NewX +
         LorBoost[1][2] * NewY + LorBoost[1][3] * NewZ) *
        GEVFM;
    Plist[PartCount][4] =
        (LorBoost[2][0] * new_quark_energy + LorBoost[2][1] * NewX +
         LorBoost[2][2] * NewY + LorBoost[2][3] * NewZ) *
        GEVFM;
    Plist[PartCount][5] =
        (LorBoost[3][0] * new_quark_energy + LorBoost[3][1] * NewX +
         LorBoost[3][2] * NewY + LorBoost[3][3] * NewZ) *
        GEVFM;
    // Additional information
    Plist[PartCount][0] = 1;   // Event ID, to match jet formatting
    Plist[PartCount][2] = 0;   // Origin, to match jet formatting
    Plist[PartCount][11] = 0;  // Status - identifies as thermal quark
    PartCount++;
  }

  nL_tot += GeneratedParticles;

  // Generate strange quarks
  GeneratedParticles = poisson_s(getRandomGenerator());

  nS_tot += GeneratedParticles;
  // List of particles ( pos(x,y,z,t), mom(px,py,pz,E), species)
  for (int partic = 0; partic < GeneratedParticles; partic++) {
    // adding space to PList for output quarks
    std::vector<double> temp = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
    Plist.push_back(temp);

    // Species - S,Sbar are equally likely
    double SpecRoll = ran();
    if (SpecRoll <= 0.5) {  // SBar
      Plist[PartCount][1] = -3;
    } else {  // S
      Plist[PartCount][1] = 3;
    }

    // Position
    // Located at x,y pos of area element
    double XRoll = ran() - 0.5;  // center at x=0
    double YRoll = ran() - 0.5;  // center at y=0
    double ZRoll = ran() - 0.5;  // center at z=0

    Plist[PartCount][7] = XRoll * L;
    Plist[PartCount][8] = YRoll * W;
    Plist[PartCount][9] = ZRoll * W;

    // Time
    Plist[PartCount][10] = Time;  // Tau = L/2.: assume jet at light speed

    // Momentum
    // Sample rest frame momentum given T and mass of s quark
    MCSampler(T, 2);  // NewP=P, NewX=Px, ...

    // PLab^u = g^u^t Lambda_t ^v Pres^w g_w _v  = Lambda ^u _w Pres_w (with
    // velocity -v) (Lambda _u ^t with velocity v) == (Lambda ^u _t with
    // velocity -v) This is boost as if particle sampled in rest frame
    new_quark_energy = sqrt(xms * xms + NewP * NewP);
    Plist[PartCount][6] =
        (LorBoost[0][0] * new_quark_energy + LorBoost[0][1] * NewX +
         LorBoost[0][2] * NewY + LorBoost[0][3] * NewZ) *
        GEVFM;
    Plist[PartCount][3] =
        (LorBoost[1][0] * new_quark_energy + LorBoost[1][1] * NewX +
         LorBoost[1][2] * NewY + LorBoost[1][3] * NewZ) *
        GEVFM;
    Plist[PartCount][4] =
        (LorBoost[2][0] * new_quark_energy + LorBoost[2][1] * NewX +
         LorBoost[2][2] * NewY + LorBoost[2][3] * NewZ) *
        GEVFM;
    Plist[PartCount][5] =
        (LorBoost[3][0] * new_quark_energy + LorBoost[3][1] * NewX +
         LorBoost[3][2] * NewY + LorBoost[3][3] * NewZ) *
        GEVFM;

    // Additional information
    Plist[PartCount][0] = 1;   // Event ID, to match jet formatting
    Plist[PartCount][2] = 0;   // Origin, to match jet formatting
    Plist[PartCount][11] = 0;  // Status - identifies as thermal quark
    PartCount++;
  }

  JSDEBUG << "Light particles: " << nL_tot;
  JSDEBUG << "Strange particles: " << nS_tot;

  num_ud = nL_tot;
  num_s = nS_tot;

  // Shuffling PList
  if (ShuffleList) {
    // Shuffle the Plist using the random engine
    std::shuffle(&Plist[0], &Plist[PartCount], getRandomGenerator());
  }

  // print Plist for testing
  /*std::cout << std::setprecision(5);
  std::ofstream thermalP;
  thermalP.open("thermal_partons.dat", std::ios::app);
  for(int p=0; p < Plist.size(); p++){
          thermalP << Plist[p][0] << " " << Plist[p][1] << " " << Plist[p][2]
          << " " << Plist[p][3] << " " << Plist[p][4] << " " << Plist[p][5]
          << " " << Plist[p][6] << " " << Plist[p][7] << " " << Plist[p][8]
          << " " << Plist[p][9] << " " << Plist[p][10] << " " << Plist[p][11]
          << "\n";
  }
  thermalP.close();*/
}

void ThermalPartonSampler::sample_3p1d(bool Cartesian_hydro) {
  // Input read from cells
  double CPos[4];  // Position of the current cell (tau/t, x, y , eta/z=0)
  double
      LFSigma[4];  // LabFrame hypersurface (tau/t,x,y,eta/z=0), expect Sigma_mu
  double CMSigma[4];  // Center of mass hypersurface (tau/t,x,y,eta/z=0)
  double TRead;       // Temperature of cell
  double Vel[4];      // Gamma & 3-Velocity of cell (gamma, Vx, Vy, Vz) NOT FOUR
                      // VELOCITY
  double tau_pos;     // proper time from position of cell
  double eta_pos;     // eta from position of cell
  double tau_sur;     // proper time from normal vector of surface
  double eta_sur;     // eta from normal vector of surface

  // Calculated global quantities
  int PartCount;            // Total count of particles over ALL cells
  double NumLight;          // Number DENSITY of light quarks at set T
  double NumStrange;        // Number DENSITY of squarks at set T
  double UDDeg = 4. * 6.;   // Degeneracy of UD quarks
  double OddDeg = 2. * 6.;  // Degeneracy of S quarks

  // Calculated quantities in cells
  double LorBoost[4][4];   // Lorentz boost defined as used - Form is always
                           // Lambda_u^v
  int GeneratedParticles;  // Number of particles to be generated this cell

  // Define parton densities
  double cut;             // Each coordinate of P is integrated this far
  double E_light;         // Store energy of light quarks
  double E_strange;       // Store energy of strange quarks
  double GWeightProd;     // Needed for integration
  double pSpatialdSigma;  // Needed for integration
  PartCount = 0;
  NumLight = 0;  // Initialize density for light and strange quarks
  NumStrange = 0;
  GeneratedParticles = 0;
  double new_quark_energy;  // store new quark energy for boost

  // counter for total number of light and strange quarks
  int nL_tot = 0;
  int nS_tot = 0;

  // define the accuracy range for temperature values in which the CDFGenerator
  // is executed (for the case the value is not in the cache yet)
  const double accuracyRange =
      0.003 / GEVFM;  // for a usual hypersurface at const temperature there
                      // should be only one entry in the cache
  // precompute CDFs and store them in a cache
  std::unordered_map<double, std::vector<std::vector<double>>> cdfLightCache;
  std::unordered_map<double, std::vector<std::vector<double>>> cdfStrangeCache;

  // lambda function to check if a temperature value is within the accuracy
  // range of a cached temperature
  auto withinAccuracyRange = [accuracyRange](double targetTemp,
                                             double cachedTemp) {
    return std::fabs(targetTemp - cachedTemp) <= accuracyRange;
  };

  // lambda function to retrieve the closest temperature from the cache within
  // the accuracy range
  auto getClosestCachedTemp =
      [](const std::unordered_map<double, std::vector<std::vector<double>>>&
             cache,
         double targetTemp) {
        double closestTemp = -1.0;
        double minDistance = std::numeric_limits<double>::max();
        for (const auto& entry : cache) {
          double cachedTemp = entry.first;
          double distance = std::fabs(targetTemp - cachedTemp);
          if (distance < minDistance) {
            minDistance = distance;
            closestTemp = cachedTemp;
          }
        }
        return closestTemp;
      };

  // Define the lambda function to compute the energy
  auto computeEnergy = [](double mass, double cut, double GAbsL, double GAbsM,
                          double GAbsK) {
    return sqrt(mass * mass + (cut * GAbsL) * (cut * GAbsL) +
                (cut * GAbsM) * (cut * GAbsM) + (cut * GAbsK) * (cut * GAbsK));
  };

  // Define the lambda function for the Fermi distribution
  auto FermiDistribution = [](double energy, double temperature) {
    return 1. / (exp(energy / temperature) + 1.);
  };

  for (int iS = 0; iS < surface.size(); ++iS) {
    tau_pos = surface[iS][0];
    CPos[1] = surface[iS][1];
    CPos[2] = surface[iS][2];
    eta_pos = surface[iS][3];  // we need t,x,y,z
    // this is also tau, x, y, eta
    tau_sur = surface[iS][4];
    LFSigma[1] = surface[iS][5];
    LFSigma[2] = surface[iS][6];
    eta_sur = surface[iS][7];
    TRead = surface[iS][8] / GEVFM;
    Vel[1] = surface[iS][9];
    Vel[2] = surface[iS][10];
    Vel[3] = surface[iS][11];

    cut = 10 * TRead;

    // check if the CDFs for light quarks for this temperature are already in
    // the cache and within the accuracy range
    auto cdfLightIter = cdfLightCache.find(TRead);
    if (cdfLightIter == cdfLightCache.end()) {
      // if not found, find the closest temperature in the cache within the
      // accuracy range
      double closestTemp = getClosestCachedTemp(cdfLightCache, TRead);

      if (closestTemp >= 0. && withinAccuracyRange(TRead, closestTemp)) {
        // use the CDFs from the closest temperature in the cache
        CDFTabLight = cdfLightCache[closestTemp];
        // set the current temperature value to the closest cached temperature
        TRead = closestTemp;
      } else {
        // if no suitable entry is found, compute the CDFs and store them in the
        // cache
        CDFGenerator(TRead, xmq, 1);  // for light quarks
        cdfLightCache[TRead] = CDFTabLight;
      }
    } else {
      // if found, use the precomputed CDFs from the cache
      CDFTabLight = cdfLightIter->second;
    }

    // check if the CDFs for strange quarks for this temperature are already in
    // the cache and within the accuracy range
    auto cdfStrangeIter = cdfStrangeCache.find(TRead);
    if (cdfStrangeIter == cdfStrangeCache.end()) {
      // if not found, find the closest temperature in the cache within the
      // accuracy range
      double closestTemp = getClosestCachedTemp(cdfStrangeCache, TRead);

      if (closestTemp >= 0 && withinAccuracyRange(TRead, closestTemp)) {
        // use the CDFs from the closest temperature in the cache
        CDFTabStrange = cdfStrangeCache[closestTemp];
        // set the current temperature value to the closest cached temperature
        TRead = closestTemp;
      } else {
        // if no suitable entry is found, compute the CDFs and store them in the
        // cache
        CDFGenerator(TRead, xms, 2);  // for strange quarks
        cdfStrangeCache[TRead] = CDFTabStrange;
      }
    } else {
      // if found, use the precomputed CDFs from the cache
      CDFTabStrange = cdfStrangeIter->second;
    }

    if (Cartesian_hydro == false) {
      // getting t,z from tau, eta
      CPos[0] = tau_pos * std::cosh(eta_pos);
      CPos[3] = tau_pos * std::sinh(eta_pos);
      // transform surface vector to Minkowski coordinates
      double cosh_eta_pos = std::cosh(eta_pos);
      double sinh_eta_pos = std::sinh(eta_pos);
      LFSigma[0] = cosh_eta_pos * tau_sur - (sinh_eta_pos / tau_pos) * eta_sur;
      LFSigma[3] = -sinh_eta_pos * tau_sur + (cosh_eta_pos / tau_pos) * eta_sur;
    } else {  // check later!
      CPos[0] = tau_pos;
      CPos[3] = eta_pos;
      LFSigma[0] = tau_sur;
      LFSigma[3] = eta_sur;
    }

    double vsquare = Vel[1] * Vel[1] + Vel[2] * Vel[2] + Vel[3] * Vel[3];
    if (vsquare < 10e-16) {
      vsquare = 10e-16;
    }

    // Deduced info
    Vel[0] = 1. / sqrt(1 - vsquare);  // gamma - Vel is not four velocity

    // Lambda_u ^v
    LorBoost[0][0] = Vel[0];
    LorBoost[0][1] = Vel[0] * Vel[1];
    LorBoost[0][2] = Vel[0] * Vel[2];
    LorBoost[0][3] = Vel[0] * Vel[3];
    LorBoost[1][0] = Vel[0] * Vel[1];
    LorBoost[1][1] = (Vel[0] - 1.) * Vel[1] * Vel[1] / vsquare + 1.;
    LorBoost[1][2] = (Vel[0] - 1.) * Vel[1] * Vel[2] / vsquare;
    LorBoost[1][3] = (Vel[0] - 1.) * Vel[1] * Vel[3] / vsquare;
    LorBoost[2][0] = Vel[0] * Vel[2];
    LorBoost[2][1] = (Vel[0] - 1.) * Vel[1] * Vel[2] / vsquare;
    LorBoost[2][2] = (Vel[0] - 1.) * Vel[2] * Vel[2] / vsquare + 1.;
    LorBoost[2][3] = (Vel[0] - 1.) * Vel[2] * Vel[3] / vsquare;
    LorBoost[3][0] = Vel[0] * Vel[3];
    LorBoost[3][1] = (Vel[0] - 1.) * Vel[1] * Vel[3] / vsquare;
    LorBoost[3][2] = (Vel[0] - 1.) * Vel[2] * Vel[3] / vsquare;
    LorBoost[3][3] = (Vel[0] - 1.) * Vel[3] * Vel[3] / vsquare + 1.;

    if (vsquare == 0) {
      LorBoost[1][1] = 1.;
      LorBoost[1][2] = 0;
      LorBoost[1][3] = 0;
      LorBoost[2][1] = 0;
      LorBoost[2][2] = 1.;
      LorBoost[2][3] = 0;
      LorBoost[3][1] = 0;
      LorBoost[3][2] = 0;
      LorBoost[3][3] = 1.;
    }
    // Lambda_u^v Sigma_v = CMSigma_u
    CMSigma[0] = (LorBoost[0][0] * LFSigma[0] + LorBoost[0][1] * LFSigma[1] +
                  LorBoost[0][2] * LFSigma[2] + LorBoost[0][3] * LFSigma[3]);
    CMSigma[1] = (LorBoost[1][0] * LFSigma[0] + LorBoost[1][1] * LFSigma[1] +
                  LorBoost[1][2] * LFSigma[2] + LorBoost[1][3] * LFSigma[3]);
    CMSigma[2] = (LorBoost[2][0] * LFSigma[0] + LorBoost[2][1] * LFSigma[1] +
                  LorBoost[2][2] * LFSigma[2] + LorBoost[2][3] * LFSigma[3]);
    CMSigma[3] = (LorBoost[3][0] * LFSigma[0] + LorBoost[3][1] * LFSigma[1] +
                  LorBoost[3][2] * LFSigma[2] + LorBoost[3][3] * LFSigma[3]);

    // GAUSSIAN INTEGRALS <n> = int f(p)d3p
    NumLight = 0.;
    NumStrange = 0.;

    for (int l = 0; l < GPoints; l++) {
      for (int m = 0; m < GPoints; m++) {
        for (int k = 0; k < GPoints; k++) {
          GWeightProd = GWeight[l] * GWeight[m] * GWeight[k];
          pSpatialdSigma = (cut * GAbs[l]) * CMSigma[1] +
                           (cut * GAbs[m]) * CMSigma[2] +
                           (cut * GAbs[k]) * CMSigma[3];

          // For UD quarks
          E_light = computeEnergy(xmq, cut, GAbs[l], GAbs[m], GAbs[k]);
          NumLight += GWeightProd * FermiDistribution(E_light, TRead) *
                      (E_light * CMSigma[0] + pSpatialdSigma) / E_light;

          // For S quarks
          E_strange = computeEnergy(xms, cut, GAbs[l], GAbs[m], GAbs[k]);
          NumStrange += GWeightProd * FermiDistribution(E_strange, TRead) *
                        (E_strange * CMSigma[0] + pSpatialdSigma) / E_strange;
        }
      }
    }

    // U, D, UBAR, DBAR QUARKS
    // <N> = V <n>
    double NumHere = NumLight * UDDeg * cut * cut * cut / (8. * PI * PI * PI);
    std::poisson_distribution<int> poisson_ud(NumHere);

    // Generating light quarks
    GeneratedParticles = poisson_ud(
        getRandomGenerator());  // Initialize particles created in this cell

    // List of particles ( pos(x,y,z,t), mom(px,py,pz,E), species)
    for (int partic = 0; partic < GeneratedParticles; partic++) {
      // adding space to PList for output quarks
      std::vector<double> temp = {0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0.};
      Plist.push_back(temp);

      // Species - U,D,UBar,Dbar are equally likely
      double SpecRoll = ran();  // Probability of species die roll
      if (SpecRoll <= 0.25) {   // UBar
        Plist[PartCount][1] = -2;
      } else if (SpecRoll <= 0.50) {  // DBar
        Plist[PartCount][1] = -1;
      } else if (SpecRoll <= 0.75) {  // D
        Plist[PartCount][1] = 1;
      } else {  // U
        Plist[PartCount][1] = 2;
      }

      // Position
      // Located at x,y pos of area element
      Plist[PartCount][10] = CPos[0] + (ran() - 0.5) * CellDT;  // Tau
      Plist[PartCount][7] = CPos[1] + (ran() - 0.5) * CellDX;
      Plist[PartCount][8] = CPos[2] + (ran() - 0.5) * CellDY;
      Plist[PartCount][9] = CPos[3] + (ran() - 0.5) * CellDZ;

      // Momentum
      // Sample rest frame momentum given T and mass of light quark
      MCSampler(TRead, 1);  // NewP=P, NewX=Px, ...

      // USE THE SAME BOOST AS BEFORE
      // PLab^u = g^u^t Lambda_t ^v pres^w g_w _v
      // Returns P in GeV
      new_quark_energy = sqrt(xmq * xmq + NewP * NewP);
      Plist[PartCount][6] =
          (LorBoost[0][0] * new_quark_energy + LorBoost[0][1] * NewX +
           LorBoost[0][2] * NewY + LorBoost[0][3] * NewZ) *
          GEVFM;
      Plist[PartCount][3] =
          (LorBoost[1][0] * new_quark_energy + LorBoost[1][1] * NewX +
           LorBoost[1][2] * NewY + LorBoost[1][3] * NewZ) *
          GEVFM;
      Plist[PartCount][4] =
          (LorBoost[2][0] * new_quark_energy + LorBoost[2][1] * NewX +
           LorBoost[2][2] * NewY + LorBoost[2][3] * NewZ) *
          GEVFM;
      Plist[PartCount][5] =
          (LorBoost[3][0] * new_quark_energy + LorBoost[3][1] * NewX +
           LorBoost[3][2] * NewY + LorBoost[3][3] * NewZ) *
          GEVFM;

      // Additional information
      Plist[PartCount][0] = 1;   // Event ID, to match jet formatting
      Plist[PartCount][2] = 0;   // Origin, to match jet formatting
      Plist[PartCount][11] = 0;  // Status - identifies as thermal quark
      PartCount++;
    }

    // S, SBAR QUARKS
    // <N> = V <n>
    double NumOddHere =
        NumStrange * OddDeg * cut * cut * cut / (8. * PI * PI * PI);
    std::poisson_distribution<int> poisson_s(NumOddHere);

    // Generate s quarks
    int nL = GeneratedParticles;
    nL_tot += nL;

    GeneratedParticles = poisson_s(
        getRandomGenerator());  // Initialize particles created in this cell

    int nS = GeneratedParticles;
    nS_tot += nS;
    // List of particles ( pos(x,y,z,t), mom(px,py,pz,E), species)
    for (int partic = 0; partic < GeneratedParticles; partic++) {
      // adding space to PList for output quarks
      std::vector<double> temp = {0., 0., 0., 0., 0., 0.,
                                  0., 0., 0., 0., 0., 0.};
      Plist.push_back(temp);

      // Species - S,Sbar are equally likely
      double SpecRoll = ran();
      if (SpecRoll <= 0.5) {  // SBar
        Plist[PartCount][1] = -3;
      } else {  // S
        Plist[PartCount][1] = 3;
      }

      // Position
      // Located at x,y pos of area element
      Plist[PartCount][10] = CPos[0] + (ran() - 0.5) * CellDT;  // Tau
      Plist[PartCount][7] = CPos[1] + (ran() - 0.5) * CellDX;
      Plist[PartCount][8] = CPos[2] + (ran() - 0.5) * CellDY;
      Plist[PartCount][9] = CPos[3] + (ran() - 0.5) * CellDZ;

      // Momentum
      // Sample rest frame momentum given T and mass of s quark
      MCSampler(TRead, 2);  // NewP=P, NewX=Px, ...

      // USE THE SAME BOOST AS BEFORE
      // PLab^u = g^u^t Lambda_t ^v pres^w g_w _v
      // Returns P in GeV
      new_quark_energy = sqrt(xms * xms + NewP * NewP);
      Plist[PartCount][6] =
          (LorBoost[0][0] * new_quark_energy + LorBoost[0][1] * NewX +
           LorBoost[0][2] * NewY + LorBoost[0][3] * NewZ) *
          GEVFM;
      Plist[PartCount][3] =
          (LorBoost[1][0] * new_quark_energy + LorBoost[1][1] * NewX +
           LorBoost[1][2] * NewY + LorBoost[1][3] * NewZ) *
          GEVFM;
      Plist[PartCount][4] =
          (LorBoost[2][0] * new_quark_energy + LorBoost[2][1] * NewX +
           LorBoost[2][2] * NewY + LorBoost[2][3] * NewZ) *
          GEVFM;
      Plist[PartCount][5] =
          (LorBoost[3][0] * new_quark_energy + LorBoost[3][1] * NewX +
           LorBoost[3][2] * NewY + LorBoost[3][3] * NewZ) *
          GEVFM;

      // Additional information
      Plist[PartCount][0] = 1;   // Event ID, to match jet formatting
      Plist[PartCount][2] = 0;   // Origin, to match jet formatting
      Plist[PartCount][11] = 0;  // Status - identifies as thermal quark
      PartCount++;
    }
  }

  JSDEBUG << "Light particles: " << nL_tot;
  JSDEBUG << "Strange particles: " << nS_tot;

  num_ud = nL_tot;
  num_s = nS_tot;

  // Shuffling PList
  if (ShuffleList) {
    // Shuffle the Plist using the random engine
    std::shuffle(&Plist[0], &Plist[PartCount], getRandomGenerator());
  }

  // print Plist for testing
  /*std::cout << std::setprecision(5);
  std::ofstream thermalP;
  thermalP.open("thermal_partons.dat", std::ios::app);
  for(int p=0; p < Plist.size(); p++){
          thermalP << Plist[p][0] << " " << Plist[p][1] << " " << Plist[p][2]
          << " " << Plist[p][3] << " " << Plist[p][4] << " " << Plist[p][5]
          << " " << Plist[p][6] << " " << Plist[p][7] << " " << Plist[p][8]
          << " " << Plist[p][9] << " " << Plist[p][10] << " " << Plist[p][11]
          << "\n";
  }
  thermalP.close();*/
}

void ThermalPartonSampler::sample_2p1d(double eta_max) {
  // Input read from cells
  double CPos[4];     // Position of the current cell (tau/t, x, y , eta/z=0)
  double LFSigma[4];  // LabFrame hypersurface (tau/t,x,y,eta/z=0)
  double CMSigma[4];  // Center of mass hypersurface (tau/t,x,y,eta/z=0)
  double TRead;       // Temperature of cell
  double Vel[4];      // Gamma & 3-Velocity of cell (gamma, Vx, Vy, Vz) NOT FOUR
                      // VELOCITY
  double tau_pos;     // proper time from position of cell
  double eta_pos;     // eta from position of cell
  double tau_sur;     // proper time from normal vector of surface
  double eta_sur;     // eta from normal vector of surface

  // Calculated global quantities
  int PartCount;            // Total count of particles over ALL cells
  double NumLight;          // Number DENSITY of light quarks at set T
  double NumStrange;        // Number DENSITY of squarks at set T
  double UDDeg = 4. * 6.;   // Degeneracy of UD quarks
  double OddDeg = 2. * 6.;  // Degeneracy of S quarks

  // Calculated quantities in cells
  double LorBoost[4][4];   // Lorentz boost defined as used - Form is always
                           // Lambda_u^v
  int GeneratedParticles;  // Number of particles to be generated this cell

  // Define parton densities
  double cut;             // Each coordinate of P is integrated this far
  double E_light;         // Store energy of light quarks
  double E_strange;       // Store energy of strange quarks
  double GWeightProd;     // Needed for integration
  double pSpatialdSigma;  // Needed for integration
  PartCount = 0;
  NumLight = 0;  // Initialize density for light and strange quarks
  NumStrange = 0;
  GeneratedParticles = 0;
  double new_quark_energy;  // store new quark energy for boost

  // counter for total number of light and strange quarks
  int nL_tot = 0;
  int nS_tot = 0;

  // define the accuracy range for temperature values in which the CDFGenerator
  // is executed (for the case the value is not in the cache yet)
  const double accuracyRange =
      0.003 / GEVFM;  // for a usual hypersurface at const temperature there
                      // should be only one entry in the cache
  // precompute CDFs and store them in a cache
  std::unordered_map<double, std::vector<std::vector<double>>> cdfLightCache;
  std::unordered_map<double, std::vector<std::vector<double>>> cdfStrangeCache;

  // lambda function to check if a temperature value is within the accuracy
  // range of a cached temperature
  auto withinAccuracyRange = [accuracyRange](double targetTemp,
                                             double cachedTemp) {
    return std::fabs(targetTemp - cachedTemp) <= accuracyRange;
  };

  // lambda function to retrieve the closest temperature from the cache within
  // the accuracy range
  auto getClosestCachedTemp =
      [](const std::unordered_map<double, std::vector<std::vector<double>>>&
             cache,
         double targetTemp) {
        double closestTemp = -1.0;
        double minDistance = std::numeric_limits<double>::max();
        for (const auto& entry : cache) {
          double cachedTemp = entry.first;
          double distance = std::fabs(targetTemp - cachedTemp);
          if (distance < minDistance) {
            minDistance = distance;
            closestTemp = cachedTemp;
          }
        }
        return closestTemp;
      };

  // Define the lambda function to compute the energy
  auto computeEnergy = [](double mass, double cut, double GAbsL, double GAbsM,
                          double GAbsK) {
    return sqrt(mass * mass + (cut * GAbsL) * (cut * GAbsL) +
                (cut * GAbsM) * (cut * GAbsM) + (cut * GAbsK) * (cut * GAbsK));
  };

  // Define the lambda function for the Fermi distribution
  auto FermiDistribution = [](double energy, double temperature) {
    return 1. / (exp(energy / temperature) + 1.);
  };

  std::vector<double> NumLightList;
  std::vector<double> NumStrangeList;

  double d_eta = CellDZ;

  int N_slices = std::floor((eta_max / d_eta) - 0.5) + 1;
  for (int slice = 1; slice <= (2 * N_slices + 1); slice++) {
    double eta_slice = (slice - N_slices - 1) * d_eta;

    JSINFO << "Sampling thermal partons for pseudorapidity slice " << slice
           << " (eta_min = " << eta_slice - (d_eta / 2.)
           << ", eta_max = " << eta_slice + (d_eta / 2.) << ")";

    for (int iS = 0; iS < surface.size(); ++iS) {
      tau_pos = surface[iS][0];
      CPos[1] = surface[iS][1];
      CPos[2] = surface[iS][2];
      eta_pos = surface[iS][3];  // we need t,x,y,z
      // this is also tau, x, y, eta
      tau_sur = surface[iS][4];
      LFSigma[1] = surface[iS][5];
      LFSigma[2] = surface[iS][6];
      eta_sur = surface[iS][7];
      TRead = surface[iS][8] / GEVFM;
      Vel[1] = surface[iS][9];
      Vel[2] = surface[iS][10];
      Vel[3] = surface[iS][11];

      cut = 10 * TRead;

      // check if the CDFs for light quarks for this temperature are already in
      // the cache and within the accuracy range
      auto cdfLightIter = cdfLightCache.find(TRead);
      if (cdfLightIter == cdfLightCache.end()) {
        // if not found, find the closest temperature in the cache within the
        // accuracy range
        double closestTemp = getClosestCachedTemp(cdfLightCache, TRead);

        if (closestTemp >= 0. && withinAccuracyRange(TRead, closestTemp)) {
          // use the CDFs from the closest temperature in the cache
          CDFTabLight = cdfLightCache[closestTemp];
          // set the current temperature value to the closest cached temperature
          TRead = closestTemp;
        } else {
          // if no suitable entry is found, compute the CDFs and store them in
          // the cache
          CDFGenerator(TRead, xmq, 1);  // for light quarks
          cdfLightCache[TRead] = CDFTabLight;
        }
      } else {
        // if found, use the precomputed CDFs from the cache
        CDFTabLight = cdfLightIter->second;
      }

      // check if the CDFs for strange quarks for this temperature are already
      // in the cache and within the accuracy range
      auto cdfStrangeIter = cdfStrangeCache.find(TRead);
      if (cdfStrangeIter == cdfStrangeCache.end()) {
        // if not found, find the closest temperature in the cache within the
        // accuracy range
        double closestTemp = getClosestCachedTemp(cdfStrangeCache, TRead);

        if (closestTemp >= 0 && withinAccuracyRange(TRead, closestTemp)) {
          // use the CDFs from the closest temperature in the cache
          CDFTabStrange = cdfStrangeCache[closestTemp];
          // set the current temperature value to the closest cached temperature
          TRead = closestTemp;
        } else {
          // if no suitable entry is found, compute the CDFs and store them in
          // the cache
          CDFGenerator(TRead, xms, 2);  // for strange quarks
          cdfStrangeCache[TRead] = CDFTabStrange;
        }
      } else {
        // if found, use the precomputed CDFs from the cache
        CDFTabStrange = cdfStrangeIter->second;
      }

      // getting t,z from tau, eta
      CPos[0] = tau_pos * std::cosh(eta_pos);
      CPos[3] = tau_pos * std::sinh(eta_pos);
      // transform surface vector to Minkowski coordinates
      double cosh_eta_pos = std::cosh(eta_pos);
      double sinh_eta_pos = std::sinh(eta_pos);
      LFSigma[0] = cosh_eta_pos * tau_sur - (sinh_eta_pos / tau_pos) * eta_sur;
      LFSigma[3] = -sinh_eta_pos * tau_sur + (cosh_eta_pos / tau_pos) * eta_sur;

      CellDZ = CPos[0] * 2. * std::sinh(d_eta / 2.);

      double vsquare = Vel[1] * Vel[1] + Vel[2] * Vel[2] + Vel[3] * Vel[3];
      if (vsquare < 10e-16) {
        vsquare = 10e-16;
      }

      // Deduced info
      Vel[0] = 1. / sqrt(1 - vsquare);  // gamma - Vel is not four velocity

      // Lambda_u ^v
      LorBoost[0][0] = Vel[0];
      LorBoost[0][1] = Vel[0] * Vel[1];
      LorBoost[0][2] = Vel[0] * Vel[2];
      LorBoost[0][3] = Vel[0] * Vel[3];
      LorBoost[1][0] = Vel[0] * Vel[1];
      LorBoost[1][1] = (Vel[0] - 1.) * Vel[1] * Vel[1] / vsquare + 1.;
      LorBoost[1][2] = (Vel[0] - 1.) * Vel[1] * Vel[2] / vsquare;
      LorBoost[1][3] = (Vel[0] - 1.) * Vel[1] * Vel[3] / vsquare;
      LorBoost[2][0] = Vel[0] * Vel[2];
      LorBoost[2][1] = (Vel[0] - 1.) * Vel[1] * Vel[2] / vsquare;
      LorBoost[2][2] = (Vel[0] - 1.) * Vel[2] * Vel[2] / vsquare + 1.;
      LorBoost[2][3] = (Vel[0] - 1.) * Vel[2] * Vel[3] / vsquare;
      LorBoost[3][0] = Vel[0] * Vel[3];
      LorBoost[3][1] = (Vel[0] - 1.) * Vel[1] * Vel[3] / vsquare;
      LorBoost[3][2] = (Vel[0] - 1.) * Vel[2] * Vel[3] / vsquare;
      LorBoost[3][3] = (Vel[0] - 1.) * Vel[3] * Vel[3] / vsquare + 1.;

      if (vsquare == 0) {
        LorBoost[1][1] = 1.;
        LorBoost[1][2] = 0;
        LorBoost[1][3] = 0;
        LorBoost[2][1] = 0;
        LorBoost[2][2] = 1.;
        LorBoost[2][3] = 0;
        LorBoost[3][1] = 0;
        LorBoost[3][2] = 0;
        LorBoost[3][3] = 1.;
      }

      double NumHere;
      double NumOddHere;
      if (slice == 1) {
        // Lambda_u^v Sigma_v = CMSigma_u
        CMSigma[0] =
            (LorBoost[0][0] * LFSigma[0] + LorBoost[0][1] * LFSigma[1] +
             LorBoost[0][2] * LFSigma[2] + LorBoost[0][3] * LFSigma[3]);
        CMSigma[1] =
            (LorBoost[1][0] * LFSigma[0] + LorBoost[1][1] * LFSigma[1] +
             LorBoost[1][2] * LFSigma[2] + LorBoost[1][3] * LFSigma[3]);
        CMSigma[2] =
            (LorBoost[2][0] * LFSigma[0] + LorBoost[2][1] * LFSigma[1] +
             LorBoost[2][2] * LFSigma[2] + LorBoost[2][3] * LFSigma[3]);
        CMSigma[3] =
            (LorBoost[3][0] * LFSigma[0] + LorBoost[3][1] * LFSigma[1] +
             LorBoost[3][2] * LFSigma[2] + LorBoost[3][3] * LFSigma[3]);

        // GAUSSIAN INTEGRALS <n> = int f(p)d3p
        NumLight = 0.;
        NumStrange = 0.;

        for (int l = 0; l < GPoints; l++) {
          for (int m = 0; m < GPoints; m++) {
            for (int k = 0; k < GPoints; k++) {
              GWeightProd = GWeight[l] * GWeight[m] * GWeight[k];
              pSpatialdSigma = (cut * GAbs[l]) * CMSigma[1] +
                               (cut * GAbs[m]) * CMSigma[2] +
                               (cut * GAbs[k]) * CMSigma[3];

              // For UD quarks
              E_light = computeEnergy(xmq, cut, GAbs[l], GAbs[m], GAbs[k]);
              NumLight += GWeightProd * FermiDistribution(E_light, TRead) *
                          (E_light * CMSigma[0] + pSpatialdSigma) / E_light;

              // For S quarks
              E_strange = computeEnergy(xms, cut, GAbs[l], GAbs[m], GAbs[k]);
              NumStrange += GWeightProd * FermiDistribution(E_strange, TRead) *
                            (E_strange * CMSigma[0] + pSpatialdSigma) /
                            E_strange;
            }
          }
        }
        // U, D, UBAR, DBAR QUARKS
        // <N> = V <n>
        NumHere = NumLight * UDDeg * cut * cut * cut / (8. * PI * PI * PI);

        // S, SBAR QUARKS
        // <N> = V <n>
        NumOddHere =
            NumStrange * OddDeg * cut * cut * cut / (8. * PI * PI * PI);

        NumLightList.push_back(NumHere);
        NumStrangeList.push_back(NumOddHere);
      } else {
        NumHere = NumLightList[iS];
        NumOddHere = NumStrangeList[iS];
      }

      std::poisson_distribution<int> poisson_ud(NumHere);

      // Generating light quarks
      GeneratedParticles = poisson_ud(
          getRandomGenerator());  // Initialize particles created in this cell

      // List of particles ( pos(x,y,z,t), mom(px,py,pz,E), species)
      for (int partic = 0; partic < GeneratedParticles; partic++) {
        // adding space to PList for output quarks
        std::vector<double> temp = {0., 0., 0., 0., 0., 0.,
                                    0., 0., 0., 0., 0., 0.};
        Plist.push_back(temp);

        // Species - U,D,UBar,Dbar are equally likely
        double SpecRoll = ran();  // Probability of species die roll
        if (SpecRoll <= 0.25) {   // UBar
          Plist[PartCount][1] = -2;
        } else if (SpecRoll <= 0.50) {  // DBar
          Plist[PartCount][1] = -1;
        } else if (SpecRoll <= 0.75) {  // D
          Plist[PartCount][1] = 1;
        } else {  // U
          Plist[PartCount][1] = 2;
        }

        // Position
        // Located at x,y pos of area element
        Plist[PartCount][10] = CPos[0] + (ran() - 0.5) * CellDT;  // Tau
        Plist[PartCount][7] = CPos[1] + (ran() - 0.5) * CellDX;
        Plist[PartCount][8] = CPos[2] + (ran() - 0.5) * CellDY;
        Plist[PartCount][9] = CPos[3] + (ran() - 0.5) * CellDZ;

        if (std::abs(Plist[PartCount][9]) >= std::abs(Plist[PartCount][10])) {
          Plist[PartCount][10] = std::abs(Plist[PartCount][9]) + 10e-3;
        }
        double temp_t =
            std::sqrt(Plist[PartCount][10] * Plist[PartCount][10] -
                      Plist[PartCount][9] * Plist[PartCount][9]) *
            std::cosh(
                eta_slice +
                0.5 * std::log((Plist[PartCount][10] + Plist[PartCount][9]) /
                               (Plist[PartCount][10] - Plist[PartCount][9])));
        double temp_z =
            std::sqrt(Plist[PartCount][10] * Plist[PartCount][10] -
                      Plist[PartCount][9] * Plist[PartCount][9]) *
            std::sinh(
                eta_slice +
                0.5 * std::log((Plist[PartCount][10] + Plist[PartCount][9]) /
                               (Plist[PartCount][10] - Plist[PartCount][9])));
        Plist[PartCount][10] = temp_t;
        Plist[PartCount][9] = temp_z;

        // Momentum
        // Sample rest frame momentum given T and mass of light quark
        MCSampler(TRead, 1);  // NewP=P, NewX=Px, ...

        Vel[1] /= cosh(eta_slice);
        Vel[2] /= cosh(eta_slice);
        Vel[3] = tanh(eta_slice);
        double vsquare = Vel[1] * Vel[1] + Vel[2] * Vel[2] + Vel[3] * Vel[3];
        if (vsquare < 10e-16) {
          vsquare = 10e-16;
        }

        // Deduced info
        Vel[0] = 1. / sqrt(1 - vsquare);  // gamma - Vel is not four velocity

        // Lambda_u ^v
        LorBoost[0][0] = Vel[0];
        LorBoost[0][1] = Vel[0] * Vel[1];
        LorBoost[0][2] = Vel[0] * Vel[2];
        LorBoost[0][3] = Vel[0] * Vel[3];
        LorBoost[1][0] = Vel[0] * Vel[1];
        LorBoost[1][1] = (Vel[0] - 1.) * Vel[1] * Vel[1] / vsquare + 1.;
        LorBoost[1][2] = (Vel[0] - 1.) * Vel[1] * Vel[2] / vsquare;
        LorBoost[1][3] = (Vel[0] - 1.) * Vel[1] * Vel[3] / vsquare;
        LorBoost[2][0] = Vel[0] * Vel[2];
        LorBoost[2][1] = (Vel[0] - 1.) * Vel[1] * Vel[2] / vsquare;
        LorBoost[2][2] = (Vel[0] - 1.) * Vel[2] * Vel[2] / vsquare + 1.;
        LorBoost[2][3] = (Vel[0] - 1.) * Vel[2] * Vel[3] / vsquare;
        LorBoost[3][0] = Vel[0] * Vel[3];
        LorBoost[3][1] = (Vel[0] - 1.) * Vel[1] * Vel[3] / vsquare;
        LorBoost[3][2] = (Vel[0] - 1.) * Vel[2] * Vel[3] / vsquare;
        LorBoost[3][3] = (Vel[0] - 1.) * Vel[3] * Vel[3] / vsquare + 1.;

        if (vsquare == 0) {
          LorBoost[1][1] = 1.;
          LorBoost[1][2] = 0;
          LorBoost[1][3] = 0;
          LorBoost[2][1] = 0;
          LorBoost[2][2] = 1.;
          LorBoost[2][3] = 0;
          LorBoost[3][1] = 0;
          LorBoost[3][2] = 0;
          LorBoost[3][3] = 1.;
        }

        // USE THE SAME BOOST AS BEFORE
        // PLab^u = g^u^t Lambda_t ^v pres^w g_w _v
        // Returns P in GeV
        new_quark_energy = sqrt(xmq * xmq + NewP * NewP);
        Plist[PartCount][6] =
            (LorBoost[0][0] * new_quark_energy + LorBoost[0][1] * NewX +
             LorBoost[0][2] * NewY + LorBoost[0][3] * NewZ) *
            GEVFM;
        Plist[PartCount][3] =
            (LorBoost[1][0] * new_quark_energy + LorBoost[1][1] * NewX +
             LorBoost[1][2] * NewY + LorBoost[1][3] * NewZ) *
            GEVFM;
        Plist[PartCount][4] =
            (LorBoost[2][0] * new_quark_energy + LorBoost[2][1] * NewX +
             LorBoost[2][2] * NewY + LorBoost[2][3] * NewZ) *
            GEVFM;
        Plist[PartCount][5] =
            (LorBoost[3][0] * new_quark_energy + LorBoost[3][1] * NewX +
             LorBoost[3][2] * NewY + LorBoost[3][3] * NewZ) *
            GEVFM;

        // Additional information
        Plist[PartCount][0] = 1;   // Event ID, to match jet formatting
        Plist[PartCount][2] = 0;   // Origin, to match jet formatting
        Plist[PartCount][11] = 0;  // Status - identifies as thermal quark
        PartCount++;
      }

      std::poisson_distribution<int> poisson_s(NumOddHere);

      // Generate s quarks
      int nL = GeneratedParticles;
      nL_tot += nL;

      GeneratedParticles = poisson_s(
          getRandomGenerator());  // Initialize particles created in this cell

      int nS = GeneratedParticles;
      nS_tot += nS;
      // List of particles ( pos(x,y,z,t), mom(px,py,pz,E), species)
      for (int partic = 0; partic < GeneratedParticles; partic++) {
        // adding space to PList for output quarks
        std::vector<double> temp = {0., 0., 0., 0., 0., 0.,
                                    0., 0., 0., 0., 0., 0.};
        Plist.push_back(temp);

        // Species - S,Sbar are equally likely
        double SpecRoll = ran();
        if (SpecRoll <= 0.5) {  // SBar
          Plist[PartCount][1] = -3;
        } else {  // S
          Plist[PartCount][1] = 3;
        }

        // Position
        // Located at x,y pos of area element
        Plist[PartCount][10] = CPos[0] + (ran() - 0.5) * CellDT;  // Tau
        Plist[PartCount][7] = CPos[1] + (ran() - 0.5) * CellDX;
        Plist[PartCount][8] = CPos[2] + (ran() - 0.5) * CellDY;
        Plist[PartCount][9] = CPos[3] + (ran() - 0.5) * CellDZ;

        if (std::abs(Plist[PartCount][9]) >= std::abs(Plist[PartCount][10])) {
          Plist[PartCount][10] = std::abs(Plist[PartCount][9]) + 10e-3;
        }
        double temp_t =
            std::sqrt(Plist[PartCount][10] * Plist[PartCount][10] -
                      Plist[PartCount][9] * Plist[PartCount][9]) *
            std::cosh(
                eta_slice +
                0.5 * std::log((Plist[PartCount][10] + Plist[PartCount][9]) /
                               (Plist[PartCount][10] - Plist[PartCount][9])));
        double temp_z =
            std::sqrt(Plist[PartCount][10] * Plist[PartCount][10] -
                      Plist[PartCount][9] * Plist[PartCount][9]) *
            std::sinh(
                eta_slice +
                0.5 * std::log((Plist[PartCount][10] + Plist[PartCount][9]) /
                               (Plist[PartCount][10] - Plist[PartCount][9])));
        Plist[PartCount][10] = temp_t;
        Plist[PartCount][9] = temp_z;

        // Momentum
        // Sample rest frame momentum given T and mass of s quark
        MCSampler(TRead, 2);  // NewP=P, NewX=Px, ...

        Vel[1] /= cosh(eta_slice);
        Vel[2] /= cosh(eta_slice);
        Vel[3] = tanh(eta_slice);
        double vsquare = Vel[1] * Vel[1] + Vel[2] * Vel[2] + Vel[3] * Vel[3];
        if (vsquare < 10e-16) {
          vsquare = 10e-16;
        }

        // Deduced info
        Vel[0] = 1. / sqrt(1 - vsquare);  // gamma - Vel is not four velocity

        // Lambda_u ^v
        LorBoost[0][0] = Vel[0];
        LorBoost[0][1] = Vel[0] * Vel[1];
        LorBoost[0][2] = Vel[0] * Vel[2];
        LorBoost[0][3] = Vel[0] * Vel[3];
        LorBoost[1][0] = Vel[0] * Vel[1];
        LorBoost[1][1] = (Vel[0] - 1.) * Vel[1] * Vel[1] / vsquare + 1.;
        LorBoost[1][2] = (Vel[0] - 1.) * Vel[1] * Vel[2] / vsquare;
        LorBoost[1][3] = (Vel[0] - 1.) * Vel[1] * Vel[3] / vsquare;
        LorBoost[2][0] = Vel[0] * Vel[2];
        LorBoost[2][1] = (Vel[0] - 1.) * Vel[1] * Vel[2] / vsquare;
        LorBoost[2][2] = (Vel[0] - 1.) * Vel[2] * Vel[2] / vsquare + 1.;
        LorBoost[2][3] = (Vel[0] - 1.) * Vel[2] * Vel[3] / vsquare;
        LorBoost[3][0] = Vel[0] * Vel[3];
        LorBoost[3][1] = (Vel[0] - 1.) * Vel[1] * Vel[3] / vsquare;
        LorBoost[3][2] = (Vel[0] - 1.) * Vel[2] * Vel[3] / vsquare;
        LorBoost[3][3] = (Vel[0] - 1.) * Vel[3] * Vel[3] / vsquare + 1.;

        if (vsquare == 0) {
          LorBoost[1][1] = 1.;
          LorBoost[1][2] = 0;
          LorBoost[1][3] = 0;
          LorBoost[2][1] = 0;
          LorBoost[2][2] = 1.;
          LorBoost[2][3] = 0;
          LorBoost[3][1] = 0;
          LorBoost[3][2] = 0;
          LorBoost[3][3] = 1.;
        }

        // USE THE SAME BOOST AS BEFORE
        // PLab^u = g^u^t Lambda_t ^v pres^w g_w _v
        // Returns P in GeV
        new_quark_energy = sqrt(xms * xms + NewP * NewP);
        Plist[PartCount][6] =
            (LorBoost[0][0] * new_quark_energy + LorBoost[0][1] * NewX +
             LorBoost[0][2] * NewY + LorBoost[0][3] * NewZ) *
            GEVFM;
        Plist[PartCount][3] =
            (LorBoost[1][0] * new_quark_energy + LorBoost[1][1] * NewX +
             LorBoost[1][2] * NewY + LorBoost[1][3] * NewZ) *
            GEVFM;
        Plist[PartCount][4] =
            (LorBoost[2][0] * new_quark_energy + LorBoost[2][1] * NewX +
             LorBoost[2][2] * NewY + LorBoost[2][3] * NewZ) *
            GEVFM;
        Plist[PartCount][5] =
            (LorBoost[3][0] * new_quark_energy + LorBoost[3][1] * NewX +
             LorBoost[3][2] * NewY + LorBoost[3][3] * NewZ) *
            GEVFM;

        // Additional information
        Plist[PartCount][0] = 1;   // Event ID, to match jet formatting
        Plist[PartCount][2] = 0;   // Origin, to match jet formatting
        Plist[PartCount][11] = 0;  // Status - identifies as thermal quark
        PartCount++;
      }
    }
  }

  JSDEBUG << "Light particles: " << nL_tot;
  JSDEBUG << "Strange particles: " << nS_tot;

  num_ud = nL_tot;
  num_s = nS_tot;

  // Shuffling PList
  if (ShuffleList) {
    // Shuffle the Plist using the random engine
    std::shuffle(&Plist[0], &Plist[PartCount], getRandomGenerator());
  }

  // print Plist for testing
  /*std::cout << std::setprecision(5);
  std::ofstream thermalP;
  thermalP.open("thermal_partons.dat", std::ios::app);
  for(int p=0; p < Plist.size(); p++){
          thermalP << Plist[p][0] << " " << Plist[p][1] << " " << Plist[p][2]
          << " " << Plist[p][3] << " " << Plist[p][4] << " " << Plist[p][5]
          << " " << Plist[p][6] << " " << Plist[p][7] << " " << Plist[p][8]
          << " " << Plist[p][9] << " " << Plist[p][10] << " " << Plist[p][11]
          << "\n";
  }
  thermalP.close();*/
}
