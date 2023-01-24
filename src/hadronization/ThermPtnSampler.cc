
#include "ThermPtnSampler.h"
#include "JetScapeXML.h"
#include "JetScapeLogger.h"
#include "JetScapeConstants.h"
#include <vector>
#include <random>

using namespace Jetscape;

thermptnsampler::thermptnsampler(){
	
	/* Static Parameters, Do not change */
	PIE = 3.141592653589793;
	GEVFM = 0.1975;
	muPi0 = 0.;
	
	/* Gaussian Weights for Integration */
	GWeight[0]  = 0.0621766166553473; GWeight[1]  = 0.0619360674206832; GWeight[2]  = 0.0614558995903167; GWeight[3]  = 0.0607379708417702; GWeight[4]  = 0.0597850587042655;
	GWeight[5]  = 0.0586008498132224; GWeight[6]  = 0.0571899256477284; GWeight[7]  = 0.0555577448062125; GWeight[8]  = 0.0537106218889962; GWeight[9]  = 0.0516557030695811;
	GWeight[10] = 0.0494009384494663; GWeight[11] = 0.0469550513039484; GWeight[12] = 0.0443275043388033; GWeight[13] = 0.0415284630901477; GWeight[14] = 0.0385687566125877;
	GWeight[15] = 0.0354598356151462; GWeight[16] = 0.0322137282235780; GWeight[17] = 0.0288429935805352; GWeight[18] = 0.0253606735700124; GWeight[19] = 0.0217802431701248;
	GWeight[20] = 0.0181155607134894; GWeight[21] = 0.0143808227614856; GWeight[22] = 0.0105905483836510; GWeight[23] = 0.0067597991957454; GWeight[24] = 0.0029086225531551;
	GWeight[25] = 0.0621766166553473; GWeight[26] = 0.0619360674206832; GWeight[27] = 0.0614558995903167; GWeight[28] = 0.0607379708417702; GWeight[29] = 0.0597850587042655;
	GWeight[30] = 0.0586008498132224; GWeight[31] = 0.0571899256477284; GWeight[32] = 0.0555577448062125; GWeight[33] = 0.0537106218889962; GWeight[34] = 0.0516557030695811;
	GWeight[35] = 0.0494009384494663; GWeight[36] = 0.0469550513039484; GWeight[37] = 0.0443275043388033; GWeight[38] = 0.0415284630901477; GWeight[39] = 0.0385687566125877;
	GWeight[40] = 0.0354598356151462; GWeight[41] = 0.0322137282235780; GWeight[42] = 0.0288429935805352; GWeight[43] = 0.0253606735700124; GWeight[44] = 0.0217802431701248;
	GWeight[45] = 0.0181155607134894; GWeight[46] = 0.0143808227614856; GWeight[47] = 0.0105905483836510; GWeight[48] = 0.0067597991957454; GWeight[49] = 0.0029086225531551;
	
	/* Gaussian Abscissas for Integration */
	GAbs[0]  =  0.0310983383271889; GAbs[1]  =  0.0931747015600861; GAbs[2]  =  0.1548905899981459; GAbs[3]  =  0.2160072368760418; GAbs[4]  =  0.2762881937795320;
	GAbs[5]  =  0.3355002454194373; GAbs[6]  =  0.3934143118975651; GAbs[7]  =  0.4498063349740388; GAbs[8]  =  0.5044581449074642; GAbs[9]  =  0.5571583045146501;
	GAbs[10] =  0.6077029271849502; GAbs[11] =  0.6558964656854394; GAbs[12] =  0.7015524687068222; GAbs[13] =  0.7444943022260685; GAbs[14] =  0.7845558329003993;
	GAbs[15] =  0.8215820708593360; GAbs[16] =  0.8554297694299461; GAbs[17] =  0.8859679795236131; GAbs[18] =  0.9130785566557919; GAbs[19] =  0.9366566189448780;
	GAbs[20] =  0.9566109552428079; GAbs[21] =  0.9728643851066920; GAbs[22] =  0.9853540840480058; GAbs[23] =  0.9853540840480058; GAbs[24] =  0.9988664044200710;
	GAbs[25] = -0.0310983383271889; GAbs[26] = -0.0931747015600861; GAbs[27] = -0.1548905899981459; GAbs[28] = -0.2160072368760418; GAbs[29] = -0.2762881937795320;
	GAbs[30] = -0.3355002454194373; GAbs[31] = -0.3934143118975651; GAbs[32] = -0.4498063349740388; GAbs[33] = -0.5044581449074642; GAbs[34] = -0.5571583045146501;
	GAbs[35] = -0.6077029271849502; GAbs[36] = -0.6558964656854394; GAbs[37] = -0.7015524687068222; GAbs[38] = -0.7444943022260685; GAbs[39] = -0.7845558329003993;
	GAbs[40] = -0.8215820708593360; GAbs[41] = -0.8554297694299461; GAbs[42] = -0.8859679795236131; GAbs[43] = -0.9130785566557919; GAbs[44] = -0.9366566189448780;
	GAbs[45] = -0.9566109552428079; GAbs[46] = -0.9728643851066920; GAbs[47] = -0.9853540840480058; GAbs[48] = -0.9853540840480058; GAbs[49] = -0.9988664044200710;

	/* Adjustable Parameters */
	ML = 0.3/GEVFM;
	MS = 0.45/GEVFM;
	T = 0.165/GEVFM;  // 165 MeV into fm^-1
	//NUMSTEP = 1025;  // 2^10+1, for steps of CDF Table, changes coarseness of momentum sampling
	NUMSTEP = 1048577;  // 2^20+1, for steps of CDF Table, changes coarseness of momentum sampling
	
	/* Adjustable params for 3+1d */ //maybe DEta here instead of above? since needed for 2+1d
	//DEta = 0.5;
	CellDX = 0.1; CellDY = 0.1;	CellDZ = 0.1; CellDT = 0.1; //may need to make these read in from hypersurface instead

	/* Flags */
	SetNum = false; // Set 'true' to set number of particles by hand- !!!Statistics use above temperature!!!
	SetNumLight = 1000000; //If SetNum == true, this many UD quarks are generated
	SetNumStrange = 0 ; //If SetNum == true, this many S quarks are generated
	ShuffleList = true; // Should list of particles be shuffled at the end?
	
	/* Brick Info */
	L = 4.0; // Thickness from box edge
	W = 4.0; // x/y width of box
	offset = 4.0; // Edge of box location

	//Vx = 0.8; // Uniform flow in x-dir
	//Vy = 0.4; // Uniform flow in y-dir
	//Vz = 0.1; // Uniform flow in z-dir

	Vx = 0.; // 'Uniform' no flow in x-dir
	Vy = 0.; // 'Uniform' no flow in y-dir
	Vz = 0.; // 'Uniform' no flow in z-dir
	
	surface.clear();
	
	num_ud = 0; num_s = 0;
	
	for(int iCDF=0; iCDF<NUMSTEP; ++iCDF){std::vector<double> temp = {0., 0.}; CDFTabLight.push_back(temp); CDFTabStrange.push_back(temp);}
	
	// random seed
	if(ran_seed != 0){eng.seed(ran_seed);}
	//else{eng.seed(std::random_device{}());}
	else{ //seeding the mt19937_64 object 'eng' PROPERLY!
		std::random_device rd;
		std::array<int,std::mt19937_64::state_size> seedarray; std::generate_n(seedarray.data(), seedarray.size(), std::ref(rd));
		std::seed_seq seeds(std::begin(seedarray), std::end(seedarray)); eng.seed(seeds);
	}
	
}

double thermptnsampler::ran() {std::uniform_real_distribution<double> uniran(0.0,1.0); return uniran(eng);}

double thermptnsampler::FermiPDF (double P, double M, double T, double mu){return 1./( exp( (sqrt(M*M + P*P) - mu)/T ) + 1. );}
//double BosePDF (double P, double M, double T, double mu){return 1./( exp( (sqrt(M*M + P*P) - mu)/T ) - 1. );}
//double BoltzPDF (double P, double M, double T, double mu){return exp( -(sqrt(M*M + P*P) - mu)/T );} //Returns ProbDensity from Boltzman (Normalization handled separate)

bool thermptnsampler::SplitSort (double goal, int floor, int ceiling, int quark){
	if(quark==1){
		int TargetPoint = ( (floor + ceiling)/2 );
		double TargetVal = CDFTabLight[TargetPoint][1];
		
		if(goal > TargetVal){return (true);}
		else{return (false);}
	}
	else if(quark==2){
		int TargetPoint = ( (floor + ceiling)/2 );
		double TargetVal = CDFTabStrange[TargetPoint][1];
		
		if(goal > TargetVal){return (true);}
		else{return (false);}
	}
	return false;
}

void thermptnsampler::CDFGen(double Temp, double M, int quark){
	double PStep;
	double PMax = 10.*Temp;  //CutOff for Integration
	PStep = PMax/(NUMSTEP-1); //Stepsize in P
		
	if(quark==1){
		/*** Initialize tabulated resuls for CDF ***/
		CDFTabLight[0][0] = 0; // For zero momentum or less...        
		CDFTabLight[0][1] = 0;// There is zero chance
		
		/*** Tabulate CDF(x) = int(0->x) PDF ***/
		for(int i =1; i<NUMSTEP; i++){
			CDFTabLight[i][0] = CDFTabLight[i-1][0] + PStep; // Calculate Momentum of next step
			double Fermi0 = FermiPDF(CDFTabLight[i-1][0], (M), Temp, muPi0); // PDF
			double Fermi1 = FermiPDF(CDFTabLight[i][0], (M), Temp, muPi0);
			CDFTabLight[i][1] = CDFTabLight[i-1][1] + (PStep/2)*( Fermi0*CDFTabLight[i-1][0]*CDFTabLight[i-1][0] + Fermi1*CDFTabLight[i][0]*CDFTabLight[i][0] );
		}
		
		for(int i =0; i<NUMSTEP; i++){CDFTabLight[i][1] = CDFTabLight[i][1]/CDFTabLight[NUMSTEP-1][1];}
	}
	else if(quark==2){
		/*** Initialize tabulated resuls for CDF ***/
		CDFTabStrange[0][0] = 0; // For zero momentum or less...        
		CDFTabStrange[0][1] = 0;// There is zero chance

		/*** Tabulate CDF(x) = int(0->x) PDF ***/
		for(int i =1; i<NUMSTEP; i++){
			CDFTabStrange[i][0] = CDFTabStrange[i-1][0] + PStep; // Calculate Momentum of next step
			double Fermi0 = FermiPDF(CDFTabStrange[i-1][0], (M), Temp, muPi0); // PDF
			double Fermi1 = FermiPDF(CDFTabStrange[i][0], (M), Temp, muPi0);
			CDFTabStrange[i][1] = CDFTabStrange[i-1][1] + (PStep/2)*( Fermi0*CDFTabStrange[i-1][0]*CDFTabStrange[i-1][0] + Fermi1*CDFTabStrange[i][0]*CDFTabStrange[i][0] );
		} 
		
		for(int i =0; i<NUMSTEP; i++){CDFTabStrange[i][1] = CDFTabStrange[i][1]/CDFTabStrange[NUMSTEP-1][1];}
	}
}

void thermptnsampler::MCSampler(double Temp, double M, int quark){ // Input in fm^-1
	double PMag, CosT, Phi, PStep;
	double PMax = 10.*Temp;  //CutOff for Integration
	PStep = PMax/(NUMSTEP-1); //Stepsize in P
	
	if(quark==1){
		int flr = 0;
		int clg = NUMSTEP - 1;
		double PRoll = ran();

		for(int i = 0; i<20;i++){
			if(SplitSort(PRoll, flr, clg, quark)){flr = ( (flr + clg)/2 );}
			else{clg = ( (flr + clg)/2 );}
		}
		
		PMag = PStep*(PRoll - CDFTabLight[flr][1])/(CDFTabLight[clg][1] - CDFTabLight[flr][1]) + CDFTabLight[flr][0];      
		CosT =( ran() - 0.5)*2.0;
		Phi = ran()*2*PIE;

		NewX = PMag*sqrt(1-CosT*CosT)*cos(Phi);
		NewY = PMag*sqrt(1-CosT*CosT)*sin(Phi);
		NewZ = PMag*CosT;
		NewP = PMag;
	}
	else if(quark==2){
		int flr = 0;
		int clg = NUMSTEP - 1;
		double PRoll = ran();
		
		for(int i = 0; i<10;i++){
			if(SplitSort(PRoll, flr, clg, quark)){flr = ( (flr + clg)/2 );}
			else{clg = ( (flr + clg)/2 );}
		}
		
		PMag = PStep*(PRoll - CDFTabStrange[flr][1])/(CDFTabStrange[clg][1] - CDFTabStrange[flr][1]) + CDFTabStrange[flr][0];      
		CosT =( ran() - 0.5)*2.0;
		Phi = ran()*2*PIE;
		
		NewX = PMag*sqrt(1-CosT*CosT)*cos(Phi);
		NewY = PMag*sqrt(1-CosT*CosT)*sin(Phi);
		NewZ = PMag*CosT;
		NewP = PMag;
	}
}

//void hypersurface(std::vector<SurfaceCellInfo> surf_in){
//	surface = surf_in;
//}

void thermptnsampler::samplebrick(){
	
	//preliminary parameter checks
	if(L < 0.){L = -L; JSWARN << "Negative brick length - setting to positive " << L << " fm.";}
	if(W < 0.){W = -W; JSWARN << "Negative brick width - setting to positive " << W << " fm.";}
	if(offset < 0.){ JSWARN << "Negative offset - Are you sure that you want to do this?";} //I guess you could, but why?
	
	/* Input Read from Cells */
	//double CPos[4];		// Position of the current cell (tau/t, x, y , eta/z=0)
	double LFSigma[4];		// LabFrame hypersurface (tau/t,x,y,eta/z=0)
	double CMSigma[4];		// Center of Mass hypersurface (tau/t,x,y,eta/z=0)
	//double TRead;			// Temperature of cell (from file) Not used.
	double Vel[4];			// Gamma & 3-Velocity of cell (gamma, Vx, Vy, Vz) NOT FOUR VELOCITY
	
	/* Calculated Global Quantities */
	int PartCount;			// Total Count of Particles over ALL cells
	double NumLight;		// Number DENSITY of light quarks at set T
	double NumStrange;		// Number DENSITY of squarks at set T
	double UDDeg = 4.*6.;	// Degeneracy of UD quarks
	double OddDeg = 2.*6.;	// Degeneracy of Squarks
	
	/* Calculated Quantities in Cells */
	double LorBoost[4][4];	// Lorentz boost defined as used- Form is always Lambda_u^v
	int GeneratedParticles;	// Number of particles to be generated this cell
	
	/* DebugTools */
	int RollCounter = 0;		// How many rolls?
	
	/*** End Definition of Static Variables ***/
	
	/*For parity in coding style for Accept/Reject code*/
	LFSigma[0] = 1.;
	LFSigma[1] = 0.;
	LFSigma[2] = 0.;
	LFSigma[3] = 0.;

	Vel[1] = Vx; Vel[2] = Vy; Vel[3] = Vz;

	if(Vel[1]*Vel[1] + Vel[2]*Vel[2] + Vel[3]*Vel[3] > 1. ){
		JSWARN << "v^2 = " << Vel[1]*Vel[1] + Vel[2]*Vel[2] + Vel[3]*Vel[3];
		JSWARN << "Unphysical velocity (brick flow)! Set to \"No Flow\" case";
		
		Vel[1] = 0.; Vel[2] = 0.; Vel[3] = 0.;
	}
	
	Vel[0] = 1. / std::sqrt(1. - Vel[1]*Vel[1] - Vel[2]*Vel[2] - Vel[3]*Vel[3] );   // gamma - Vel is not four velocity
	
	double vsquare = ( Vel[1]*Vel[1] + Vel[2]*Vel[2] + Vel[3]*Vel[3]  );
	
	/*** Lambda_u ^v from (lab frame to rest frame with flow velocity) ***/
	
	if(vsquare == 0){
		LorBoost[0][0]= Vel[0];
		LorBoost[0][1]= Vel[0]*Vel[1];
		LorBoost[0][2]= Vel[0]*Vel[2];
		LorBoost[0][3]= Vel[0]*Vel[3];
		LorBoost[1][0]= Vel[0]*Vel[1];
		LorBoost[1][1]= 1.;  
		LorBoost[1][2]= 0.;
		LorBoost[1][3]= 0.;
		LorBoost[2][0]= Vel[0]*Vel[2];
		LorBoost[2][1]= 0.;
		LorBoost[2][2]= 1.;
		LorBoost[2][3]= 0.;
		LorBoost[3][0]= Vel[0]*Vel[3];
		LorBoost[3][1]= 0.;
		LorBoost[3][2]= 0.;
		LorBoost[3][3]= 1.;
	}
	else{
		LorBoost[0][0]= Vel[0];
		LorBoost[0][1]= Vel[0]*Vel[1];
		LorBoost[0][2]= Vel[0]*Vel[2];
		LorBoost[0][3]= Vel[0]*Vel[3];
		LorBoost[1][0]= Vel[0]*Vel[1];
		LorBoost[1][1]= (Vel[0] - 1.)*Vel[1]*Vel[1]/vsquare + 1.;  
		LorBoost[1][2]= (Vel[0] - 1.)*Vel[1]*Vel[2]/vsquare;
		LorBoost[1][3]= (Vel[0] - 1.)*Vel[1]*Vel[3]/vsquare;
		LorBoost[2][0]= Vel[0]*Vel[2];
		LorBoost[2][1]= (Vel[0] - 1.)*Vel[1]*Vel[2]/vsquare;
		LorBoost[2][2]= (Vel[0] - 1.)*Vel[2]*Vel[2]/vsquare + 1.;
		LorBoost[2][3]= (Vel[0] - 1.)*Vel[2]*Vel[3]/vsquare;
		LorBoost[3][0]= Vel[0]*Vel[3];
		LorBoost[3][1]= (Vel[0] - 1.)*Vel[1]*Vel[3]/vsquare;
		LorBoost[3][2]= (Vel[0] - 1.)*Vel[2]*Vel[3]/vsquare;
		LorBoost[3][3]= (Vel[0] - 1.)*Vel[3]*Vel[3]/vsquare + 1.;
	}
	
	/*** Code parity with hypersurface case ***/
	/*** Lambda_u^v Sigma_v = CMSigma_u ***/
	CMSigma[0] = (LorBoost[0][0]*LFSigma[0] + LorBoost[0][1]*LFSigma[1] + LorBoost[0][2]*LFSigma[2] + LorBoost[0][3]*LFSigma[3]);
	CMSigma[1] = (LorBoost[1][0]*LFSigma[0] + LorBoost[1][1]*LFSigma[1] + LorBoost[1][2]*LFSigma[2] + LorBoost[1][3]*LFSigma[3]);
	CMSigma[2] = (LorBoost[2][0]*LFSigma[0] + LorBoost[2][1]*LFSigma[1] + LorBoost[2][2]*LFSigma[2] + LorBoost[2][3]*LFSigma[3]);
	CMSigma[3] = (LorBoost[3][0]*LFSigma[0] + LorBoost[3][1]*LFSigma[1] + LorBoost[3][2]*LFSigma[2] + LorBoost[3][3]*LFSigma[3]);
	
	/* Define Parton Densities */
	
	double cut = 10.*T; // Each coordinate of P is integrated this far
	PartCount = 0;
	NumLight = 0;       // Initialize density for Light and strange quarks
	NumStrange = 0;
	
	/*** GAUSSIAN INTEGRALS <n> = int f(p)d3p ***/
	
	for (int l=0; l<GPoints; l++){
		for (int m=0; m<GPoints; m++){
			for(int k=0; k<GPoints; k++){
				/*** For UD quarks ***/
				NumLight = NumLight + (GWeight[l]*GWeight[m]*GWeight[k]*1./(exp(sqrt(ML*ML + (cut*GAbs[l])*(cut*GAbs[l]) + (cut*GAbs[m])*(cut*GAbs[m]) + (cut*GAbs[k])*(cut*GAbs[k]))/T) + 1.))*
				  (sqrt(ML*ML + (cut*GAbs[l])*(cut*GAbs[l]) + (cut*GAbs[m])*(cut*GAbs[m]) + (cut*GAbs[k])*(cut*GAbs[k]))*CMSigma[0] + (cut*GAbs[l])*CMSigma[1] + (cut*GAbs[m])*CMSigma[2] + (cut*GAbs[k])*CMSigma[3])/
				  (sqrt(ML*ML + (cut*GAbs[l])*(cut*GAbs[l]) + (cut*GAbs[m])*(cut*GAbs[m]) + (cut*GAbs[k])*(cut*GAbs[k])));
					   
				   /*** For Squarks ***/
				NumStrange = NumStrange + (GWeight[l]*GWeight[m]*GWeight[k] * 1./(exp(((sqrt(MS*MS + (cut*GAbs[l])*(cut*GAbs[l]) + (cut*GAbs[m])*(cut*GAbs[m]) + (cut*GAbs[k])*(cut*GAbs[k]))))/T) + 1.))*
				  (sqrt(MS*MS + (cut*GAbs[l])*(cut*GAbs[l]) + (cut*GAbs[m])*(cut*GAbs[m]) + (cut*GAbs[k])*(cut*GAbs[k]))*CMSigma[0] + (cut*GAbs[l])*CMSigma[1] + (cut*GAbs[m])*CMSigma[2] + (cut*GAbs[k])*CMSigma[3])/
				  (sqrt(MS*MS + (cut*GAbs[l])*(cut*GAbs[l]) + (cut*GAbs[m])*(cut*GAbs[m]) + (cut*GAbs[k])*(cut*GAbs[k])));
			}
		}
	}
	
	/*** Normalization factors: degeneracy, gaussian integration, and 1/(2Pi)^3 ***/
	NumLight = NumLight*UDDeg*cut*cut*cut/(8.*PIE*PIE*PIE);
	NumStrange = NumStrange*OddDeg*cut*cut*cut/(8.*PIE*PIE*PIE);
	
	/** Define Rest-frame Cumulative functions**/
	CDFGen(T, ML, 1); // for light quarks
	CDFGen(T, MS, 2); // for squarks
	//cout << CDFTabLight[512][0] << "  " << CDFTabLight[512][1] << "  " << CDFTabLight[1024][0] << "  " << CDFTabLight[1024][1] << endl;
	
	/*** U, D, UBAR, DBAR QUARKS ***/
	/*** <N> = V <n> ***/
	double NumHere=NumLight*L*W*W;

	/*** S, SBAR QUARKS ***/
	/*** <N> = V <n> ***/
	double NumOddHere=NumStrange*L*W*W;
	
	/*~~~~~~~~~~~~~If Overwriting~~~~~~~~~~*/
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	if(SetNum){NumHere = SetNumLight; NumOddHere = SetNumStrange;}
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	
	/* Generating Light Quarks */
	GeneratedParticles = (int) std::round(NumHere); //UD quarks to be Generated in the brick; Currently, non fluctuating
	
	//adding space to PList for output quarks
	for(int iPl=0; iPl<int(std::round(NumHere)+std::round(NumOddHere)); ++iPl){std::vector<double> temp = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}; Plist.push_back(temp);}
	
	//List of Particles ( pos(x,y,z,t), mom(px,py,pz,E), species)
	for(int partic=0; partic<GeneratedParticles; partic++){
		/* Species - U,D,UBar,Dbar are equally likely */
		double SpecRoll = ran(); //Probability of species die roll
		SpecRoll = ran();
		
		if(SpecRoll <=0.25){Plist[PartCount][1] = -2;} // DBar
		else{
			if(SpecRoll <= 0.50){Plist[PartCount][1] = -1;} // UBar
			else{
				if(SpecRoll <=0.75){Plist[PartCount][1] = 1;} // U
				else{Plist[PartCount][1] = 2;} // D
			}
		}
		
		/* Position */
		/*** Located at x,y pos of area element ***/
		double XRoll = ran()-0.5; // center at x=0
		double YRoll = ran()-0.5; // center at y=0
		//double ZRoll = ran(); // particle within L of offset
		double ZRoll = ran()-0.5; // particle within L of offset
		
		Plist[PartCount][7] = XRoll*L;
		Plist[PartCount][8] = YRoll*W;
		//Plist[PartCount][9] = offset-ZRoll*L;
		Plist[PartCount][9] = ZRoll*W;
		
		/*** Time ***/
		//Plist[PartCount][10] = offset; // Tau = offset: assume jet at light speed
		Plist[PartCount][10] = L/2.; // Tau = L/2.: assume jet at light speed

		
		/* Momentum */
		/*** Sample rest frame momentum given T and Mass of light quark ***/
		MCSampler(T, ML, 1); // NewP= P, NewX=Px, ...
		
		bool IsAccept = false;
		while(!IsAccept){
			RollCounter++;
			double Chance = ran()*(CMSigma[0] + sqrt(CMSigma[1]*CMSigma[1] + CMSigma[2]*CMSigma[2] + CMSigma[3]*CMSigma[3]));
			double Post = ( sqrt(ML*ML + NewP*NewP)*CMSigma[0] + NewX*CMSigma[1] + NewY*CMSigma[2] + NewZ*CMSigma[3])/sqrt(ML*ML + NewP*NewP);
			
			if(Chance>Post){MCSampler(T, ML, 1);} //REJECTED  // NewP= P, NewX=Px, ...
			else{IsAccept = true;}
		}
		
		/*** PLab^u = g^u^t Lambda_t ^v pres^w g_w _v  ***/
		/*** This is boost as if particle sampled in rest frame ***/
		Plist[PartCount][6]= (LorBoost[0][0]*sqrt(ML*ML + NewP*NewP) + LorBoost[0][1]*NewX + LorBoost[0][2]*NewY + LorBoost[0][3]*NewZ)*GEVFM;
		Plist[PartCount][3]= (LorBoost[1][0]*sqrt(ML*ML + NewP*NewP) + LorBoost[1][1]*NewX + LorBoost[1][2]*NewY + LorBoost[1][3]*NewZ)*GEVFM;
		Plist[PartCount][4]= (LorBoost[2][0]*sqrt(ML*ML + NewP*NewP) + LorBoost[2][1]*NewX + LorBoost[2][2]*NewY + LorBoost[2][3]*NewZ)*GEVFM;
		Plist[PartCount][5]= (LorBoost[3][0]*sqrt(ML*ML + NewP*NewP) + LorBoost[3][1]*NewX + LorBoost[3][2]*NewY + LorBoost[3][3]*NewZ)*GEVFM;

		/*** Above Returns P in GeV ***/
		/** Additional information **/
		Plist[PartCount][0]  = 1; //Event ID, to match jet formatting
		Plist[PartCount][2]  = 0; //Origin, to match jet formatting
		Plist[PartCount][11] = 0; //Status- Identifies as thermal quark
		PartCount++;
	}
	
	/* Generate Squarks */
	
	GeneratedParticles = (int) std::round(NumOddHere); //Squarks to be Generated in the brick; Currently, non fluctuating
	
	//List of Particles ( pos(x,y,z,t), mom(px,py,pz,E), species)
	for(int partic=0; partic<GeneratedParticles; partic++){
		/* Species - S,Sbar are equally likely */
		double SpecRoll = ran();
		if(SpecRoll <=0.50000){Plist[PartCount][1] = -3;} // SBar
		else{Plist[PartCount][1] = 3;} //S
		
		/* Position */
		/*** Located at x,y pos of area element ***/
		double XRoll = ran()-0.5; // center at x=0
		double YRoll = ran()-0.5; // center at y=0
		//double ZRoll = ran(); // particle within L of offset
		double ZRoll = ran()-0.5; // particle within L of offset
		
		Plist[PartCount][7] = XRoll*L;
		Plist[PartCount][8] = YRoll*W;
		//Plist[PartCount][9] = offset-ZRoll*L;
		Plist[PartCount][9] = ZRoll*W;
		
		/*** Time ***/
		//Plist[PartCount][10] = L; // Tau = L: assume jet at light speed
		Plist[PartCount][10] = L/2.; // Tau = L/2.: assume jet at light speed

		
		/* Momentum */
		/*** Sample rest frame momentum given T and Mass of squark ***/
		MCSampler(T, MS, 2); // NewP= P, NewX=Px, ...
		
		bool IsAccept = false;
		while(!IsAccept){
			RollCounter++;
			double Chance = ran()*(CMSigma[0] + sqrt(CMSigma[1]*CMSigma[1] + CMSigma[2]*CMSigma[2] + CMSigma[3]*CMSigma[3]));
			double Post = ( sqrt(MS*MS + NewP*NewP)*CMSigma[0] + NewX*CMSigma[1] + NewY*CMSigma[2] + NewZ*CMSigma[3])/sqrt(MS*MS + NewP*NewP);

			if(Chance>Post){MCSampler(T, MS, 1);} //REJECTED  // NewP= P, NewX=Px, ...
			else{IsAccept = true;}
		}
		
		/*** PLab^u = g^u^t Lambda_t ^v Pres^w g_w _v  = Lambda ^u _w Pres_w (with velocity -v)***/
		/*** (Lambda _u ^t with velocity v) == (Lambda ^u _t with velocity -v) ***/
		/*** This is boost as if particle sampled in rest frame ***/
		Plist[PartCount][6]= (LorBoost[0][0]*sqrt(MS*MS + NewP*NewP) + LorBoost[0][1]*NewX + LorBoost[0][2]*NewY + LorBoost[0][3]*NewZ)*GEVFM;
		Plist[PartCount][3]= (LorBoost[1][0]*sqrt(MS*MS + NewP*NewP) + LorBoost[1][1]*NewX + LorBoost[1][2]*NewY + LorBoost[1][3]*NewZ)*GEVFM;
		Plist[PartCount][4]= (LorBoost[2][0]*sqrt(MS*MS + NewP*NewP) + LorBoost[2][1]*NewX + LorBoost[2][2]*NewY + LorBoost[2][3]*NewZ)*GEVFM;
		Plist[PartCount][5]= (LorBoost[3][0]*sqrt(MS*MS + NewP*NewP) + LorBoost[3][1]*NewX + LorBoost[3][2]*NewY + LorBoost[3][3]*NewZ)*GEVFM;
		/*** Above Returns P in GeV ***/
		/** Additional information **/
		Plist[PartCount][0]  = 1; //Event ID, to match jet formatting
		Plist[PartCount][2]  = 0; //Origin, to match jet formatting
		Plist[PartCount][11] = 0; //Status- Identifies as thermal quark
		PartCount++;
	}
	
	JSDEBUG << "Light Particles: " << (int) std::round(NumHere);
	JSDEBUG << "Strange Particles: " << (int) std::round(NumOddHere);
	
	num_ud = (int) std::round(NumHere); num_s = (int) std::round(NumOddHere);
	
	//Shuffling PList
	if(ShuffleList){
	for(int i=0;i<PartCount-1;++i){
		int ranelement = i+floor((PartCount-i)*ran());
		double temp[12]; 
		for(int x=0; x<12;x++){temp[x]= Plist[i][x];}
		for(int x=0; x<12;x++){Plist[i][x] = Plist[ranelement][x];}
		for(int x=0; x<12;x++){Plist[ranelement][x] = temp[x];}
	}}

	JSDEBUG << "For " << NumHere << " Light Quarks and " << NumOddHere << " strange quarks, we needed:";
	JSDEBUG << RollCounter << " Total Rolls";
}

void thermptnsampler::sample_3p1d(){
	
	/* Input Read from Cells */
	double CPos[4];		// Position of the current cell (tau/t, x, y , eta/z=0)
	double LFSigma[4];		// LabFrame hypersurface (tau/t,x,y,eta/z=0)
	double CMSigma[4];		// Center of Mass hypersurface (tau/t,x,y,eta/z=0)
	double TRead;			// Temperature of cell (from file) Not used.
	double Vel[4];			// Gamma & 3-Velocity of cell (gamma, Vx, Vy, Vz) NOT FOUR VELOCITY
	
	/* Calculated Global Quantities */
	int PartCount;			// Total Count of Particles over ALL cells
	double NumLight;		// Number DENSITY of light quarks at set T
	double NumStrange;		// Number DENSITY of squarks at set T
	double UDDeg = 4.*6.;	// Degeneracy of UD quarks
	double OddDeg = 2.*6.;	// Degeneracy of Squarks
	
	/* Calculated Quantities in Cells */
	double LorBoost[4][4];	// Lorentz boost defined as used- Form is always Lambda_u^v
	int GeneratedParticles;	// Number of particles to be generated this cell
	
	/* Define Parton Densities */
	double cut = 10*T; // Each coordinate of P is integrated this far
	PartCount = 0;
	NumLight = 0;       // Initialize density for Light and strange quarks
	NumStrange = 0;
	GeneratedParticles = 0;
	
	/** Define Rest-frame Cumulative functions**/
	/** See MCSampler.cpp **/
	CDFGen(T, ML, 1); // for light quarks
	CDFGen(T, MS, 2); // for squarks
	//cout << CDFTabLight[512][0] << "  " << CDFTabLight[512][1] << "  " << CDFTabLight[1024][0] << "  " << CDFTabLight[1024][1] << endl;
   
	/*** Open File with Hypersurface data ***/
	/*** Format is presumed t, x, y, sigma_t, sigma_x, sigma_y, temperature, v_x, v_y  ***/
	//HyperSurface read here!
	//Surface.open("hyper_surface_3+1d.dat",ios::in);
    
	//while(Surface.good()){ //File has more data
	for(int iS=0; iS<surface.size(); ++iS){
		
		/*** Position of surface element ***/
//		Surface >> CPos[0]; Surface >> CPos[1]; Surface >> CPos[2]; Surface >> CPos[3];
		/*** Surface element area vector ***/
//		Surface >> LFSigma[0]; Surface >> LFSigma[1]; Surface >> LFSigma[2]; Surface >> LFSigma[3];
		/*** Temperature ***/
//		Surface >> TRead;
		/*** THREE Velocity ***/
//		Surface >> Vel[1]; Surface >> Vel[2]; Surface >> Vel[3];

		double tau_pos, eta_pos, tau_sur, eta_sur;

//		tau_pos = surface[iS].tau; CPos[1] = surface[iS].x; CPos[2] = surface[iS].y; eta_pos = surface[iS].eta;  //we need t,x,y,z
//		//this is also tau, x, y, eta
//		tau_sur = surface[iS].d3sigma_mu[0]; LFSigma[1] = surface[iS].d3sigma_mu[1]; LFSigma[2] = surface[iS].d3sigma_mu[2]; eta_sur = surface[iS].d3sigma_mu[3];
//		TRead = surface[iS].temperature; Vel[1] = surface[iS].vx; Vel[2] = surface[iS].vy; Vel[3] = surface[iS].vz;


		tau_pos = surface[iS][0]; CPos[1] = surface[iS][1]; CPos[2] = surface[iS][2]; eta_pos = surface[iS][3];  //we need t,x,y,z
		//this is also tau, x, y, eta
		tau_sur = surface[iS][4]; LFSigma[1] = surface[iS][5]; LFSigma[2] = surface[iS][6]; eta_sur = surface[iS][7];
		TRead = surface[iS][8]; Vel[1] = surface[iS][9]; Vel[2] = surface[iS][10]; Vel[3] = surface[iS][11];
		
		
		//getting t,z from tau,eta
		CPos[0] = tau_pos*std::cosh(eta_pos); CPos[3] = tau_pos*std::sinh(eta_pos);
		LFSigma[0] = tau_sur*std::cosh(eta_sur); LFSigma[3] = tau_sur*std::sinh(eta_sur);
		
		
		/*** Deduced info ***/
		Vel[0] = 1. / sqrt(1- Vel[1]*Vel[1] - Vel[2]*Vel[2] - Vel[3]*Vel[3] );   // gamma - Vel is not four velocity
		
		//if(Surface.good()){ //Required in while loop to NOT double count last entry
			
			double vsquare = ( Vel[1]*Vel[1] + Vel[2]*Vel[2] + Vel[3]*Vel[3]  );
			
			/*** Lambda_u ^v ***/
			LorBoost[0][0]= Vel[0];
			LorBoost[0][1]= Vel[0]*Vel[1];
			LorBoost[0][2]= Vel[0]*Vel[2];
			LorBoost[0][3]= Vel[0]*Vel[3];
			LorBoost[1][0]= -Vel[0]*Vel[1];
			LorBoost[1][1]= -(Vel[0] - 1.)*Vel[1]*Vel[1]/vsquare - 1.;  
			LorBoost[1][2]= -(Vel[0] - 1.)*Vel[1]*Vel[2]/vsquare;
			LorBoost[1][3]= -(Vel[0] - 1.)*Vel[1]*Vel[3]/vsquare;
			LorBoost[2][0]= -Vel[0]*Vel[2];
			LorBoost[2][1]= -(Vel[0] - 1.)*Vel[1]*Vel[2]/vsquare;
			LorBoost[2][2]= -(Vel[0] - 1.)*Vel[2]*Vel[2]/vsquare - 1.;
			LorBoost[2][3]= -(Vel[0] - 1.)*Vel[2]*Vel[3]/vsquare;
			LorBoost[3][0]= -Vel[0]*Vel[3];
			LorBoost[3][1]= -(Vel[0] - 1.)*Vel[1]*Vel[3]/vsquare;
			LorBoost[3][2]= -(Vel[0] - 1.)*Vel[2]*Vel[3]/vsquare;
			LorBoost[3][3]= -(Vel[0] - 1.)*Vel[3]*Vel[3]/vsquare - 1.;
			
			if(vsquare == 0){
				LorBoost[1][1]= - 1.;  
				LorBoost[1][2]= 0;
				LorBoost[1][3]= 0;
				LorBoost[2][1]= 0;
				LorBoost[2][2]= - 1.;
				LorBoost[2][3]= 0;
				LorBoost[3][1]= 0;
				LorBoost[3][2]= 0;
				LorBoost[3][3]= - 1.;
			}
			/*** Lambda_u^v Sigma_v = CMSigma_u ***/
			CMSigma[0] = (LorBoost[0][0]*LFSigma[0] + LorBoost[0][1]*LFSigma[1] + LorBoost[0][2]*LFSigma[2] + LorBoost[0][3]*LFSigma[3]);
			CMSigma[1] = (LorBoost[1][0]*LFSigma[0] + LorBoost[1][1]*LFSigma[1] + LorBoost[1][2]*LFSigma[2] + LorBoost[1][3]*LFSigma[3]);
			CMSigma[2] = (LorBoost[2][0]*LFSigma[0] + LorBoost[2][1]*LFSigma[1] + LorBoost[2][2]*LFSigma[2] + LorBoost[2][3]*LFSigma[3]);
			CMSigma[3] = (LorBoost[3][0]*LFSigma[0] + LorBoost[3][1]*LFSigma[1] + LorBoost[3][2]*LFSigma[2] + LorBoost[3][3]*LFSigma[3]);
			
			/*** GAUSSIAN INTEGRALS <n> = int f(p)d3p ***/
			NumLight = 0.;
			NumStrange = 0.;
			
			for (int l=0; l<GPoints; l++){
				for (int m=0; m<GPoints; m++){
					for(int k=0; k<GPoints; k++){
						/*** For UD quarks ***/
						NumLight=NumLight+(GWeight[l]*GWeight[m]*GWeight[k]*1./(exp(sqrt(ML*ML+(cut*GAbs[l])*(cut*GAbs[l])+(cut*GAbs[m])*(cut*GAbs[m])+(cut*GAbs[k])*(cut*GAbs[k]))/T)+1.))*
						(sqrt(ML*ML+(cut*GAbs[l])*(cut*GAbs[l])+(cut*GAbs[m])*(cut*GAbs[m])+(cut*GAbs[k])*(cut*GAbs[k]))*CMSigma[0]
						+(cut*GAbs[l])*CMSigma[1]+(cut*GAbs[m])*CMSigma[2]+(cut*GAbs[k])*CMSigma[3])/
						(sqrt(ML*ML+(cut*GAbs[l])*(cut*GAbs[l])+(cut*GAbs[m])*(cut*GAbs[m])+(cut*GAbs[k])*(cut*GAbs[k])));
						
						/*** For Squarks ***/
						NumStrange=NumStrange+(GWeight[l]*GWeight[m]*GWeight[k]*1./(exp(sqrt(MS*MS+(cut*GAbs[l])*(cut*GAbs[l])+(cut*GAbs[m])*(cut*GAbs[m])+(cut*GAbs[k])*(cut*GAbs[k]))/T)+1.))*
						(sqrt(MS*MS+(cut*GAbs[l])*(cut*GAbs[l])+(cut*GAbs[m])*(cut*GAbs[m])+(cut*GAbs[k])*(cut*GAbs[k]))*CMSigma[0]
						+(cut*GAbs[l])*CMSigma[1]+(cut*GAbs[m])*CMSigma[2]+(cut*GAbs[k])*CMSigma[3])/
						(sqrt(MS*MS+(cut*GAbs[l])*(cut*GAbs[l])+(cut*GAbs[m])*(cut*GAbs[m])+(cut*GAbs[k])*(cut*GAbs[k])));
					}
				}
			}
			
			/*** U, D, UBAR, DBAR QUARKS ***/
			/*** <N> = V <n> ***/
			double NumHere=NumLight*UDDeg*cut*cut*cut/(8.*PIE*PIE*PIE);
			
			/* Generating Light Quarks */
			
			GeneratedParticles = 0; //Initialize particle screated this cell
			double ThisRoll = ran(); // Probability of existence die roll
			
			if( ThisRoll >= (exp(-NumHere))){ //Probability rolled > prob of 1 particle
				if( ThisRoll >= ((exp(-NumHere))+ (NumHere*exp(-NumHere)))){ //Probability rolled > prob of 2 particles
					if(ThisRoll>=((exp(-NumHere))+(NumHere*exp(-NumHere))+(0.5*NumHere*NumHere*exp(-NumHere)))){GeneratedParticles+=3;} //Probability rolled > prob of 3 particles
					else{GeneratedParticles+=2;} // 2, but not 3 particles
				}
				else{GeneratedParticles++;} // 1, but not 2 particles
			}
			else{/*** No Particles here ***/}

			//List of Particles ( pos(x,y,z,t), mom(px,py,pz,E), species)
			for(int partic=0; partic<GeneratedParticles; partic++){
				
				//adding space to PList for output quarks
				std::vector<double> temp = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}; Plist.push_back(temp);
				
				/* Species - U,D,UBar,Dbar are equally likely */
				double SpecRoll = ran(); //Probability of species die roll
				if(SpecRoll <=0.25){Plist[PartCount][1] = -2;} // DBar
				else{
					if(SpecRoll <= 0.50){Plist[PartCount][1] = -1;} // UBar
					else{
						if(SpecRoll <=0.75){Plist[PartCount][1] = 1;} // U
						else{Plist[PartCount][1] = 2;} // D
					}
				}
				
				/* Position */
				/*** Located at x,y pos of area element ***/
				Plist[PartCount][10] = CPos[0] + ( ran() - 0.5)*CellDT; // Tau
				Plist[PartCount][7] = CPos[1] + ( ran() - 0.5)*CellDX;
				Plist[PartCount][8] = CPos[2] + ( ran() - 0.5)*CellDY;
				Plist[PartCount][9] = CPos[3] + ( ran() - 0.5)*CellDZ;

				/* Momentum */
				/*** Sample rest frame momentum given T and Mass of light quark ***/
				MCSampler(T, ML, 1); // NewP= P, NewX=Px, ...

				/* USE THE SAME LORBOOST AS BEFORE */
				/*** PLab^u = g^u^t Lambda_t ^v pres^w g_w _v  ***/
				Plist[PartCount][6]= (LorBoost[0][0]*sqrt(ML*ML + NewP*NewP) - LorBoost[0][1]*NewX - LorBoost[0][2]*NewY - LorBoost[0][3]*NewZ)*GEVFM;
				Plist[PartCount][3]=(-LorBoost[1][0]*sqrt(ML*ML + NewP*NewP) + LorBoost[1][1]*NewX + LorBoost[1][2]*NewY + LorBoost[1][3]*NewZ)*GEVFM;
				Plist[PartCount][4]=(-LorBoost[2][0]*sqrt(ML*ML + NewP*NewP) + LorBoost[2][1]*NewX + LorBoost[2][2]*NewY + LorBoost[2][3]*NewZ)*GEVFM;
				Plist[PartCount][5]=(-LorBoost[3][0]*sqrt(ML*ML + NewP*NewP) + LorBoost[3][1]*NewX + LorBoost[3][2]*NewY + LorBoost[3][3]*NewZ)*GEVFM;
				/*** Above Returns P in GeV ***/
				/** Additional information **/
				Plist[PartCount][0]  = 1; //Event ID, to match jet formatting
				Plist[PartCount][2]  = 0; //Origin, to match jet formatting
				Plist[PartCount][11] = 0; //Status- Identifies as thermal quark
				PartCount++;
			}

			/*** S, SBAR QUARKS ***/
			/*** <N> = V <n> ***/
			double NumOddHere=NumStrange*OddDeg*cut*cut*cut/(8.*PIE*PIE*PIE);
			/* Generate Squarks */
			int nL = GeneratedParticles;
			GeneratedParticles = 0; //Initialize particle screated this cell
			ThisRoll = ran();// Probability of existence die roll

			if( ThisRoll >= (exp(-NumOddHere))){
				if( ThisRoll >= ((exp(-NumOddHere))+ (NumOddHere*exp(-NumOddHere)))){ //Probability rolled > prob of 2 particles
					if(ThisRoll>=((exp(-NumOddHere))+(NumOddHere*exp(-NumOddHere))+(0.5*NumOddHere*NumOddHere*exp(-NumOddHere)))){GeneratedParticles+=3;} //Probability rolled > prob of 3 particles
					else{GeneratedParticles+=2;} // 2, but not 3 particles
				}
				else{GeneratedParticles++;} // 1, but not 2 particles
			}
			else{/*Nothing to do here*/}
			int nS = GeneratedParticles;
			//List of Particles ( pos(x,y,z,t), mom(px,py,pz,E), species)
			for(int partic=0; partic<GeneratedParticles; partic++){
				
				//adding space to PList for output quarks
				std::vector<double> temp = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}; Plist.push_back(temp);
				
				/* Species - S,Sbar are equally likely */
				double SpecRoll = ran();
				if(SpecRoll <=0.50000){Plist[PartCount][1] = -3;} // SBar
				else{Plist[PartCount][1] = 3;} //S

				/* Position */
				/*** Located at x,y pos of area element ***/
				Plist[PartCount][10] = CPos[0];
				Plist[PartCount][7] = CPos[1];
				Plist[PartCount][8] = CPos[2];
				Plist[PartCount][9] = CPos[3];

				/* Momentum */
				/*** Sample rest frame momentum given T and Mass of squark ***/
				MCSampler(T, MS, 2); // NewP= P, NewX=Px, ...

				/*** Velocity with which to boost this particle from IT'S frame (off eta=0) ***/ /* USE THE SAME LORBOOST AS BEFORE */

				/*** PLab^u = g^u^t Lambda_t ^v pres^w g_w _v  ***/
				Plist[PartCount][6]=( LorBoost[0][0]*sqrt(MS*MS + NewP*NewP) - LorBoost[0][1]*NewX - LorBoost[0][2]*NewY - LorBoost[0][3]*NewZ)*GEVFM;
				Plist[PartCount][3]=(-LorBoost[1][0]*sqrt(MS*MS + NewP*NewP) + LorBoost[1][1]*NewX + LorBoost[1][2]*NewY + LorBoost[1][3]*NewZ)*GEVFM;
				Plist[PartCount][4]=(-LorBoost[2][0]*sqrt(MS*MS + NewP*NewP) + LorBoost[2][1]*NewX + LorBoost[2][2]*NewY + LorBoost[2][3]*NewZ)*GEVFM;
				Plist[PartCount][5]=(-LorBoost[3][0]*sqrt(MS*MS + NewP*NewP) + LorBoost[3][1]*NewX + LorBoost[3][2]*NewY + LorBoost[3][3]*NewZ)*GEVFM;
				/*** Above Returns P in GeV ***/
				/** Additional information **/
				Plist[PartCount][0]  = 1; //Event ID, to match jet formatting
				Plist[PartCount][2]  = 0; //Origin, to match jet formatting
				Plist[PartCount][11] = 0; //Status- Identifies as thermal quark
				PartCount++;
			}
			
		//	std::cout << iS << " " << nL+nS << " " << Vel[1] << " " << Vel[2] << " " << Vel[3] << " " << CPos[0] << " " << CPos[1] << " " << CPos[2] << " " << CPos[3] << " " << LFSigma[0] << " " << LFSigma[1] << " " << LFSigma[2] << " " << LFSigma[3] << "\n";
		//}
	}
    
//	Surface.close();
	
	//Shuffling PList
	if(ShuffleList){
	for(int i=0;i<PartCount-1;++i){
		int ranelement = i+floor((PartCount-i)*ran());
		double temp[12]; 
		for(int x=0; x<12;x++){temp[x]= Plist[i][x];}
		for(int x=0; x<12;x++){Plist[i][x] = Plist[ranelement][x];}
		for(int x=0; x<12;x++){Plist[ranelement][x] = temp[x];}
	}}
}
