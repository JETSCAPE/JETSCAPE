/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
 * 
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#include "Martini.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include <string>

#include "tinyxml2.h"
#include<iostream>

#include "FluidDynamics.h"
#include "MartiniMutex.h"
#define hbarc 0.197327053

#define MAGENTA "\033[35m"

using namespace Jetscape;

// Register the module with the base class
RegisterJetScapeModule<Martini> Martini::reg("Martini");

using std::ofstream;
using std::ifstream;
using std::ostream;
using std::ios;


Martini::Martini()
{
  SetId("Martini");
  VERBOSE(8);

  //vectors for elastic rates:
  dGamma_qq = new vector<double>;
  dGamma_qg = new vector<double>;
  dGamma_qq_q = new vector<double>;
  dGamma_qg_q = new vector<double>;

  // create and set Martini Mutex
  auto martini_mutex = make_shared<MartiniMutex>();
  SetMutex(martini_mutex);
}

Martini::~Martini()
{
  VERBOSE(8);
}

void Martini::Init()
{
  JSINFO<<"Intialize Martini ...";

  double deltaT = 0.0;
  double Martini_deltaT_Max = 0.01 + rounding_error;

  deltaT = GetXMLElementDouble({"Eloss", "deltaT"});

  if ( deltaT > Martini_deltaT_Max ) {
    JSWARN << "Timestep for Martini ( deltaT = " << deltaT << " ) is too large. "
	 << "Please choose a detaT smaller than or equal to 0.01 in the XML file.";
    throw std::runtime_error("Martini not properly initialized in XML file ...");
  }

  string s = GetXMLElementText({"Eloss", "Martini", "name"});
  JSDEBUG << s << " to be initilizied ...";

  Q0 = GetXMLElementDouble({"Eloss", "Martini", "Q0"});
  alpha_s = GetXMLElementDouble({"Eloss", "Martini", "alpha_s"});
  pcut = GetXMLElementDouble({"Eloss", "Martini", "pcut"});
  hydro_Tc = GetXMLElementDouble({"Eloss", "Martini", "hydro_Tc"});

  g = sqrt(4.*M_PI*alpha_s);
  alpha_em = 1./137.;
  hydro_tStart = 0.6;

  // Path to additional data
  PathToTables = GetXMLElementText({"Eloss", "Martini", "path"});

  // Initialize random number distribution
  ZeroOneDistribution = uniform_real_distribution<double> { 0.0, 1.0 };
  
  readRadiativeRate(&dat, &Gam);
  readElasticRateOmega();
  readElasticRateQ();
}

void Martini::DoEnergyLoss(double deltaT, double Time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut)
{
  VERBOSESHOWER(5)<< MAGENTA << "SentInPartons Signal received : "<<deltaT<<" "<<Q2<<" "<< pIn.size();

  // particle info
  int Id, newId;
  double pAbs, px, py, pz, ee;  // momentum for initial parton (pIn)
  double pRest, pxRest;         // momentum in the rest frame of fluid cell (pIn)
  double pyRest, pzRest;
  double k, kRest;              // momentum for radiated parton (pOut)
  double pNew, pxNew;           // momentum for final parton (pOut)
  double pyNew, pzNew;
  double pNewRest;              // momentum in the rest frame of fluid cell (pOut)
  double omega, q;              // transferred energy/momentum for scattering
  double xx, yy, zz, tt;        // position of initial parton (pIn)
  FourVector pVec, pVecNew;     // 4 vectors for momenta before & after process
  FourVector pVecRest;          // 4 vector in the rest frame of fluid cell
  FourVector pVecNewRest;
  FourVector kVec;              // 4 vector for momentum of radiated particle
  FourVector xVec;              // 4 vector for position (for next time step!)
  double velocity_jet[4];       // jet velocity for MATTER
  double eta;                   // pseudo-rapidity
  
  // flow info
  double vx, vy, vz;            // 3 components of flow velocity
  double T;                     // Temperature of fluid cell
  double beta, gamma;           // flow velocity & gamma factor
  double cosPhi;                // angle between flow and particle
  double cosPhiRest;            // angle between flow and particle in rest frame
  double boostBack;             // factor for boosting back to lab frame
  double cosPhiRestEl;          // angle between flow and scat. particle in rest frame
  double boostBackEl;
  
  for (int i=0;i<pIn.size();i++) {

    // Particle infomration
    Id = pIn[i].pid();

    px = pIn[i].px();
    py = pIn[i].py();
    pz = pIn[i].pz();

    ee = sqrt(px*px+py*py+pz*pz);

    // In MARTINI, particles are all massless and on-shell
    pAbs = sqrt(px*px+py*py+pz*pz);
    pVec = FourVector ( px, py, pz, pAbs );

    tt = pIn[i].x_in().t();
    xx = pIn[i].x_in().x() + (Time-tt)*px/ee;
    yy = pIn[i].x_in().y() + (Time-tt)*py/ee;
    zz = pIn[i].x_in().z() + (Time-tt)*pz/ee;

    eta = pIn[i].eta();

    // Extract fluid properties
    std::unique_ptr<FluidCellInfo> check_fluid_info_ptr;
    GetHydroCellSignal(Time, xx, yy, zz, check_fluid_info_ptr);
    VERBOSE(8)<< MAGENTA<<"Temperature from Brick (Signal) = "
	      <<check_fluid_info_ptr->temperature;

    vx = check_fluid_info_ptr->vx;
    vy = check_fluid_info_ptr->vy;
    vz = check_fluid_info_ptr->vz;
    T = check_fluid_info_ptr->temperature;

    beta = sqrt( vx*vx + vy*vy + vz*vz );

    // Only accept low t particles
    if (pIn[i].t() > Q0*Q0 + rounding_error || Time <=hydro_tStart || T < hydro_Tc)
      continue;
    TakeResponsibilityFor ( pIn[i] ); // Generate error if another module already has responsibility.

    // Set momentum in fluid cell's frame
    // 1: for brick
    if (beta < 1e-10)
      {
	gamma = 1.;
	cosPhi = 1.;
	pRest = pAbs;
	pVecRest = pVec;

	cosPhiRest = 1.;
	boostBack = 1.;
      }
    // 2: for evolving medium
    else
      {
	gamma  = 1./sqrt( 1. - beta*beta );
	cosPhi = ( px*vx + py*vy + pz*vz )/( pAbs*beta );

	// boost particle to the local rest frame of fluid cell
	pRest  = pAbs*gamma*( 1. - beta*cosPhi );

	pxRest = -vx*gamma*pAbs
	  + (1.+(gamma-1.)*vx*vx/(beta*beta))*px
	  + (gamma-1.)*vx*vy/(beta*beta)*py
	  + (gamma-1.)*vx*vz/(beta*beta)*pz;
	pyRest = -vy*gamma*pAbs
	  + (1.+(gamma-1.)*vy*vy/(beta*beta))*py
	  + (gamma-1.)*vx*vy/(beta*beta)*px
	  + (gamma-1.)*vy*vz/(beta*beta)*pz;
	pzRest = -vz*gamma*pAbs
	  + (1.+(gamma-1.)*vz*vz/(beta*beta))*pz
	  + (gamma-1.)*vx*vz/(beta*beta)*px
	  + (gamma-1.)*vy*vz/(beta*beta)*py;

	pVecRest = FourVector ( pxRest, pyRest, pzRest, pRest );

	cosPhiRest = ( pxRest*vx + pyRest*vy + pzRest*vz )/( pRest*beta );
	boostBack = gamma*( 1. + beta*cosPhiRest );
      }

    if (pRest < pcut) continue;

    xVec = FourVector( xx+px/pAbs*deltaT, yy+py/pAbs*deltaT, zz+pz/pAbs*deltaT,
		       Time+deltaT );

    velocity_jet[0]=1.0;
    velocity_jet[1]=pIn[i].jet_v().x();
    velocity_jet[2]=pIn[i].jet_v().y();
    velocity_jet[3]=pIn[i].jet_v().z();

    double deltaTRest = deltaT/gamma;
    int process = DetermineProcess(pRest, T, deltaTRest, Id);
    VERBOSE(8)<< MAGENTA
	      << "Time = " << Time << " Id = " << Id
	      << " process = " << process << " T = " << T
	      << " pAbs = " << pAbs << " " << px << " " << py << " " << pz 
	      << " | position = " << xx << " " << yy << " " << zz;


    // Do nothing for this parton at this timestep
    if (process == 0) 
      {
	pOut.push_back(Parton(0, Id, 0, pVec, xVec));
	pOut[pOut.size()-1].set_form_time(0.);
	pOut[pOut.size()-1].set_jet_v(velocity_jet); // use initial jet velocity

	return;
      }
    if (std::abs(Id) == 1 || std::abs(Id) == 2 || std::abs(Id) == 3)
      {
	// quark radiating gluon (q->qg)
	if (process == 1)
	  {
	    if (pRest/T < AMYpCut) return;

	    // sample radiated parton's momentum
	    kRest = getNewMomentumRad(pRest, T, process);
	    if(kRest > pRest) return;

	    // final state parton's momentum
	    pNewRest = pRest - kRest;

	    // if pNew is smaller than pcut, final state parton is
	    // absorbed into medium
	    if (pNewRest > pcut)
	      {
		pNew = pNewRest*boostBack;
		pVecNew.Set( (px/pAbs)*pNew, (py/pAbs)*pNew, (pz/pAbs)*pNew, pNew );
		pOut.push_back(Parton(0, Id, 0, pVecNew, xVec));
		pOut[pOut.size()-1].set_form_time(0.);
		pOut[pOut.size()-1].set_jet_v(velocity_jet); // use initial jet velocity

	      }

	    if (kRest > pcut)
	      {
		k = kRest*boostBack;
		kVec.Set( (px/pAbs)*k, (py/pAbs)*k, (pz/pAbs)*k, k );
		pOut.push_back(Parton(0, 21, 0, kVec, xVec));
		pOut[pOut.size()-1].set_form_time(0.);
		pOut[pOut.size()-1].set_jet_v(velocity_jet); // use initial jet velocity
	      }

	    return;
	  }
	// quark radiating photon (q->qgamma)
	else if (process == 2)
	  {
	    if (pRest/T < AMYpCut) return;

	    // sample radiated parton's momentum
	    kRest = getNewMomentumRad(pRest, T, process);
	    if(kRest > pRest) return;

	    // final state parton's momentum
	    pNewRest = pRest - kRest;

	    // if pNew is smaller than pcut, final state parton is
	    // absorbed into medium
	    if (pNewRest > pcut)
	      {
		pNew = pNewRest*boostBack;
		pVecNew.Set( (px/pAbs)*pNew, (py/pAbs)*pNew, (pz/pAbs)*pNew, pNew );
		pOut.push_back(Parton(0, Id, 0, pVecNew, xVec));
		pOut[pOut.size()-1].set_form_time(0.);
		pOut[pOut.size()-1].set_jet_v(velocity_jet); // use initial jet velocity
	      }

	    //// photon doesn't have energy threshold; No absorption into medium
	    //// However, we only keep positive energy photons
	    //if (kRest > 0.)
	    //  {
	    //    k = kRest*boostBack;
	    //    kVec.Set( (px/pAbs)*k, (py/pAbs)*k, (pz/pAbs)*k, k );
	    //    pOut.push_back(Parton(0, 22, 0, kVec, xVec));
	    //    pOut[pOut.size()-1].set_form_time(0.);
	    //    pOut[pOut.size()-1].set_jet_v(velocity_jet); // use initial jet velocity
	    //  }

	    return;
	  }
	// quark scattering with either quark (qq->qq) or gluon (qg->qg)
	else if (process == 5 || process == 6)
	  {
	    omega = getEnergyTransfer(pRest, T, process);
	    q = getMomentumTransfer(pRest, omega, T, process);

	    // momentum transfer is always space-like
	    if(q < fabs(omega)) return;

	    pVecNewRest = getNewMomentumElas(pVecRest, omega, q);

	    pNewRest = pVecNewRest.t();

	    // if pNew is smaller than pcut, final state parton is
	    // absorbed into medium
	    if (pNewRest > pcut)
	      {
		// Boost scattered particle to lab frame
		// 1: for brick
		if (beta < 1e-10)
		  {
		    pVecNew = pVecNewRest;
		  }
		// 2: for evolving medium
		else
		  {
		    pxNew = vx*gamma*pVecNewRest.t() 
		      + (1.+(gamma-1.)*vx*vx/(beta*beta))*pVecNewRest.x()
		      + (gamma-1.)*vx*vy/(beta*beta)*pVecNewRest.y()
		      + (gamma-1.)*vx*vz/(beta*beta)*pVecNewRest.z();
		    pyNew = vy*gamma*pVecNewRest.t() 
		      + (1.+(gamma-1.)*vy*vy/(beta*beta))*pVecNewRest.y()
		      + (gamma-1.)*vx*vy/(beta*beta)*pVecNewRest.x()
		      + (gamma-1.)*vy*vz/(beta*beta)*pVecNewRest.z();
		    pzNew = vz*gamma*pVecNewRest.t() 
		      + (1.+(gamma-1.)*vz*vz/(beta*beta))*pVecNewRest.z()
		      + (gamma-1.)*vx*vz/(beta*beta)*pVecNewRest.x()
		      + (gamma-1.)*vy*vz/(beta*beta)*pVecNewRest.y();
			  
		    pNew = sqrt( pxNew*pxNew + pyNew*pyNew + pzNew*pzNew );
		    pVecNew.Set( pxNew, pyNew, pzNew, pNew );
		  }

		pOut.push_back(Parton(0, Id, 0, pVecNew, xVec));
		pOut[pOut.size()-1].set_form_time(0.);
		pOut[pOut.size()-1].set_jet_v(velocity_jet); // use initial jet velocity
	      }

	    return;
	  }
	// quark converting to gluon
	else if (process == 9)
	  {
	    pOut.push_back(Parton(0, 21, 0, pVec, xVec));
	    pOut[pOut.size()-1].set_form_time(0.);
	    pOut[pOut.size()-1].set_jet_v(velocity_jet); // use initial jet velocity

	    return;
	  }
	// quark converting to photon
	else if (process == 10)
	  {
	    pOut.push_back(Parton(0, 22, 0, pVec, xVec));
	    pOut[pOut.size()-1].set_form_time(0.);
	    pOut[pOut.size()-1].set_jet_v(velocity_jet); // use initial jet velocity

	    return;
	  }
      }
    else if (Id == 21)
      {
	// gluon radiating gluon (g->gg)
	if (process == 3)
	  {
	    if (pRest/T < AMYpCut) return;

	    // sample radiated parton's momentum
	    kRest = getNewMomentumRad(pRest, T, process);
	    if(kRest > pRest) return;

	    // final state parton's momentum
	    pNewRest = pRest - kRest;

	    // if pNew is smaller than pcut, final state parton is
	    // absorbed into medium
	    if (pNewRest > pcut)
	      {
		pNew = pNewRest*boostBack;
		pVecNew.Set( (px/pAbs)*pNew, (py/pAbs)*pNew, (pz/pAbs)*pNew, pNew );
		pOut.push_back(Parton(0, Id, 0, pVecNew, xVec));
		pOut[pOut.size()-1].set_form_time(0.);
		pOut[pOut.size()-1].set_jet_v(velocity_jet); // use initial jet velocity
	      }

	    if (kRest > pcut)
	      {
		k = kRest*boostBack;
		kVec.Set( (px/pAbs)*k, (py/pAbs)*k, (pz/pAbs)*k, k );
		pOut.push_back(Parton(0, Id, 0, kVec, xVec));
		pOut[pOut.size()-1].set_form_time(0.);
		pOut[pOut.size()-1].set_jet_v(velocity_jet); // use initial jet velocity
	      }

	    return;
	  }
	// gluon split into quark-antiquark pair (g->qqbar)
	if (process == 4)
	  {
	    if (pRest/T < AMYpCut) return;

	    // sample radiated parton's momentum
	    kRest = getNewMomentumRad(pRest, T, process);
	    if(kRest > pRest) return;

	    // final state parton's momentum
	    pNewRest = pRest - kRest;

	    // choose the Id of new qqbar pair. Note that we only deal with nf = 3
	    double r = ZeroOneDistribution(*GetMt19937Generator());
	    if (r < 1./3.) newId = 1;
	    else if (r < 2./3.) newId = 2;
	    else newId = 3;

	    // if pNew is smaller than pcut, final state parton is
	    // absorbed into medium
	    if (pNewRest > pcut)
	      {
		// *momentum of quark is usually larger than that of anti-quark
		pNew = pNewRest*boostBack;
		pVecNew.Set( (px/pAbs)*pNew, (py/pAbs)*pNew, (pz/pAbs)*pNew, pNew );
		pOut.push_back(Parton(0, newId, 0, pVecNew, xVec));
		pOut[pOut.size()-1].set_form_time(0.);
		pOut[pOut.size()-1].set_jet_v(velocity_jet); // use initial jet velocity
	      }

	    if (kRest > pcut)
	      {
		k = kRest*boostBack;
		kVec.Set( (px/pAbs)*k, (py/pAbs)*k, (pz/pAbs)*k, k );
		pOut.push_back(Parton(0, -newId, 0, kVec, xVec));
		pOut[pOut.size()-1].set_form_time(0.);
		pOut[pOut.size()-1].set_jet_v(velocity_jet); // use initial jet velocity
	      }

	    return;
	  }
	// gluon scattering with either quark (gq->gq) or gluon (gg->gg)
	else if (process == 7 || process == 8)
	  {
	    omega = getEnergyTransfer(pRest, T, process);
	    q = getMomentumTransfer(pRest, omega, T, process);

	    // momentum transfer is always space-like
	    if(q < fabs(omega)) return;

	    pVecNewRest = getNewMomentumElas(pVecRest, omega, q);

	    pNewRest = pVecNewRest.t();

	    // if pNew is smaller than pcut, final state parton is
	    // absorbed into medium
	    if (pNewRest > pcut)
	      {
		// Boost scattered particle to lab frame
		// 1: for brick
		if (beta < 1e-10)
		  {
		    pVecNew = pVecNewRest;
		  }
		// 2: for evolving medium
		else
		  {
		    pxNew = vx*gamma*pVecNewRest.t()
		      + (1.+(gamma-1.)*vx*vx/(beta*beta))*pVecNewRest.x()
		      + (gamma-1.)*vx*vy/(beta*beta)*pVecNewRest.y()
		      + (gamma-1.)*vx*vz/(beta*beta)*pVecNewRest.z();
		    pyNew = vy*gamma*pVecNewRest.t()
		      + (1.+(gamma-1.)*vy*vy/(beta*beta))*pVecNewRest.y()
		      + (gamma-1.)*vx*vy/(beta*beta)*pVecNewRest.x()
		      + (gamma-1.)*vy*vz/(beta*beta)*pVecNewRest.z();
		    pzNew = vz*gamma*pVecNewRest.t() 
		      + (1.+(gamma-1.)*vz*vz/(beta*beta))*pVecNewRest.z()
		      + (gamma-1.)*vx*vz/(beta*beta)*pVecNewRest.x()
		      + (gamma-1.)*vy*vz/(beta*beta)*pVecNewRest.y();
			  
		    pNew = sqrt( pxNew*pxNew + pyNew*pyNew + pzNew*pzNew );
		    pVecNew.Set( pxNew, pyNew, pzNew, pNew );
		  }

		pOut.push_back(Parton(0, Id, 0, pVecNew, xVec));
		pOut[pOut.size()-1].set_form_time(0.);
		pOut[pOut.size()-1].set_jet_v(velocity_jet); // use initial jet velocity
	      }

	    return;
	  }
	// gluon converting to quark
	else if (process == 11)
	  {
	    // choose the Id of new qqbar pair. Note that we only deal with nf = 3
	    double r = ZeroOneDistribution(*GetMt19937Generator());
	    if (r < 1./3.) newId = 1;
	    else if (r < 2./3.) newId = 2;
	    else newId = 3;

	    double antiquark = ZeroOneDistribution(*GetMt19937Generator());
            if(antiquark < 0.5) newId *= -1;

	    pOut.push_back(Parton(0, newId, 0, pVec, xVec));
	    pOut[pOut.size()-1].set_form_time(0.);
	    pOut[pOut.size()-1].set_jet_v(velocity_jet); // use initial jet velocity

	    return;
	  }
      } // Id==21
  } // particle loop
}

int Martini::DetermineProcess(double pRest, double T, double deltaTRest, int Id)
{

  double dT = deltaTRest/hbarc;   // time step in [GeV^(-1)]

  // get the rates for each process
  // total Probability = dT*Rate
  RateRadiative rateRad;
  rateRad = getRateRadTotal(pRest, T);
  RateElastic rateElas;
  rateElas = getRateElasTotal(pRest, T);
  RateConversion rateConv;
  rateConv = getRateConv(pRest, T);

  // evolution for quark (u, d, s)
  if (std::abs(Id) == 1 || std::abs(Id) == 2 || std::abs(Id) == 3)
    {
      // multiplying by (ef/e)^2
      if (abs(Id) == 1)
	{
	  rateRad.qqgamma *= 4./9.;
	  rateConv.qgamma *= 4./9.;
	}
      else
	{
	  rateRad.qqgamma *= 1./9.;
	  rateConv.qgamma *= 1./9.;
	}

      double totalQuarkProb = 0.;

      /* block the photon process at this moment */
      //if (pRest/T > AMYpCut) totalQuarkProb += (rateRad.qqg + rateRad.qqgamma)*dT;
      //totalQuarkProb += (rateElas.qq + rateElas.qg + rateConv.qg + rateConv.qgamma)*dT;

      if (pRest/T > AMYpCut) totalQuarkProb += rateRad.qqg*dT;
      totalQuarkProb += (rateElas.qq + rateElas.qg + rateConv.qg)*dT;

      // warn if total probability exceeds 1
      if (totalQuarkProb > 1.){
	JSWARN << " : Total Probability for quark processes exceeds 1 ("
	     << totalQuarkProb << "). "
	     << " : Most likely this means you should choose a smaller deltaT in the xml"
	     << " (e.g. 0.01).";
	//throw std::runtime_error ("Martini probability problem.");
      }

      double accumProb = 0.;
      double nextProb = 0.;
      double Prob = 0.;

      if (ZeroOneDistribution(*GetMt19937Generator()) < totalQuarkProb)
	{
	  /* label for process
	     [1-4  : Radiation ] 1: q->qg , 2 : q->qgamma, 3 : g->gg , 4: g->qqbar
	     [5-8  : Elastic   ] 5: qq->qq, 6 : qg->qg   , 7 : gq->gq, 8: gg->gg
	     [9-11 : Conversion] 9: q->g  , 10: q->gamma , 11: g->q                */
	  double randProb = ZeroOneDistribution(*GetMt19937Generator());

	  // AMY radiation only happens if energy scale is above certain threshold.
	  // but elastic/conversion processes doesn't have threshold.
	  if(pRest/T > AMYpCut)
	    {
	      Prob = rateRad.qqg*dT/totalQuarkProb;
	      if (accumProb <= randProb && randProb < (accumProb + Prob))
		return 1;

	      //accumProb += Prob;
	      //Prob = rateRad.qqgamma*dT/totalQuarkProb;
	      //if (accumProb <= randProb && randProb < (accumProb + Prob))
	      //  return 2;
	    }

	  accumProb += Prob;
	  Prob = rateElas.qq*dT/totalQuarkProb;
	  if (accumProb <= randProb && randProb < (accumProb + Prob))
	    return 5;

	  accumProb += Prob;
	  Prob = rateElas.qg*dT/totalQuarkProb;
	  if (accumProb <= randProb && randProb < (accumProb + Prob))
	    return 6;

	  accumProb += Prob;
	  Prob = rateConv.qg*dT/totalQuarkProb;
	  if (accumProb <= randProb && randProb < (accumProb + Prob))
	    return 9;

	  //accumProb += Prob;
	  //Prob = rateConv.qgamma*dT/totalQuarkProb;
	  //if (accumProb <= randProb && randProb < (accumProb + Prob))
	  //  return 10;
	}
      else
	{
	  // nothing happens to quark
	  return 0;
	}
    }
  // evolution for gluon
  else if (Id == 21)
    {
      double totalGluonProb = 0.;

      if (pRest/T > AMYpCut) totalGluonProb += (rateRad.ggg + rateRad.gqq)*dT;
      totalGluonProb += (rateElas.gq + rateElas.gg + rateConv.gq)*dT;

      // warn if total probability exceeds 1
      if (totalGluonProb > 1.){
	JSWARN << " : Total Probability for gluon processes exceeds 1 ("
	     << totalGluonProb << "). "
	     << " : Most likely this means you should choose a smaller deltaT in the xml"
	     << " (e.g. 0.01).";
	//throw std::runtime_error ("Martini probability problem.");
      }

      double accumProb = 0.;
      double nextProb = 0.;
      double Prob = 0.;

      if (ZeroOneDistribution(*GetMt19937Generator()) < totalGluonProb)
	{
	  /* label for process
	     [1-4  : Radiation ] 1: q->qg, 2 : q->qgamma, 3 : g->gg, 4: g->qq
	     [5-8  : Elastic   ] 5: q->q , 6 : q->g     , 7 : g->q , 8: g->g
	     [9-11 : Conversion] 9: q->g , 10: q->gamma , 11: g->q            */
	  double randProb = ZeroOneDistribution(*GetMt19937Generator());

	  // AMY radiation only happens if energy scale is above certain threshold.
	  // but elastic/conversion processes doesn't have threshold.
	  if (pRest/T > AMYpCut)
	    {
	      Prob = rateRad.ggg*dT/totalGluonProb;
	      if (accumProb <= randProb && randProb < (accumProb + Prob))
		return 3;

	      accumProb += Prob;
	      Prob = rateRad.gqq*dT/totalGluonProb;
	      if (accumProb <= randProb && randProb < (accumProb + Prob))
		return 4;
	    }

	  accumProb += Prob;
	  Prob = rateElas.gq*dT/totalGluonProb;
	  if (accumProb <= randProb && randProb < (accumProb + Prob))
	    return 7;

	  accumProb += Prob;
	  Prob = rateElas.gg*dT/totalGluonProb;
	  if (accumProb <= randProb && randProb < (accumProb + Prob))
	    return 8;

	  accumProb += Prob;
	  Prob = rateConv.gq*dT/totalGluonProb;
	  if (accumProb <= randProb && randProb < (accumProb + Prob))
	    return 11;
	}
      else
	{
	  // nothing happens to gluon
	  return 0;
	}
    }

  // if parton is other than u,d,s,g, do nothing
  return 0;
}

RateRadiative Martini::getRateRadTotal(double pRest, double T)
{
  RateRadiative rate;

  double u = pRest/T;  // making arguments in log to be dimensionless

  if(u > AMYpCut) 
    {
      rate.qqg     = (0.8616 - 3.2913/(u*u) + 2.1102/u - 0.9485/sqrt(u))*pow(g, 4.)*T;
      rate.ggg     = (1.9463 + 61.7856/(u*u*u) - 30.7877/(u*u) + 8.0409/u - 2.6249/sqrt(u))
            	 *pow(g, 4.)*T;
      rate.gqq     = (2.5830/(u*u*u) - 1.7010/(u*u) + 1.4977/u - 1.1961/pow(u,0.8)
            	+ 0.1807/sqrt(u))*pow(g, 4.)*T*nf; 
      rate.qqgamma = (0.0053056 + 2.3279/pow(u,3.) - 0.6676/u + 0.3223/sqrt(u))
            	 *pow(g, 4.)*alpha_em/alpha_s*T;

      double runningFactor = log(g*T*pow(10., 0.25)/.175)/log(g*T*pow(u, 0.25)/.175);
      if (runningFactor < 1.)
        {
          rate.qqg *= runningFactor;
          rate.gqq *= runningFactor;
          rate.ggg *= runningFactor;
          rate.qqgamma *= runningFactor;
        }
    }
  else
    {
      rate.qqg = 0.;
      rate.ggg = 0.;
      rate.gqq = 0.;
      rate.qqgamma = 0.;
    }

  return rate;
}

RateRadiative Martini::getRateRadPos(double u, double T)
{
  RateRadiative rate;

  rate.qqg = (0.5322 - 3.1037/(u*u) + 2.0139/u - 0.9417/sqrt(u))*pow(g, 4.)*T;
  rate.ggg = (1.1923 - 11.5250/(u*u*u) + 3.3010/u - 1.9049/sqrt(u))*pow(g, 4.)*T;
  rate.gqq = (0.0004656 - 0.04621/(u*u) + 0.0999/u - 0.08171/pow(u,0.8) 
	      + 0.008090/pow(u,0.2) - 1.2525*pow(10.,-8.)*u)*pow(g, 4.)*T*nf; 
  rate.qqgamma = 0.;

  return rate;
}

RateRadiative Martini::getRateRadNeg(double u, double T)
{
  RateRadiative rate;

  rate.qqg = (0.3292 - 0.6759/(u*u) + 0.4871/pow(u,1.5) - 0.05393/u + 0.007878/sqrt(u))
	     *pow(g, 4.)*T;
  rate.ggg = (0.7409 + 1.8608/(u*u*u) - 0.1353/(u*u) + 0.1401/u)*pow(g, 4.)*T;
  rate.gqq = (0.03215/(u*u*u) + 0.01419/(u*u) + 0.004338/u - 0.00001246/sqrt(u))
	     *pow(g, 4.)*T*nf;
  rate.qqgamma = 0.;

  return rate;
}

double Martini::getNewMomentumRad(double pRest, double T, int process)
{
  double u = pRest/T; // making arguments in log to be dimensionless

  double kNew = 0.;   // momentum of radiated gluon (dimentionless)
  double y;           // kNew candidate
  double x;           // random number, uniform on [0,1]
  double randA;       // uniform random number on [0, Area under the envelop function]
  double fy;          // total area under the envelop function
  double fyAct;       // actual rate

  RateRadiative Pos, Neg;
  Pos = getRateRadPos(u, T);
  Neg = getRateRadNeg(u, T);

  // this switch will hold the decision whether k is positive or negative:
  // 0 : negative, 1 : positive
  int posNegSwitch = 1; 

  /* process == 1 : quark radiating gluon
     process == 2 : quark radiating photon
     process == 3 : gluon radiating gluon
     process == 4 : gluon split into quark-antiquark pair */

  if (process == 1)
    {
      // decide whether k shall be positive or negative 
      // if x (uniform on [0,1]) < area(k<0)/area(all k) then k < 0
      if (ZeroOneDistribution(*GetMt19937Generator()) < Neg.qqg/(Neg.qqg+Pos.qqg))
	posNegSwitch = 0;

      if (posNegSwitch == 1) // if k > 0
	{
	  do
	    {
	      randA = ZeroOneDistribution(*GetMt19937Generator())
		     *area(u+12., u, posNegSwitch, 1);
	      y = 2.5/(LambertW(2.59235*pow(10.,23.)*exp(-100.*randA)));

	      fy = 0.025/(y*y)+0.01/y;
	      fyAct = function(u, y, process);

	      x = ZeroOneDistribution(*GetMt19937Generator());

	    } while (x > fyAct/fy); 
	  // reject if x is larger than the ratio fyAct/fy
	  kNew = y;
	}
      else // if k < 0
	{
	  do
	    {
	      randA = ZeroOneDistribution(*GetMt19937Generator())
		     *area(-0.05, u, posNegSwitch, 1);
	      y = -12./(1.+480.*randA);

	      fy = 0.025/(y*y);
	      fyAct = function(u, y, process);

	      x = ZeroOneDistribution(*GetMt19937Generator());

	    } while (x > fyAct/fy); 
	  // reject if x is larger than the ratio fyAct/fy
	  kNew = y;
	}
    }
  else if (process == 2)
    {
      do
	{
	  randA = ZeroOneDistribution(*GetMt19937Generator())
		 *area(1.15*u, u, posNegSwitch, 2);
	  y = 83895.3*pow(pow(u, 0.5)*randA, 10./3.);

	  fy = (0.01/(pow(y, 0.7)))/pow(u, 0.5);
	  fyAct = function(u, y, process);

	  x = ZeroOneDistribution(*GetMt19937Generator());

	} while (x > fyAct/fy); 
      // reject if x is larger than the ratio fyAct/fy
      kNew = y;
    }
  else if (process == 3)
    {
      // decide whether k shall be positive or negative 
      // if x (uniform on [0,1]) < area(k<0)/area(all k) then k < 0
      if (ZeroOneDistribution(*GetMt19937Generator()) < Neg.ggg/(Neg.ggg+Pos.ggg))
	posNegSwitch = 0;

      if( posNegSwitch == 1 ) // if k > 0
	{
	  do
	    {
	      randA = ZeroOneDistribution(*GetMt19937Generator())
		     *area(u/2., u, posNegSwitch, 3);
	      y = 5./(LambertW(2.68812*pow(10., 45.)*exp(-50.*randA)));

	      fy = 0.1/(y*y)+0.02/y;
	      fyAct = function(u, y, process);

	      x = ZeroOneDistribution(*GetMt19937Generator());

	    } while (x > fyAct/fy); 
	  // reject if x is larger than the ratio fyAct/fy
	  kNew = y;
	}
      else // if k < 0
	{
	  do
	    {
	      randA = ZeroOneDistribution(*GetMt19937Generator())
		     *area(-0.05, u, posNegSwitch, 3);
	      y = -12./(1. + 120.*randA);

	      fy = 0.1/(y*y);
	      fyAct = function(u, y, process);

	      x = ZeroOneDistribution(*GetMt19937Generator());

	    } while(x > fyAct/fy); 
	  // reject if x is larger than the ratio fyAct/fy
	  kNew = y;
	}
    }
  else if (process == 4)
    {
      // decide whether k shall be positive or negative 
      // if x (uniform on [0,1]) < area(k<0)/area(all k) then k < 0
      if (ZeroOneDistribution(*GetMt19937Generator()) < Neg.gqq/(Neg.gqq+Pos.gqq))
	posNegSwitch = 0;

      if( posNegSwitch == 1 ) // if k > 0
	{
	  do
	    {
	      randA = ZeroOneDistribution(*GetMt19937Generator())
		     *area(u/2., u, posNegSwitch, 4);
	      y = 0.83333*(0.06*function(u, 0.05, process)+randA)
		 /function(u, 0.05, process);

	      fy = 1.2*function(u, 0.05, process);
	      fyAct = function(u, y, process);

	      x = ZeroOneDistribution(*GetMt19937Generator());

	    } while (x > fyAct/fy); 
	  // reject if x is larger than the ratio fyAct/fy
	  kNew = y;
	}
      else // if k < 0
	{
	  do
	    {
	      randA = ZeroOneDistribution(*GetMt19937Generator())
		     *area(-0.05, u, posNegSwitch, 4);
	      y = (2.5-u*log(7.81082*pow(10., -6.)*exp(14.5/u)
			  +(-115.883+113.566*u)*randA))/(1.-0.98*u);

	      fy = 0.98*exp((1.-1./u)*(-2.5+y))/u;
	      fyAct = function(u, y, process);

	      x = ZeroOneDistribution(*GetMt19937Generator());

	    } while (x > fyAct/fy); 
	  // reject if x is larger than the ratio fyAct/fy
	  kNew = y;
	}
    }
  else
    {
      JSWARN << "Invalid process number (" << process << ")";
    }

  return kNew*T;  // kNew*T is in [GeV]
}

// calculates the area under the envelope function when using the rejection method
// integrals had been solved analytically before
double Martini::area(double y, double u, int posNegSwitch, int process)
{
  if (process == 1)
    {
      if (posNegSwitch == 1)
	return (0.5299 - 0.025/y + 0.01*log(y));
      else 
	return (-0.002083-0.025/y);
    }
  else if (process == 2)
    {
      return ((0.03333*pow(y,0.3))/pow(u,0.5));
    }
  else if (process == 3)
    {
      if (posNegSwitch == 1)
	return (2.05991 - 0.1/y + 0.02*log(y));
      else 
	return (-0.008333 - 0.1/y);
    }      
  else if (process == 4)
    {
      if (posNegSwitch == 1)
	return (1.2*function(u, 0.05, process)*(y-0.05));
      else 
	return ((6.8778*pow(10., -8.)*exp(14.5/u)
		 -0.008805*exp((2.5-y+0.98*u*y)/u))/(1.0204-u));
    }      

  return 0.;
}

double Martini::function(double u, double y, int process)
{
  if (process == 1)      return getRate_qqg(u, y);
  else if (process == 2) return getRate_qqgamma(u, y);
  else if (process == 3) return getRate_ggg(u, y);
  else if (process == 4) return getRate_gqq(u, y);

  return 0.;
}

RateElastic Martini::getRateElasTotal(double pRest, double T)
{
  // compute the total transition rate in GeV, integrated over k, omega
  // and q and angles for fixed E=p and temperature T
  // using parametrization of numerically computed integral
  // then interpolate to get right value for used alpha_s 
  // IMPORTANT: all computed values below are for a minimal omega of 0.05*T
  // so this is the cutoff to use in the calculation also!
  // also Nf=3 was used ... scales out though

  RateElastic rate;

  double u = pRest/T;  // making arguments in log to be dimensionless

  double alpha0 = 0.15;
  double deltaAlpha = 0.03;
  double iAlpha = floor((alpha_s-alpha0)/deltaAlpha+0.001);
  double alphaFrac = (alpha_s-alpha0)/deltaAlpha - iAlpha;
  double rateLower, rateUpper;
  double c[2][6], d[2][6];

  for (int i=0; i<2; i++)
    for (int j=0; j<6; j++)
      {
	c[i][j] = 0.;
	d[i][j] = 0.;
      }

  if ( alpha_s >= 0.15 && alpha_s < 0.18 )
    {
      c[0][0] = 0.18172488396136807;   c[1][0] = 0.224596478395945;
      c[0][1] = 0.6004740049060965;    c[1][1] = 1.0874259848101948;
      c[0][2] = 0.36559627257898347;   c[1][2] = 0.6436398538984057;
      c[0][3] = 0.10607576568373664;   c[1][3] = 0.11585154613692052;
      c[0][4] = 0.004322466954618182;  c[1][4] = -0.001719701730785056;
      c[0][5] = 0.04731599462749122;   c[1][5] = 0.06734745496415469;
    }
  else if ( alpha_s >= 0.18 && alpha_s < 0.21 )
    {
      c[0][0] = 0.224596478395945;     c[1][0] = 0.2686436092048326;
      c[0][1] = 1.0874259848101948;    c[1][1] = 1.7286136256785387;
      c[0][2] = 0.6436398538984057;    c[1][2] = 0.9826325498183079;
      c[0][3] = 0.11585154613692052;   c[1][3] = 0.13136670133029682;
      c[0][4] = -0.001719701730785056; c[1][4] = -0.004876376882437649;
      c[0][5] = 0.06734745496415469;   c[1][5] = 0.09140316977554151;
    }
  else if ( alpha_s >= 0.21 && alpha_s < 0.24 )
    {
      c[0][0] = 0.2686436092048326;    c[1][0] = 0.3137234778163784;
      c[0][1] = 1.7286136256785387;    c[1][1] = 2.445764079999846;
      c[0][2] = 0.9826325498183079;    c[1][2] = 1.3083241146035964;
      c[0][3] = 0.13136670133029682;   c[1][3] = 0.18341717903923757;
      c[0][4] = -0.004876376882437649; c[1][4] = 0.006098371807040589;
      c[0][5] = 0.09140316977554151;   c[1][5] = 0.12054238276023879;
    }
  else if ( alpha_s >= 0.24 && alpha_s < 0.27 )
    {
      c[0][0] = 0.3137234778163784;    c[1][0] = 0.3597255453974444;
      c[0][1] = 2.445764079999846;     c[1][1] = 3.140669321831845;
      c[0][2] = 1.3083241146035964;    c[1][2] = 1.535549334026633;
      c[0][3] = 0.18341717903923757;   c[1][3] = 0.30505450230754705;
      c[0][4] = 0.006098371807040589;  c[1][4] = 0.04285103618362223;
      c[0][5] = 0.12054238276023879;   c[1][5] = 0.1558288379712527;
    }
  else if ( alpha_s >= 0.27 && alpha_s < 0.3 )
    {
      c[0][0] = 0.3597255453974444;    c[1][0] = 0.40656130602563223;
      c[0][1] = 3.140669321831845;     c[1][1] = 3.713430971987352;
      c[0][2] = 1.535549334026633;     c[1][2] = 1.5818298058630476;
      c[0][3] = 0.30505450230754705;   c[1][3] = 0.5269042544852683;
      c[0][4] = 0.04285103618362223;   c[1][4] = 0.11594975218839362;
      c[0][5] = 0.1558288379712527;    c[1][5] = 0.1982063104156748;
    }
  else if ( alpha_s >= 0.3 && alpha_s < 0.33 )
    {
      c[0][0] = 0.40656130602563223;   c[1][0] = 0.45415805200862863;
      c[0][1] = 3.713430971987352;     c[1][1] = 4.0758813206143785;
      c[0][2] = 1.5818298058630476;    c[1][2] = 1.3775134184861555;
      c[0][3] = 0.5269042544852683;    c[1][3] = 0.873527536823307;
      c[0][4] = 0.11594975218839362;   c[1][4] = 0.23371456949506658;
      c[0][5] = 0.1982063104156748;    c[1][5] = 0.24840524848507203;
    }
  else if ( alpha_s >= 0.33 && alpha_s < 0.36 )
    {
      c[0][0] = 0.45415805200862863;   c[1][0] = 0.5024541413891354;
      c[0][1] = 4.0758813206143785;    c[1][1] = 4.159425815179756;
      c[0][2] = 1.3775134184861555;    c[1][2] = 0.8719749565879445;
      c[0][3] = 0.873527536823307;     c[1][3] = 1.3606690530660879;
      c[0][4] = 0.23371456949506658;   c[1][4] = 0.4010658149846402;
      c[0][5] = 0.24840524848507203;   c[1][5] = 0.3067901992139913;
    }
  else if ( alpha_s >= 0.36 && alpha_s < 0.39 )
    {
      c[0][0] = 0.5024541413891354;    c[1][0] = 0.5513999693402064;
      c[0][1] = 4.159425815179756;     c[1][1] = 3.893153859527746;
      c[0][2] = 0.8719749565879445;    c[1][2] = 0.009578762778659829;
      c[0][3] = 1.3606690530660879;    c[1][3] = 2.0095157488463244;
      c[0][4] = 0.4010658149846402;    c[1][4] = 0.6260756501912864;
      c[0][5] = 0.3067901992139913;    c[1][5] = 0.37424991045026396;
    }
  else if ( alpha_s >= 0.39 && alpha_s <= 0.42 )
    {
      c[0][0] = 0.5513999693402064;    c[1][0] = 0.600941593540798;
      c[0][1] = 3.893153859527746;     c[1][1] = 3.293344337592684;
      c[0][2] = 0.009578762778659829;  c[1][2] = -1.1764805445298645;
      c[0][3] = 2.0095157488463244;    c[1][3] = 2.792180001243466;
      c[0][4] = 0.6260756501912864;    c[1][4] = 0.8949534049225013;
      c[0][5] = 0.37424991045026396;   c[1][5] = 0.44878529934031575;
    }

  rateLower = T*(c[0][0]             + c[0][1]/pow(u, 4.) - 
		 c[0][2]/pow(u, 3.)  - c[0][3]/pow(u, 2.) + 
		 c[0][4]/pow(u, 1.5) - c[0][5]/u);
  rateUpper = T*(c[1][0]             + c[1][1]/pow(u, 4.) - 
		 c[1][2]/pow(u, 3.)  - c[1][3]/pow(u, 2.) + 
		 c[1][4]/pow(u, 1.5) - c[1][5]/u);

  rate.qq = (1.-alphaFrac)*rateLower+alphaFrac*rateUpper;
  rate.qq *= nf/3.; // adjust number of flavors

  rate.gq = rate.qq*9./4.;

  if ( alpha_s >= 0.15 && alpha_s < 0.18 )
    {
      d[0][0] = 0.9364689080337059;    d[1][0] = 1.1485486950080581;
      d[0][1] = 2.626076478553979;     d[1][1] = 4.993647646894147;
      d[0][2] = 2.1171556605834274;    d[1][2] = 3.7295251994302876;
      d[0][3] = 0.13123339226210134;   d[1][3] = -0.0017620287506503757;
      d[0][4] = 0.02875811664147147;   d[1][4] = 0.010598257485913224;
      d[0][5] = 0.27736469898722244;   d[1][5] = 0.3949856219367327;
    }
  else if ( alpha_s >= 0.18 && alpha_s < 0.21 )
    {
      d[0][0] = 1.1485486950080581;    d[1][0] = 1.3645568637616001;
      d[0][1] = 4.993647646894147;     d[1][1] = 8.174225869366722;
      d[0][2] = 3.7295251994302876;    d[1][2] = 5.732101892684938;
      d[0][3] = -0.001762028750650376; d[1][3] = -0.1416811579957863;
      d[0][4] = 0.010598257485913224;  d[1][4] = 0.011703596451947428;
      d[0][5] = 0.3949856219367327;    d[1][5] = 0.5354757997870718;
    }
  else if ( alpha_s >= 0.21 && alpha_s < 0.24 )
    {
      d[0][0] = 1.3645568637616001;    d[1][0] = 1.5839378568555678;
      d[0][1] = 8.174225869366722;     d[1][1] = 11.785897000063443;
      d[0][2] = 5.732101892684938;     d[1][2] = 7.758388282689373;
      d[0][3] = -0.1416811579957863;   d[1][3] = -0.13163385415183002;
      d[0][4] = 0.011703596451947428;  d[1][4] = 0.09016386041913003;
      d[0][5] = 0.5354757997870718;    d[1][5] = 0.7042577279136836;
    }
  else if ( alpha_s >= 0.24 && alpha_s < 0.27 )
    {
      d[0][0] = 1.5839378568555678;    d[1][0] = 1.8062676019060235;
      d[0][1] = 11.785897000063443;    d[1][1] = 15.344112642069764;
      d[0][2] = 7.758388282689373;     d[1][2] =  9.384190917330093;
      d[0][3] = -0.13163385415183002;  d[1][3] = 0.19709400976261568;
      d[0][4] = 0.09016386041913003;   d[1][4] = 0.30577623140224813;
      d[0][5] = 0.7042577279136836;    d[1][5] = 0.9066501895009754;
    }
  else if ( alpha_s >= 0.27 && alpha_s < 0.3 )
    {
      d[0][0] = 1.8062676019060235;    d[1][0] = 2.0312125903238236;
      d[0][1] = 15.344112642069764;    d[1][1] = 18.36844006721506;
      d[0][2] = 9.384190917330093;     d[1][2] = 10.209988454804193;
      d[0][3] = 0.19709400976261568;   d[1][3] = 0.9957025988944573;
      d[0][4] = 0.30577623140224813;   d[1][4] =  0.7109302867706849;
      d[0][5] = 0.9066501895009754;    d[1][5] = 1.1472148515742653;
    }
  else if ( alpha_s >= 0.3 && alpha_s < 0.33 )
    {
      d[0][0] = 2.0312125903238236;    d[1][0] = 2.258502734110078;
      d[0][1] = 18.36844006721506;     d[1][1] = 20.43444928479894;
      d[0][2] = 10.209988454804193;    d[1][2] = 9.896928897847518;
      d[0][3] = 0.9957025988944573;    d[1][3] = 2.3867073785159003;
      d[0][4] = 0.7109302867706849;    d[1][4] = 1.3473328178504662;
      d[0][5] = 1.1472148515742653;    d[1][5] = 1.429497460496924;
    }
  else if ( alpha_s >= 0.33 && alpha_s < 0.36 )
    {
      d[0][0] = 2.258502734110078;     d[1][0] = 2.4879110920956653;
      d[0][1] = 20.43444928479894;     d[1][1] = 21.220550462966102;
      d[0][2] = 9.896928897847518;     d[1][2] = 8.20639681844989;
      d[0][3] = 2.3867073785159003;    d[1][3] = 4.445222616370339;
      d[0][4] = 1.3473328178504662;    d[1][4] = 2.2381176005506016;
      d[0][5] = 1.429497460496924;     d[1][5] = 1.7550164762706189;
    }
  else if ( alpha_s >= 0.36 && alpha_s < 0.39 )
    {
      d[0][0] = 2.4879110920956653 ;   d[1][0] = 2.7192501243929903;
      d[0][1] = 21.220550462966102;    d[1][1] = 20.470583876561985;
      d[0][2] = 8.20639681844989 ;     d[1][2] = 4.954737209403953;
      d[0][3] = 4.445222616370339;     d[1][3] = 7.227667929705693;
      d[0][4] = 2.2381176005506016;    d[1][4] = 3.401378906197122;
      d[0][5] = 1.7550164762706189;    d[1][5] = 2.1251383942923474;
    }
  else if ( alpha_s >= 0.39 && alpha_s <= 0.42 )
    {
      d[0][0] = 2.7192501243929903 ;   d[1][0] = 2.9523522354248817;
      d[0][1] = 20.470583876561985;    d[1][1] = 18.027772799078463;
      d[0][2] = 4.954737209403953;     d[1][2] = 0.050298242947981846;
      d[0][3] = 7.227667929705693;     d[1][3] = 10.747352232336384;
      d[0][4] = 3.401378906197122;     d[1][4] = 4.8378133911595285;
      d[0][5] = 2.1251383942923474;    d[1][5] = 2.5391647730624003;
    }

  rateLower = T*(d[1][0]             + d[0][1]/pow(u, 4.) -
		 d[0][2]/pow(u, 3.)  - d[0][3]/pow(u, 2.) +
		 d[0][4]/pow(u, 1.5) - d[0][5]/u);
  rateUpper = T*(d[1][0]             + d[1][1]/pow(u, 4.) -
		 d[1][2]/pow(u, 3.)  - d[1][3]/pow(u, 2.) +
		 d[1][4]/pow(u, 1.5) - d[1][5]/u);
      
  rate.qg = (1.-alphaFrac)*rateLower+alphaFrac*rateUpper;
  rate.qg /= 3.; // historic reasons
  
  rate.gg = rate.qg*9./4.;

  return rate;
}

RateElastic Martini::getRateElasPos(double u, double T)
{
  RateElastic rate;

  double alpha0 = 0.15;
  double deltaAlpha = 0.03;
  double iAlpha = floor((alpha_s-alpha0)/deltaAlpha+0.001);
  double alphaFrac = (alpha_s-alpha0)/deltaAlpha - iAlpha;
  double rateLower, rateUpper;
  double c[2][6], d[2][6];

  for (int i=0; i<2; i++)
    for (int j=0; j<6; j++)
      {
	c[i][j] = 0.;
	d[i][j] = 0.;
      }

  if ( alpha_s >= 0.15 && alpha_s < 0.18 )
    {
      c[0][0] = 0.12199410313320332;   c[1][0] = 0.15243607717720586;
      c[0][1] = 0.23732051765097376;   c[1][1] = 0.5403120875137825;
      c[0][2] = -0.03285419708803458;  c[1][2] = 0.06440920730334501;
      c[0][3] = 0.2255419254079952;    c[1][3] = 0.2881594349535524;
      c[0][4] = 0.03991522899907729;   c[1][4] = 0.04948438583750772;
      c[0][5] = 0.05022641428394594;   c[1][5] = 0.07152523367501308;
    }
  else if ( alpha_s >= 0.18 && alpha_s < 0.21 )
    {
      c[0][0] = 0.15243607717720586;   c[1][0] = 0.15243607717720586;
      c[0][1] = 0.5403120875137825;    c[1][1] = 0.5403120875137825;
      c[0][2] = 0.06440920730334501;   c[1][2] = 0.06440920730334501;
      c[0][3] = 0.2881594349535524;    c[1][3] = 0.2881594349535524;
      c[0][4] = 0.04948438583750772;   c[1][4] = 0.04948438583750772;
      c[0][5] = 0.07152523367501308;   c[1][5] = 0.07152523367501308;
    }
  else if ( alpha_s >= 0.21 && alpha_s < 0.24 )
    {
      c[0][0] = 0.15243607717720586;   c[1][0] = 0.21661000995329158;
      c[0][1] = 0.5403120875137825;    c[1][1] = 1.4087570376612657;
      c[0][2] = 0.06440920730334501;   c[1][2] = 0.2713885880193171;
      c[0][3] = 0.2881594349535524;    c[1][3] = 0.48681971936565244;
      c[0][4] = 0.04948438583750772;   c[1][4] = 0.09567346780679847;
      c[0][5] = 0.07152523367501308;   c[1][5] = 0.12780677622585393;
    }
  else if ( alpha_s >= 0.24 && alpha_s < 0.27 )
    {
      c[0][0] = 0.21661000995329158;   c[1][0] = 0.2501007467879627;
      c[0][1] = 1.4087570376612657;    c[1][1] = 1.8034683081244214;
      c[0][2] = 0.2713885880193171;    c[1][2] = 0.228092470920281;
      c[0][3] = 0.48681971936565244;   c[1][3] = 0.6841577896561725;
      c[0][4] = 0.09567346780679847;   c[1][4] = 0.15430793601338547;
      c[0][5] = 0.12780677622585393;   c[1][5] = 0.1648297331159989;
    }
  else if ( alpha_s >= 0.27 && alpha_s < 0.3 )
    {
      c[0][0] = 0.2501007467879627;    c[1][0] = 0.28440720063047276;
      c[0][1] = 1.8034683081244214;    c[1][1] = 2.0448244620634055;
      c[0][2] = 0.228092470920281;     c[1][2] = -0.018574547528236382;
      c[0][3] = 0.6841577896561725;    c[1][3] = 0.9863974758613413;
      c[0][4] = 0.15430793601338547;   c[1][4] = 0.2503738253300167;
      c[0][5] = 0.1648297331159989;    c[1][5] = 0.2090067594645225;
    }
  else if ( alpha_s >= 0.3 && alpha_s < 0.33 )
    {
      c[0][0] = 0.28440720063047276;   c[1][0] = 0.31945943548344036;
      c[0][1] = 2.0448244620634055;    c[1][1] = 2.0482495934952256;
      c[0][2] = -0.018574547528236382; c[1][2] = -0.5350999123662686;
      c[0][3] = 0.9863974758613413;    c[1][3] = 1.4169725257394696;
      c[0][4] = 0.2503738253300167;    c[1][4] = 0.3918202096574105;
      c[0][5] = 0.2090067594645225;    c[1][5] = 0.26103455441873036;
    }
  else if ( alpha_s >= 0.33 && alpha_s < 0.36 )
    {
      c[0][0] = 0.31945943548344036;   c[1][0] = 0.35519799231686516;
      c[0][1] = 2.0482495934952256;    c[1][1] = 1.7485135425544152;
      c[0][2] = -0.5350999123662686;   c[1][2] = -1.3692232011881413;
      c[0][3] = 1.4169725257394696;    c[1][3] = 1.9906086576701993;
      c[0][4] = 0.3918202096574105;    c[1][4] = 0.5832315715098879;
      c[0][5] = 0.26103455441873036;   c[1][5] = 0.32124694953933486;
    }
  else if ( alpha_s >= 0.36 && alpha_s < 0.39 )
    {
      c[0][0] = 0.35519799231686516;   c[1][0] = 0.39157507493019383;
      c[0][1] = 1.7485135425544152;    c[1][1] = 1.0778995684787331;
      c[0][2] = -1.3692232011881413;   c[1][2] = -2.5738838613236457;
      c[0][3] = 1.9906086576701993;    c[1][3] = 2.727543221296746;
      c[0][4] = 0.5832315715098879;    c[1][4] = 0.8323699786704292;
      c[0][5] = 0.32124694953933486;   c[1][5] = 0.3905055907877247;
    }
  else if ( alpha_s >= 0.39 && alpha_s <= 0.42 )
    {
      c[0][0] = 0.39157507493019383;   c[1][0] = 0.4285382777192131;
      c[0][1] = 1.0778995684787331;    c[1][1] = 0.05505396151716547;
      c[0][2] = -2.5738838613236457;   c[1][2] = -4.113979132685303;
      c[0][3] = 2.727543221296746;     c[1][3] = 3.5992808060371506;
      c[0][4] = 0.8323699786704292;    c[1][4] = 1.1252568207814462;
      c[0][5] = 0.3905055907877247;    c[1][5] = 0.4667953957378259;
    }

  rateLower = T*(c[0][0]             + c[0][1]/pow(u, 4.) -
		 c[0][2]/pow(u, 3.)  - c[0][3]/pow(u, 2.) +
		 c[0][4]/pow(u, 1.5) - c[0][5]/u);
  rateUpper = T*(c[1][0]             + c[1][1]/pow(u, 4.) -
		 c[1][2]/pow(u, 3.)  - c[1][3]/pow(u, 2.) +
		 c[1][4]/pow(u, 1.5) - c[1][5]/u);

  rate.qq = (1.-alphaFrac)*rateLower+alphaFrac*rateUpper;
  rate.qq *= nf/3.; // adjust number of flavors

  rate.gq = rate.qq*9./4.;

  if ( alpha_s >= 0.15 && alpha_s < 0.18 )
    {
      d[0][0] = 0.6197775378922895;    d[1][0] = 0.7680959463632293;
      d[0][1] = 1.5268694134079064;    d[1][1] = 3.282164035377037;
      d[0][2] = 0.6939337312845367;    d[1][2] = 1.6359849897319092;
      d[0][3] = 0.5967602676773388;    d[1][3] =  0.6770046238563808;
      d[0][4] = 0.17320784052297564;   d[1][4] = 0.22074166337990309;
      d[0][5] = 0.28964614117694565;   d[1][5] = 0.4128184793199476;
    }
  else if ( alpha_s >= 0.18 && alpha_s < 0.21 )
    {
      d[0][0] = 0.7680959463632293;    d[1][0] = 0.9206225398305536;
      d[0][1] = 3.282164035377037;     d[1][1] = 5.690562370150853;
      d[0][2] = 1.6359849897319092;    d[1][2] = 2.8341906487774318;
      d[0][3] = 0.6770046238563808;    d[1][3] = 0.7900156706763937;
      d[0][4] = 0.22074166337990309;   d[1][4] = 0.2995126102416747;
      d[0][5] = 0.4128184793199476;    d[1][5] = 0.5598645426609049;
    }
  else if ( alpha_s >= 0.21 && alpha_s < 0.24 )
    {
      d[0][0] = 0.9206225398305536;    d[1][0] = 1.0767954081327265;
      d[0][1] = 5.690562370150853;     d[1][1] = 8.378841394880034;
      d[0][2] = 2.8341906487774318;    d[1][2] = 3.9338968631891396;
      d[0][3] = 0.7900156706763937;    d[1][3] = 1.0874771229885156;
      d[0][4] = 0.2995126102416747;    d[1][4] = 0.46570985770548107;
      d[0][5] = 0.5598645426609049;    d[1][5] = 0.7360069767362173;
    }
  else if ( alpha_s >= 0.24 && alpha_s < 0.27 )
    {
      d[0][0] = 1.0767954081327265;    d[1][0] = 1.2361819653856791;
      d[0][1] = 8.378841394880034;     d[1][1] = 10.877148035367144;
      d[0][2] = 3.9338968631891396;    d[1][2] = 4.526191560392149;
      d[0][3] = 1.0874771229885156;    d[1][3] = 1.731930015138816;
      d[0][4] = 0.46570985770548107;   d[1][4] = 0.7769917594310469;
      d[0][5] = 0.7360069767362173;    d[1][5] = 0.9463662091275489;
    }
  else if ( alpha_s >= 0.27 && alpha_s < 0.3 )
    {
      d[0][0] = 1.2361819653856791;    d[1][0] = 1.3984393292278847;
      d[0][1] = 10.877148035367144;    d[1][1] = 12.72181515837248;
      d[0][2] = 4.526191560392149;     d[1][2] = 4.227297031355039;
      d[0][3] = 1.731930015138816;     d[1][3] = 2.868526983329731;
      d[0][4] = 0.7769917594310469;    d[1][4] = 1.2836917844304823;
      d[0][5] = 0.9463662091275489;    d[1][5] = 1.1953148369630755;
    }
  else if ( alpha_s >= 0.3 && alpha_s < 0.33 )
    {
      d[0][0] = 1.3984393292278847;    d[1][0] = 1.5632880021613935;
      d[0][1] = 12.72181515837248;     d[1][1] = 13.502896915302873;
      d[0][2] = 4.227297031355039;     d[1][2] = 2.7113406243010467;
      d[0][3] = 2.868526983329731;     d[1][3] = 4.615035662049938;
      d[0][4] = 1.2836917844304823;    d[1][4] = 2.0259357821768784;
      d[0][5] = 1.1953148369630755;    d[1][5] = 1.486253368704046;
    }
  else if ( alpha_s >= 0.33 && alpha_s < 0.36 )
    {
      d[0][0] = 1.5632880021613935;    d[1][0] = 1.730492163581557;
      d[0][1] = 13.502896915302873;    d[1][1] = 12.913294655478987;
      d[0][2] = 2.7113406243010467;    d[1][2] = -0.2477159937428581;
      d[0][3] = 4.615035662049938;     d[1][3] = 7.042004003229154;
      d[0][4] = 2.0259357821768784;    d[1][4] = 3.0253452576771465;
      d[0][5] = 1.486253368704046;     d[1][5] = 1.8205651561017433;
    }
  else if ( alpha_s >= 0.36 && alpha_s < 0.39 )
    {
      d[0][0] = 1.730492163581557;     d[1][0] = 1.8998560359992867;
      d[0][1] = 12.913294655478987;    d[1][1] = 10.708892844334745;
      d[0][2] = -0.2477159937428581;   d[1][2] = -4.823210983922782;
      d[0][3] = 7.042004003229154;     d[1][3] = 10.202109059054063;
      d[0][4] = 3.0253452576771465;    d[1][4] = 4.298747764427364;
      d[0][5] = 1.8205651561017433;    d[1][5] = 2.199497022778097;
    }
  else if ( alpha_s >= 0.39 && alpha_s <= 0.42 )
    {
      d[0][0] = 1.8998560359992867;    d[1][0] = 2.071204284004704;
      d[0][1] = 10.708892844334745;    d[1][1] = 6.741738604119316;
      d[0][2] = -4.823210983922782;    d[1][2] = -11.099716230158746;
      d[0][3] = 10.202109059054063;    d[1][3] = 14.106488110189458;
      d[0][4] = 4.298747764427364;     d[1][4] = 5.846203546614067;
      d[0][5] = 2.199497022778097;     d[1][5] = 2.62230136903594;
    }

  rateLower = T*(d[0][0]             + d[0][1]/pow(u, 4.) -
		 d[0][2]/pow(u, 3.)  - d[0][3]/pow(u, 2.) +
		 d[0][4]/pow(u, 1.5) - d[0][5]/u);
  rateUpper = T*(d[1][0]             + d[1][1]/pow(u, 4.) -
		 d[1][2]/pow(u, 3.)  - d[1][3]/pow(u, 2.) +
		 d[1][4]/pow(u, 1.5) - d[1][5]/u);
      
  rate.qg = (1.-alphaFrac)*rateLower+alphaFrac*rateUpper;
  rate.qg /= 3.; // historic reasons
  
  rate.gg = rate.qg*9./4.;

  return rate;
}

RateElastic Martini::getRateElasNeg(double u, double T)
{
  RateElastic rate;

  double alpha0 = 0.15;
  double deltaAlpha = 0.03;
  double iAlpha = floor((alpha_s-alpha0)/deltaAlpha+0.001);
  double alphaFrac = (alpha_s-alpha0)/deltaAlpha - iAlpha;
  double rateLower, rateUpper;
  double c[2][6], d[2][6];

  for (int i=0; i<2; i++)
    for (int j=0; j<6; j++)
      {
	c[i][j] = 0.;
	d[i][j] = 0.;
      }

  if ( alpha_s >= 0.15 && alpha_s < 0.18 )
    {
      c[0][0] = 0.059730780828164666;  c[1][0] = 0.07216040121873951;
      c[0][1] = 0.3631534872548789;    c[1][1] = 0.5471138972952214;
      c[0][2] = 0.39845046966687;      c[1][2] = 0.5792306465939813;
      c[0][3] = 0.11946615972422633;   c[1][3] = 0.1723078888161528;
      c[0][4] = 0.03559276204445307;   c[1][4] = 0.05120408756812135;
      c[0][5] = 0.00291041965645416;   c[1][5] = 0.0041777787108426695;
    }
  else if ( alpha_s >= 0.18 && alpha_s < 0.21 )
    {
      c[0][0] = 0.07216040121873951;   c[1][0] = 0.0846236909779996;
      c[0][1] = 0.5471138972952214;    c[1][1] = 0.7725791286875564;
      c[0][2] = 0.5792306465939813;    c[1][2] = 0.7931123494736929;
      c[0][3] = 0.1723078888161528;    c[1][3] = 0.23406373724706608;
      c[0][4] = 0.05120408756812135;   c[1][4] = 0.06935459958589639;
      c[0][5] = 0.0041777787108426695; c[1][5] = 0.005644055718614478;
    }
  else if ( alpha_s >= 0.21 && alpha_s < 0.24 )
    {
      c[0][0] = 0.0846236909779996;    c[1][0] = 0.09711346786308672;
      c[0][1] = 0.7725791286875564;    c[1][1] = 1.0370070423372528;
      c[0][2] = 0.7931123494736929;    c[1][2] = 1.036935526583188;
      c[0][3] = 0.23406373724706608;   c[1][3] = 0.3034025403259155;
      c[0][4] = 0.06935459958589639;   c[1][4] = 0.08957509599955729;
      c[0][5] = 0.005644055718614478;  c[1][5] = 0.007264393465593115;
    }
  else if ( alpha_s >= 0.24 && alpha_s < 0.27 )
    {
      c[0][0] = 0.09711346786308672;   c[1][0] = 0.10962479860948156;
      c[0][1] = 1.0370070423372528;    c[1][1] = 1.3372010137066646;
      c[0][2] = 1.036935526583188;     c[1][2] = 1.307456863105879;
      c[0][3] = 0.3034025403259155;    c[1][3] = 0.37910328734850873;
      c[0][4] = 0.08957509599955729;   c[1][4] = 0.111456899829735;
      c[0][5] = 0.007264393465593115;  c[1][5] = 0.009000895144744121;
    }
  else if ( alpha_s >= 0.27 && alpha_s < 0.3 )
    {
      c[0][0] = 0.10962479860948156;   c[1][0] = 0.1221541053951596;
      c[0][1] = 1.3372010137066646;    c[1][1] = 1.6686065099273535;
      c[0][2] = 1.307456863105879;     c[1][2] = 1.600404353394210;
      c[0][3] = 0.37910328734850873;   c[1][3] = 0.4594932213772782;
      c[0][4] = 0.111456899829735;     c[1][4] = 0.13442407314203592;
      c[0][5] = 0.009000895144744121;  c[1][5] = 0.010800449048880756;
    }
  else if ( alpha_s >= 0.3 && alpha_s < 0.33 )
    {
      c[0][0] = 0.1221541053951596;    c[1][0] = 0.13469861652518803;
      c[0][1] = 1.6686065099273535;    c[1][1] = 2.0276317271182074;
      c[0][2] = 1.600404353394210;     c[1][2] = 1.912613330851788;
      c[0][3] = 0.4594932213772782;    c[1][3] = 0.5434449889160747;
      c[0][4] = 0.13442407314203592;   c[1][4] = 0.15810564016236883;
      c[0][5] = 0.010800449048880756;  c[1][5] = 0.012629305933671075;
    }
  else if ( alpha_s >= 0.33 && alpha_s < 0.36 )
    {
      c[0][0] = 0.13469861652518803;   c[1][0] = 0.14725614907227047;
      c[0][1] = 2.0276317271182074;    c[1][1] = 2.4109122726272654;
      c[0][2] = 1.912613330851788;     c[1][2] = 2.241198157777867;
      c[0][3] = 0.5434449889160747;    c[1][3] = 0.6299396046048817;
      c[0][4] = 0.15810564016236883;   c[1][4] = 0.18216575652552597;
      c[0][5] = 0.012629305933671075;  c[1][5] = 0.014456750325370632;
    }
  else if ( alpha_s >= 0.36 && alpha_s < 0.39 )
    {
      c[0][0] = 0.14725614907227047;   c[1][0] = 0.15982489441001274;
      c[0][1] = 2.4109122726272654;    c[1][1] = 2.815254291049982;
      c[0][2] = 2.241198157777867;     c[1][2] = 2.583462624103292;
      c[0][3] = 0.6299396046048817;    c[1][3] = 0.7180274724508857;
      c[0][4] = 0.18216575652552597;   c[1][4] = 0.20629432847931367;
      c[0][5] = 0.014456750325370632;  c[1][5] = 0.01625568033747704;
    }
  else if ( alpha_s >= 0.39 && alpha_s <= 0.42 )
    {
      c[0][0] = 0.15982489441001274;   c[1][0] = 0.17240331582158486;
      c[0][1] = 2.815254291049982;     c[1][1] = 3.238290376079149;
      c[0][2] = 2.583462624103292;     c[1][2] = 2.9374985881586273;
      c[0][3] = 0.7180274724508857;    c[1][3] = 0.8071008047950518;
      c[0][4] = 0.20629432847931367;   c[1][4] = 0.23030341585944009;
      c[0][5] = 0.01625568033747704;   c[1][5] = 0.018010096397556033;
    }

  rateLower = T*(c[0][0]             + c[0][1]/pow(u, 4.) -
		 c[0][2]/pow(u, 3.)  - c[0][3]/pow(u, 2.) +
		 c[0][4]/pow(u, 1.5) - c[0][5]/u);
  rateUpper = T*(c[1][0]             + c[1][1]/pow(u, 4.) -
		 c[1][2]/pow(u, 3.)  - c[1][3]/pow(u, 2.) +
		 c[1][4]/pow(u, 1.5) - c[1][5]/u);

  rate.qq = (1.-alphaFrac)*rateLower+alphaFrac*rateUpper;
  rate.qq *= nf/3.; // adjust number of flavors

  rate.gq = rate.qq*9./4.;

  if (alpha_s >= 0.15 && alpha_s < 0.18)
    {
      d[0][0] = 0.3166913701414167;    d[1][0] = 0.3804527486448292;
      d[0][1] = 1.0992070651449564;    d[1][1] = 1.7114836115114735;
      d[0][2] = 1.4232219292986843;    d[1][2] = 2.093540209692791;
      d[0][3] = 0.4655268754156213;    d[1][3] = 0.6787666526042345;
      d[0][4] = 0.1444497238817506;    d[1][4] = 0.21014340589291994;
      d[0][5] = 0.012281442189758532;  d[1][5] = 0.017832857383112792;
    }
  else if (alpha_s >= 0.18 && alpha_s < 0.21)
    {
      d[0][0] = 0.3804527486448292;    d[1][0] = 0.44393432393104637;
      d[0][1] = 1.7114836115114735;    d[1][1] = 2.483663499207573;
      d[0][2] = 2.093540209692791;     d[1][2] = 2.8979112438999044;
      d[0][3] = 0.6787666526042345;    d[1][3] = 0.9316968286688833;
      d[0][4] = 0.21014340589291994;   d[1][4] = 0.28780901378857465;
      d[0][5] = 0.017832857383112792;  d[1][5] = 0.02438874287373154;
    }
  else if (alpha_s >= 0.21 && alpha_s < 0.24)
    {
      d[0][0] = 0.44393432393104637;   d[1][0] = 0.5071424487228405;
      d[0][1] = 2.483663499207573;     d[1][1] = 3.4070556051784515;
      d[0][2] = 2.8979112438999044;    d[1][2] = 3.824491419496227;
      d[0][3] = 0.9316968286688833;    d[1][3] = 1.2191109771387096;
      d[0][4] = 0.28780901378857465;   d[1][4] = 0.3755459972857442;
      d[0][5] = 0.02438874287373154;   d[1][5] = 0.03174924882247299;
    }
  else if (alpha_s >= 0.24 && alpha_s < 0.27)
    {
      d[0][0] = 0.5071424487228405;    d[1][0] = 0.5700856365203443;
      d[0][1] = 3.4070556051784515;    d[1][1] = 4.466964606692036;
      d[0][2] = 3.824491419496227;     d[1][2] = 4.857999356928031;
      d[0][3] = 1.2191109771387096;    d[1][3] = 1.5348360053714125;
      d[0][4] = 0.3755459972857442;    d[1][4] = 0.471215528026891;
      d[0][5] = 0.03174924882247299;   d[1][5] = 0.03971601962636114;
    }
  else if (alpha_s >= 0.27 && alpha_s < 0.3)
    {
      d[0][0] = 0.5700856365203443;    d[1][0] = 0.6327732610959403;
      d[0][1] = 4.466964606692036;     d[1][1] = 5.646624908846933;
      d[0][2] = 4.857999356928031;     d[1][2] = 5.982691423451806;
      d[0][3] = 1.5348360053714125;    d[1][3] = 1.8728243844356356;
      d[0][4] = 0.471215528026891;     d[1][4] = 0.572761497659723;
      d[0][5] = 0.03971601962636114;   d[1][5] = 0.04809998538877525;
    }
  else if (alpha_s >= 0.3 && alpha_s < 0.33)
    {
      d[0][0] = 0.6327732610959403;    d[1][0] = 0.6952147319486842;
      d[0][1] = 5.646624908846933;     d[1][1] = 6.931552369487635;
      d[0][2] = 5.982691423451806;     d[1][2] = 7.185588273540373;
      d[0][3] = 1.8728243844356356;    d[1][3] = 2.228328283532209;
      d[0][4] = 0.572761497659723;     d[1][4] = 0.6786029643259804;
      d[0][5] = 0.04809998538877525;   d[1][5] = 0.056755908207122875;
    }
  else if (alpha_s >= 0.33 && alpha_s < 0.36)
    {
      d[0][0] = 0.6952147319486842;    d[1][0] = 0.7574189285141091;
      d[0][1] = 6.931552369487635;     d[1][1] = 8.307255807497631;
      d[0][2] = 7.185588273540373;     d[1][2] = 8.454112812202247;
      d[0][3] = 2.228328283532209;     d[1][3] = 2.596781386863294;
      d[0][4] = 0.6786029643259804;    d[1][4] = 0.7872276571283385;
      d[0][5] = 0.056755908207122875;  d[1][5] = 0.06554867983133447;
    }
  else if (alpha_s >= 0.36 && alpha_s < 0.39)
    {
      d[0][0] = 0.7574189285141091;    d[1][0] = 0.8193940883937045;
      d[0][1] = 8.307255807497631;     d[1][1] = 9.761691032241623;
      d[0][2] = 8.454112812202247;     d[1][2] = 9.777948193339808;
      d[0][3] = 2.596781386863294;     d[1][3] = 2.9744411293541457;
      d[0][4] = 0.7872276571283385;    d[1][4] = 0.8973688582323887;
      d[0][5] = 0.06554867983133447;   d[1][5] = 0.07435862848596686;
    }
  else if (alpha_s >= 0.39 && alpha_s <= 0.42)
    {
      d[0][0] = 0.8193940883937045;    d[1][0] = 0.8811479514201789;
      d[0][1] = 9.761691032241623;     d[1][1] = 11.286034194965852;
      d[0][2] = 9.777948193339808;     d[1][2] = 11.15001447311135;
      d[0][3] = 2.9744411293541457;    d[1][3] = 3.3591358778545803;
      d[0][4] = 0.8973688582323887;    d[1][4] = 1.0083901554550654;
      d[0][5] = 0.07435862848596686;   d[1][5] = 0.08313659597360733;
    }

  rateLower = T*(d[0][0]             + d[0][1]/pow(u, 4.) -
		 d[0][2]/pow(u, 3.)  - d[0][3]/pow(u, 2.) +
		 d[0][4]/pow(u, 1.5) - d[0][5]/u);
  rateUpper = T*(d[1][0]             + d[1][1]/pow(u, 4.) -
		 d[1][2]/pow(u, 3.)  - d[1][3]/pow(u, 2.) +
		 d[1][4]/pow(u, 1.5) - d[1][5]/u);
      
  rate.qg = (1.-alphaFrac)*rateLower+alphaFrac*rateUpper;
  rate.qg /= 3.; // historic reasons
  
  rate.gg = rate.qg*9./4.;

  return rate;
}

RateConversion Martini::getRateConv(double pRest, double T)
{
  RateConversion rate;

  rate.qg     = 4./3.*2.*M_PI*alpha_s*alpha_s*T*T/(3.*pRest)
	       *(0.5*log(pRest*T/((1./6.)*pow(g*T, 2.)))-0.36149);
  rate.gq     = nf*3./8.*4./3.*2.*M_PI*alpha_s*alpha_s*T*T/(3.*pRest)
	       *(0.5*log(pRest*T/((1./6.)*pow(g*T, 2.)))-0.36149);
  rate.qgamma = 2.*M_PI*alpha_s*alpha_em*T*T/(3.*pRest)
	       *(0.5*log(pRest*T/((1./6.)*pow(g*T, 2.)))-0.36149);

  return rate;
}

double Martini::getEnergyTransfer(double pRest, double T, int process)
{
  double u = pRest/T; // making arguments in log to be dimensionless

  double omega = 0.;  // momentum of radiated gluon (dimentionless)
  double y;           // omega candidate
  double x;           // random number, uniform on [0,1]
  double randA;       // uniform random number on [0, Area under the envelop function]
  double fy;          // total area under the envelop function
  double fyAct;       // actual rate

  RateElastic Pos, Neg;
  Pos = getRateElasPos(u, T);
  Neg = getRateElasNeg(u, T);
  
  // this switch will hold the decision whether k is positive or negative:
  // 0 : negative, 1 : positive
  int posNegSwitch = 1; 
  
  /* process == 5 : qq -> qq (first q is hard, second q is soft) 
     process == 6 : qg -> qg 
     process == 7 : gq -> gq 
     process == 8 : gg -> gg  */

  if (process == 5 || process == 7) // for qq or gq
    {
      // decide whether k shall be positive or negative 
      // if x (uniform on [0,1]) < area(k<0)/area(all k) then k < 0
      if (process == 5)
	{
	  if (ZeroOneDistribution(*GetMt19937Generator()) < Neg.qq/(Neg.qq+Pos.qq))
	    posNegSwitch = 0;
	}
      else
	{
	  if (ZeroOneDistribution(*GetMt19937Generator()) < Neg.gq/(Neg.gq+Pos.gq))
	    posNegSwitch = 0;
	}

      if (posNegSwitch == 1) // if omega > 0
	{
	  do
	    {
	      randA = ZeroOneDistribution(*GetMt19937Generator())
		     *areaOmega(u, posNegSwitch, process);
	      y = exp((-1.41428*pow(10., 9.)*alpha_s
		      - 8.08158*pow(10., 8.)*alpha_s*alpha_s + 2.02327*pow(10., 9.)*randA)
		      /(alpha_s*(4.72097*pow(10., 8.) + 2.6977*pow(10., 8.)*alpha_s)));

	      fy = alpha_s/0.15*(0.035 + alpha_s*0.02)/sqrt(y*y);
	      fyAct = functionOmega(u, y, process);

	      x = ZeroOneDistribution(*GetMt19937Generator());

	    } while (x > fyAct/fy); 
	  // reject if x is larger than the ratio fyAct/fy
	  omega = y;
	}
      else // if omega < 0
	{
	  do
	    {
	      randA = ZeroOneDistribution(*GetMt19937Generator())
		     *areaOmega(-0.05, posNegSwitch, process);
	      y = -12.*exp(-30.*randA/(alpha_s*(7.+ 4.*alpha_s)));

	      fy = alpha_s/0.15*(0.035 + alpha_s*0.02)/sqrt(y*y);
	      fyAct = functionOmega(u, y, process);

	      x = ZeroOneDistribution(*GetMt19937Generator());

	    } while (x > fyAct/fy); 
	  // reject if x is larger than the ratio fyAct/fy
	  omega = y;
	}
    }
  else if (process == 6 || process == 8) // for qg or gg
    {
      // decide whether k shall be positive or negative 
      // if x (uniform on [0,1]) < area(k<0)/area(all k) then k < 0
      if (process == 6)
	{
	  if (ZeroOneDistribution(*GetMt19937Generator()) < Neg.qg/(Neg.qg+Pos.qg))
	    posNegSwitch = 0;
	}
      else
	{
	  if (ZeroOneDistribution(*GetMt19937Generator()) < Neg.gg/(Neg.gg+Pos.gg))
	    posNegSwitch = 0;
	}

      if( posNegSwitch == 1 ) // if omega > 0
	{
	  do
	    {
	      randA = ZeroOneDistribution(*GetMt19937Generator())
		     *areaOmega(u, posNegSwitch, process);
	      y = exp((-2.32591*pow(10.,17.)*alpha_s
		      - 1.32909*pow(10.,17.)*alpha_s*alpha_s + 2.2183*pow(10.,17.)*randA)
		      /(alpha_s*(7.76406*pow(10.,16.) + 4.43661*pow(10.,16.)*alpha_s)));

	      fy = 1.5*alpha_s/0.15*(0.035 + alpha_s*0.02)/sqrt(y*y);
	      fyAct = functionOmega(u, y, process);

	      x = ZeroOneDistribution(*GetMt19937Generator());

	    } while (x > fyAct/fy); 
	  // reject if x is larger than the ratio fyAct/fy
	  omega = y;
	}
      else // if omega < 0
	{
	  do
	    {
	      randA = ZeroOneDistribution(*GetMt19937Generator())
		     *areaOmega(-0.05, posNegSwitch, process);
	      y = -12.*exp(-2.81475*pow(10.,15.)*randA
			 /(alpha_s*(9.85162*pow(10.,14.) + 5.6295*pow(10.,14.)*alpha_s)));

	      fy = 1.5*alpha_s/0.15*(0.035 + alpha_s*0.02)/sqrt(y*y);
	      fyAct = functionOmega(u, y, process);

	      x = ZeroOneDistribution(*GetMt19937Generator());

	    } while (x > fyAct/fy); 
	  // reject if x is larger than the ratio fyAct/fy
	  omega = y;
	}
    }
  else
    {
      JSWARN << "Invalid process number (" << process << ")";

      omega = 0.;
    }

  return omega*T;  // omega*T is in [GeV]
}

double Martini::getMomentumTransfer(double pRest, double omega, double T, int process)
{
  double u = pRest/T; // making arguments in log to be dimensionless
  omega /= T;

  double q = 0.;      // momentum of radiated gluon (dimentionless)
  double y;           // q candidate
  double x;           // random number, uniform on [0,1]
  double randA;       // uniform random number on [0, Area under the envelop function]
  double fy;          // total area under the envelop function
  double fyAct;       // actual rate

  double A, B;

  /* process == 5 : qq -> qq (first q is hard, second q is soft)
     process == 6 : qg -> qg 
     process == 7 : gq -> gq 
     process == 8 : gg -> gg  */

  // small omega using the rejection method
  if (omega < 10. && omega > -3.)
    {
      if (process == 5 || process == 7) // for qq or gq
	{
	  A = (0.7+alpha_s)*0.0014
	     *(1000.+40./sqrt(omega*omega)+10.5*pow(omega, 4.))*alpha_s;
	  B = 2.*sqrt(omega*omega)+0.01;
	}
      else if (process == 6 || process == 8) // for qg or gg
	{
	  A = (0.7+alpha_s)*0.0022
	     *(1000.+40./sqrt(omega*omega)+10.5* pow(omega, 4.))*alpha_s;
	  B = 2.*sqrt(omega*omega)+0.002;
	}
      else
	{
	  JSWARN << "Invalid process number (" << process << ")";

	  A = 0.;
	  B = 0.;
	}

      do
	{
	  randA = ZeroOneDistribution(*GetMt19937Generator())*areaQ(u, omega, process);
	  y = pow(B, 0.25)*sqrt(tan((2.*sqrt(B)*randA+A*atan(omega*omega/sqrt(B)))/A));

	  fy = A*y/(pow(y, 4.)+B);
	  fyAct = functionQ(u, omega, y, process);

	  x = ZeroOneDistribution(*GetMt19937Generator());

	} while (x > fyAct/fy); 
      // reject if x is larger than the ratio fyAct/fy
      q = y;
    }
  // large omega using the Metropolis method
  else
    {
      double g = 0, g_new = 0;
      double ratio;
      double y_new;
    
      // the ranges in which the variables u and phi need to be sampled
      const double y_min = sqrt(omega*omega);
      const double y_max = u;

      int count = 0;
      // randomly select initial values of q=y, such that
      do
	{
	  y = y_min+ZeroOneDistribution(*GetMt19937Generator())*(y_max-y_min);
	  g = functionQ(u, omega, y, process);

	  // if no y having non-zero g is found, give up here.
	  if(count > 100) return 0.;
	  count++;
	} while (g == 0.);
    
      // number of steps in the Markov chain
      const int n_steps = 500;
    
      for(int i=0; i<n_steps; i++)
	{
	  do
	    {
	      y_new = y_min+ZeroOneDistribution(*GetMt19937Generator())*(y_max-y_min);
	    }
	  while (y_new < y_min || y_new > y_max);
	  // check that the new value is in range

	  // calculate the function at the proposed point
	  g_new = functionQ(u, omega, y_new, process);
	  // ratio of g(y_new)/g(y)
	  ratio = g_new/g;

	  // accept if probability g(y_new)/g(y) is larger than randon number
	  if (ZeroOneDistribution(*GetMt19937Generator()) < ratio)
	    {
	      y = y_new;
	      g = g_new;
	    }
	}
      q = y;
    }
  
  return q*T;  // q*T is in [GeV]
}

// calculates the area under the envelope function when using the rejection method
// integrals had been solved analytically before
double Martini::areaOmega(double u, int posNegSwitch, int process)
{
  if (process == 5 || process == 7)
    {
      if (posNegSwitch == 1)
	return (0.0333333*alpha_s*(7.+ 4.*alpha_s)*(2.99573+log(u)));
      else 
	return (-0.133333*alpha_s*(1.75+alpha_s)*log(-0.0833333*u));
    }
  else if (process == 6 || process == 8)
    {
      if (posNegSwitch == 1)
	return (0.05*alpha_s*(7.+ 4.*alpha_s)*(2.99573+log(u)));
      else 
	return (-0.2*alpha_s*(1.75+alpha_s)*log(-0.0833333*u));
    }
  else
    {
      JSWARN << "Invalid process number (" << process << ")";
    }

  return 0.;
}

double Martini::areaQ(double u, double omega, int process)
{
  double A, B;
  double areaQ;

  if (process == 5 || process == 7)
    {
      A = (0.7+alpha_s-0.3)*0.0014*alpha_s
	 *(1000.+40./sqrt(omega*omega)+10.5*pow(omega, 4.));
      B = 2.*sqrt(omega*omega)+0.01;
    }
  else if (process == 6 || process == 8)
    {
      A = (0.7+alpha_s-0.3)*0.0022*alpha_s
	 *(1000.+40./sqrt(omega*omega)+10.5*pow(omega, 4.));
      B = 2.*sqrt(omega*omega)+0.002;
    }
  else
    {
      JSWARN << "Invalid process number (" << process << ")";

      A = 0.;
      B = 0.;
    }
  areaQ = (0.5*A*(atan(u*u/sqrt(B))-atan(omega*omega/sqrt(B))))/sqrt(B);

  return areaQ;
}

FourVector Martini::getNewMomentumElas(FourVector pVec, double omega, double q)
{
  FourVector pVecNew, pVecNewTemp;
  FourVector etVec, qtVec, qlVec;
  FourVector r;
  double qt, ql;
  double cosTheta_pq;
  double pAbs=pVec.t();
  double phi;
  double M[3][3]; //rotation matrix
  double u;
  double xx, yy, zz ,tt;

  if (omega == q)
    {
      xx = pVec.x()*(pAbs-omega)/pAbs;
      yy = pVec.y()*(pAbs-omega)/pAbs;
      zz = pVec.z()*(pAbs-omega)/pAbs;
      tt = pVec.t()*(pAbs-omega)/pAbs;
    
      pVecNew.Set(xx, yy, zz ,tt);
      return pVecNew;
    }

  cosTheta_pq = (-omega*omega+2.*pAbs*omega+q*q)/(2.*pAbs*q);
  qt = q*sqrt(1.-cosTheta_pq*cosTheta_pq);          // transverse momentum transfer
  ql = q*cosTheta_pq;
  
  if (pVec.y()*pVec.y() > pVec.x()*pVec.x())
    {
      xx = 0.;
      yy = -pVec.z();
      zz = pVec.y();
    }
  else
    {
      xx = pVec.z();
      yy = 0.;
      zz = -pVec.x();
    }

  tt = sqrt(xx*xx+yy*yy+zz*zz);
  etVec.Set(xx/tt, yy/tt, zz/tt, 1.);  // normalized to 1

  // the transverse transferred momentum vector
  qtVec.Set(etVec.x()*qt, etVec.y()*qt, etVec.z()*qt, etVec.t()*qt);
  // the longuitudinal transferred momentum vector
  qlVec.Set(pVec.x()/pAbs*ql, pVec.y()/pAbs*ql,
	    pVec.z()/pAbs*ql, pVec.t()/pAbs*ql);

  pVecNewTemp = pVec;
  pVecNewTemp -= qtVec;  // change transverse momentum
  pVecNewTemp -= qlVec;  // change longitudinal momentum


  phi = 2.*M_PI*ZeroOneDistribution(*GetMt19937Generator());
  r.Set(pVec.x()/pVec.t(), pVec.y()/pVec.t(), pVec.z()/pVec.t(), 1.);
  u = 1.-cos(phi);
 
  // define the rotation matrix for rotations around pvecRest
  M[0][0]=r.x()*r.x()*u+cos(phi);
  M[1][0]=r.x()*r.y()*u-r.z()*sin(phi);
  M[2][0]=r.x()*r.z()*u+r.y()*sin(phi);

  M[0][1]=r.y()*r.x()*u+r.z()*sin(phi);
  M[1][1]=r.y()*r.y()*u+cos(phi);
  M[2][1]=r.y()*r.z()*u-r.x()*sin(phi);

  M[0][2]=r.z()*r.x()*u-r.y()*sin(phi);
  M[1][2]=r.z()*r.y()*u+r.x()*sin(phi);
  M[2][2]=r.z()*r.z()*u+cos(phi);

  xx = M[0][0]*pVecNewTemp.x()+M[0][1]*pVecNewTemp.y()+M[0][2]*pVecNewTemp.z();
  yy = M[1][0]*pVecNewTemp.x()+M[1][1]*pVecNewTemp.y()+M[1][2]*pVecNewTemp.z();
  zz = M[2][0]*pVecNewTemp.x()+M[2][1]*pVecNewTemp.y()+M[2][2]*pVecNewTemp.z();
  tt = sqrt(xx*xx+yy*yy+zz*zz);

  pVecNew.Set(xx, yy, zz, tt);
  return pVecNew;
}

// Reads in the binary stored file of dGamma values
void Martini::readRadiativeRate(Gamma_info *dat, dGammas *Gam)
{
  FILE *rfile;
  string filename;
  filename = PathToTables+"radgamma";

  JSINFO << "Reading rates of inelastic collisions from file ";
  JSINFO << filename.c_str() << " ... ";
  size_t bytes_read;

  rfile = fopen(filename.c_str(), "rb"); 
  bytes_read = fread((char *)(&dat->ddf), sizeof(double), 1, rfile);
  bytes_read = fread((char *)(&dat->dda), sizeof(double), 1, rfile);
  bytes_read = fread((char *)(&dat->dcf), sizeof(double), 1, rfile);
  bytes_read = fread((char *)(&dat->dca), sizeof(double), 1, rfile);
  bytes_read = fread((char *)(&dat->Nc), sizeof(int), 1, rfile);
  bytes_read = fread((char *)(&dat->Nf), sizeof(int), 1, rfile);
  bytes_read = fread((char *)(&dat->BetheHeitler),sizeof(int) , 1, rfile);
  bytes_read = fread((char *)(&dat->BDMPS), sizeof(int), 1, rfile);
  bytes_read = fread((char *)(&dat->include_gluons), sizeof(int), 1, rfile);
  bytes_read = fread((char *)Gam->qqg, sizeof(double), NP*NK, rfile);
  bytes_read = fread((char *)Gam->tau_qqg, sizeof(double), NP*NK, rfile);
  bytes_read = fread((char *)Gam->gqq, sizeof(double), NP*NK , rfile);
  bytes_read = fread((char *)Gam->tau_gqq, sizeof(double), NP*NK, rfile);
  bytes_read = fread((char *)Gam->ggg, sizeof(double), NP*NK, rfile);
  bytes_read = fread((char *)Gam->tau_ggg, sizeof(double), NP*NK, rfile);
  bytes_read = fread((char *)Gam->qqgamma, sizeof(double), NP*NK, rfile);
  bytes_read = fread((char *)Gam->tau_qqgamma, sizeof(double), NP*NK, rfile);
  fclose (rfile);

  dat->Nf = nf;
  dat->dp = 0.05;
  dat->p_max = 20;
  dat->p_min = 0;                 // exp(LogEmin) set to zero because my array starts at zero!...
  dat->n_p = static_cast<int>(1.001+dat->p_max/dat->dp);    // np = int(0.4 + 121 / 0.5) = 241
  dat->p_max = dat->dp*dat->n_p;                            // p_max = 0.5*241 = 120.5
  dat->n_pmin = static_cast<int>(1.001+dat->p_min/dat->dp); // np_min = int(0.4 + 3.3 / 0.5) = 7
  dat->n_p -= dat->n_pmin-1;                                // n_p = 241 - (7 - 1) = 235
  dat->p_min = dat->dp * dat->n_pmin;                       // p_min = 0.5 * 7 = 3.5
  dat->n_kmin = 1+2*(static_cast<int>(2./dat->dp));
  dat->k_min = -dat->dp*dat->n_kmin;
  dat->n_k = static_cast<int>((8+dat->p_max)/(2*dat->dp));
  dat->k_max = 2*dat->dp*(dat->n_k-1)+dat->k_min;
}

void Martini::readElasticRateOmega()
{
  ifstream fin;
  string filename[2];

  double as, omega;
  double dGamma;
      
  // open files with data to read in:
  filename[0] = PathToTables + "logEnDtrqq";
  filename[1] = PathToTables + "logEnDtrqg";
  
  JSINFO << "Reading rates of elastic collisions from files";
  JSINFO << filename[0];
  JSINFO << filename[1] << " ...";

  fin.open(filename[0].c_str(), ios::in);
  if(!fin)
    {
      JSWARN << "[readElasticRateOmega]: ERROR: Unable to open file " << filename[0];
      throw std::runtime_error("[readElasticRateQ]: ERROR: Unable to open ElasticRateOmega file");
    }

  int ik = 0;
  while (!fin.eof())
    {
      fin >> as;
      fin >> omega;
      fin >> dGamma;
      dGamma_qq->push_back(dGamma);

      ik++;
    }
  fin.close();

  fin.open(filename[1].c_str(), ios::in);
  if(!fin)
    {
      JSWARN << "[readElasticRateOmega]: ERROR: Unable to open file " << filename[1];
      throw std::runtime_error("[readElasticRateQ]: ERROR: Unable to open ElasticRateOmega file");
    }

  ik = 0;
  while (!fin.eof())
    {
      fin >> as;
      fin >> omega;
      fin >> dGamma;
      dGamma_qg->push_back(dGamma);

      ik++;
    }
  fin.close();
}

void Martini::readElasticRateQ()
{
  ifstream fin;
  string filename[2];

  double as, omega, q;
  double dGamma;
      
  // open files with data to read in:
  filename[0] = PathToTables + "logEnDqtrqq";
  filename[1] = PathToTables + "logEnDqtrqg";
  
  JSINFO << "Reading rates of elastic collisions from files";
  JSINFO << filename[0];
  JSINFO << filename[1] << " ...";

  fin.open(filename[0].c_str(), ios::in);
  if(!fin)
    {
      JSWARN << "[readElasticRateQ]: ERROR: Unable to open file " << filename[0];
      throw std::runtime_error("[readElasticRateQ]: ERROR: Unable to open ElasticRateQ file");
    }

  int ik = 0;
  while (!fin.eof())
    {
      fin >> as;
      fin >> omega;
      fin >> q;
      fin >> dGamma;
      dGamma_qq_q->push_back(dGamma);

      ik++;
    }
  fin.close();

  fin.open(filename[1].c_str(),ios::in);
  if(!fin)
    {
      JSWARN << "[readElasticRateQ]: ERROR: Unable to open file " << filename[1];
      throw std::runtime_error("[readElasticRateQ]: ERROR: Unable to open ElasticRateQ file");
    }
  
  ik = 0;
  while (!fin.eof())
    {
      fin >> as;
      fin >> omega;
      fin >> q;
      fin >> dGamma;
      dGamma_qg_q->push_back(dGamma);

      ik++;
    }
  fin.close();
}

double Martini::getRate_qqg(double p, double k)
{
  return use_table(p, k, Gam.qqg, 0);
}

double Martini::getRate_gqq(double p, double k)
{
  if (k < p/2.) return use_table(p, k, Gam.gqq, 1);
  else return 0.;
}

double Martini::getRate_ggg(double p, double k)
{
  if (k < p/2.) return use_table(p, k, Gam.ggg, 2);
  else return 0.;
}

double Martini::getRate_qqgamma(double p, double k)
{
  return use_table(p, k, Gam.qqgamma, 3);
}

double Martini::use_table(double p, double k, double dGamma[NP][NK], int which_kind)
/* Uses the lookup table and simple interpolation to get the value
   of dGamma/dkdx at some value of p,k.
   This works by inverting the relations between (p,k) and (n_p,n_k)
   used in building the table, to find out what continuous values
   of n_p, n_k should be considered; then linearly interpolates.     */
{
  double a, b, result;     // fraction of way from corner of box
  int    n_p, n_k;         // location of corner of box

  // out of range
  if ((p < 4.01) || (p > 46000.) || (k < -12.) || (k > p+12.))
    return 0.;

  if ((which_kind % 3) && (k > p/2))
    k = p - k;  // Take advantage of symmetry in these cases

  a = 24.7743737154026*log(p*0.2493765586034912718l);
  n_p = (int)a;
  a -= n_p;
  if (k < 2.)
    {
      if (k < -1)
	{
	  if (k < -2) b = 60.+5.*k;
	  else b = 70.+10.*k;
	}
      else
	{
	  if (k < 1.) b = 80. + 20.*k;
	  else b = 90.+10.*k;
	}
    }
  else if ( k < p-2. )
    { /* This is that tricky middle ground. */
      b = 190.-10.*log(1.000670700260932956l/ 
		       (0.0003353501304664781l+(k-2.)/(p-4.))-1.);
    }
  else
    {
      if (k < p+1.)
	{
	  if (k < p-1.) b = 290. + 10.*(k-p);
	  else  b = 300. + 20.*(k-p);
	}
      else
	{
	  if (k < p+2.) b = 310. + 10.*(k-p);
	  else b = 320. + 5.*(k-p);
	}
    }

  n_k = (int)b;
  b -= n_k;
  result = (1.-a)*((1.-b)*dGamma[n_p][n_k]+b*dGamma[n_p][n_k+1])
    +a*((1.-b)*dGamma[n_p+1][n_k]+b*dGamma[n_p+1][n_k+1]);

  if (std::abs(k) > 0.001) // Avoid division by 0, should never get asked for
    {
      switch (which_kind)
	{
	case 0:
	  result /= k;
	  if (k < 20.)
	    result /= 1.-exp(-k);
	  if (k > p-20.)
	    result /= 1. + exp(k-p);
	  break;
	case 1:
	  result /= p;
	  if (k < 20.)
	    result /= 1 + exp(-k);
	  if (k > p-20.)
	    result /= 1. + exp(k-p);
	  break;
	case 2:
	  result /= k*(p-k)/p;
	  if (k < 20.)
	    result /= 1.-exp(-k);
	  if (k > p-20.)
	    result /= 1.-exp(k-p);
	  break;
	case 3:
	  result /= k;
	  if (k < 0) result = 0.;
	  if (k > p-20.)
	    result /= 1. + exp(k-p);
	  break;
	}
    }

  return result;
}

double Martini::getElasticRateOmega(double u, double omega, int process)
{
  if ((omega > 0. && omega < u) || (omega < 0. && omega > -u))
    return use_elastic_table_omega(omega, process);

  return 0.;
}

double Martini::getElasticRateQ(double u, double omega, double q, int process)
{
  if (q > std::abs(omega) && ((omega > 0 && omega < u) || (omega < 0 && omega > -u)))
    return use_elastic_table_q(omega, q, process);

  return 0.;
}

double Martini::use_elastic_table_omega(double omega, int which_kind)
/* Uses the lookup table and simple interpolation to get the value
   of dGamma/domega at some value of omega and alpha_s. */
{
  double result;
  double alphaFrac, omegaFrac;
  int iOmega, iAlphas;
  int position, positionAlphaUp, positionOmegaUp, positionAlphaUpOmegaUp;
  double rate, rateAlphaUp, rateOmegaUp, rateAlphaUpOmegaUp;
  double rateOmegaAv, rateAlphaUpOmegaAv;

  if (omega > 0.) iOmega = Nomega/2+floor((log(omega)+5)/omegaStep);
  else iOmega = Nomega/2-ceil((log(-omega)+5)/omegaStep)-1;
  iAlphas = floor((alpha_s-0.15)/alphaStep);

  position = iOmega+Nomega*(iAlphas);
  positionAlphaUp = iOmega+Nomega*(iAlphas+1);
  positionOmegaUp = iOmega+1+Nomega*(iAlphas);
  positionAlphaUpOmegaUp = iOmega+1+Nomega*(iAlphas+1);

  alphaFrac = (alpha_s-(floor((alpha_s-alphaMin)/alphaStep)*alphaStep+alphaMin))/alphaStep;
  if (omega > 0.) 
    {
      if (exp(ceil((log(omega)+5)/omegaStep)*omegaStep-5)
	  != exp(floor((log(omega)+5)/omegaStep)*omegaStep-5))
	omegaFrac = (omega - (exp(floor((log(omega)+5)/omegaStep)*omegaStep-5)))
	  /((exp(ceil((log(omega)+5)/omegaStep)*omegaStep-5))
	    - exp(floor((log(omega)+5)/omegaStep)*omegaStep-5));
      else omegaFrac = 0.;
    }
  else
    {
      if (exp(ceil((log(-omega)+5)/omegaStep)*omegaStep-5)
	  != exp(floor((log(-omega)+5)/omegaStep)*omegaStep-5))
	omegaFrac = (-omega - (exp(floor((log(-omega)+5)/omegaStep)*omegaStep-5)))
	  /((exp(ceil((log(-omega)+5)/omegaStep)*omegaStep-5))
	    -exp(floor((log(-omega)+5)/omegaStep)*omegaStep-5));
      else omegaFrac = 0.;
    }

  if (which_kind == 5 || which_kind == 7)
    {
      if (position > 0 && iAlphas < Nalphas && iOmega < Nomega) 
	rate = dGamma_qq->at(position);
      else 
	rate = 0.;

      if (iAlphas+1 < Nalphas) 
	rateAlphaUp = dGamma_qq->at(positionAlphaUp);
      else 
	rateAlphaUp = rate;

      if (iOmega+1 < Nomega)
	rateOmegaUp = dGamma_qq->at(positionOmegaUp);
      else 
	rateOmegaUp = rate;

      if (iAlphas < Nalphas && iOmega < Nomega)
	rateAlphaUpOmegaUp = dGamma_qq->at(positionAlphaUpOmegaUp);
      else
	rateAlphaUpOmegaUp = rate;
    }
  else 
    {
      if (position > 0 && iAlphas<Nalphas && iOmega<Nomega)
	rate = dGamma_qg->at(position);
      else
	rate = 0.;

      if (iAlphas+1 < Nalphas)
	rateAlphaUp = dGamma_qg->at(positionAlphaUp);
      else
	rateAlphaUp = rate;

      if (iOmega+1 < Nomega)
	rateOmegaUp = dGamma_qg->at(positionOmegaUp);
      else
	rateOmegaUp = rate;

      if (iAlphas < Nalphas && iOmega < Nomega)
	rateAlphaUpOmegaUp = dGamma_qg->at(positionAlphaUpOmegaUp);
      else
	rateAlphaUpOmegaUp = rate;
    }
  
  if (omega > 0.)
    {
      rateOmegaAv = (1.-omegaFrac)*rate + omegaFrac*rateOmegaUp;
      rateAlphaUpOmegaAv = (1.-omegaFrac)*rateAlphaUp + omegaFrac*rateAlphaUpOmegaUp;
    }
  else
    {
      rateOmegaAv = (omegaFrac)*rate + (1.-omegaFrac)*rateOmegaUp;
      rateAlphaUpOmegaAv = (omegaFrac)*rateAlphaUp + (1.-omegaFrac)*rateAlphaUpOmegaUp;
    }
  result = (1.-alphaFrac)*rateOmegaAv + alphaFrac*rateAlphaUpOmegaAv;

  // leave out the *9./4. for processes 6 and 8 to use the same envelope later
  return result; 
}

double Martini::use_elastic_table_q(double omega, double q, int which_kind)
/* Uses the lookup table and simple interpolation to get the value
   of dGamma/domegadq at some value of omega, q, and alpha_s. */
{
  double result;
  double alphaFrac, omegaFrac, qFrac;
  int iOmega, iAlphas, iQ;
  int position, positionAlphaUp, positionOmegaUp, positionAlphaUpOmegaUp;
  int positionQUp, positionAlphaUpQUp, positionOmegaUpQUp, positionAlphaUpOmegaUpQUp;
  int position2QUp;
  double rate, rateAlphaUp, rateOmegaUp, rateAlphaUpOmegaUp;
  double rateQUp, rateAlphaUpQUp, rateOmegaUpQUp, rateAlphaUpOmegaUpQUp;
  double rate2QUp, rateAlphaUp2QUp, rateOmegaUp2QUp, rateAlphaUpOmegaUp2QUp;
  double rateOmegaAv, rateAlphaUpOmegaAv, rateQUpOmegaAv, rateAlphaUpQUpOmegaAv;
  double rate2QUpOmegaAv, rateAlphaUp2QUpOmegaAv;
  double rateQAv, rateAlphaUpQAv;
  double slope, slopeAlphaUp;

  rate2QUp = 0.;
  rateAlphaUp2QUp = 0.;
  rateOmegaUp2QUp = 0.;
  rateAlphaUpOmegaUp2QUp = 0.;
  rateOmegaAv = 0.;
  rateAlphaUpOmegaAv = 0.;
  rateQUpOmegaAv = 0.;
  rateAlphaUpQUpOmegaAv = 0.;
  rate2QUpOmegaAv = 0.;
  rateAlphaUp2QUpOmegaAv = 0.;

  if (omega > 0.) iOmega = Nomega/2+floor((log(omega)+5)/omegaStep);
  else iOmega = Nomega/2-ceil((log(-omega)+5)/omegaStep)-1;
  iQ = floor((log(q)+5)/qStep+0.0001);
  iAlphas = floor((alpha_s-0.15)/alphaStep+0.0001);

  position = iQ + Nq*(iOmega+Nomega*(iAlphas));
  positionAlphaUp = iQ + Nq*(iOmega+Nomega*(iAlphas+1));
  positionOmegaUp = iQ + Nq*(iOmega+1+Nomega*(iAlphas));
  positionQUp = iQ+1 + Nq*(iOmega+Nomega*(iAlphas));
  position2QUp = iQ+2 + Nq*(iOmega+Nomega*(iAlphas));
  positionAlphaUpOmegaUp = iQ + Nq*(iOmega+1+Nomega*(iAlphas+1));
  positionAlphaUpQUp = iQ+1 + Nq*(iOmega+Nomega*(iAlphas+1));
  positionOmegaUpQUp = iQ+1 + Nq*(iOmega+1+Nomega*(iAlphas));
  positionAlphaUpOmegaUpQUp = iQ+1 + Nq*(iOmega+1+Nomega*(iAlphas+1));

  alphaFrac = (alpha_s-(floor((alpha_s-alphaMin)/alphaStep)*alphaStep+alphaMin))/alphaStep;
  if (omega > 0.) 
    {
      if (exp(ceil((log(omega)+5)/omegaStep)*omegaStep-5)
	  != exp(floor((log(omega)+5)/omegaStep)*omegaStep-5))
	omegaFrac = (omega-(exp(floor((log(omega)+5)/omegaStep)*omegaStep-5)))
	  /((exp(ceil((log(omega)+5)/omegaStep)*omegaStep-5))
	    -exp(floor((log(omega)+5)/omegaStep)*omegaStep-5));
      else omegaFrac = 0.;
    }
  else
    {
      if (exp(ceil((log(-omega)+5)/omegaStep)*omegaStep-5)
	  != exp(floor((log(-omega)+5)/omegaStep)*omegaStep-5))
	omegaFrac = (-omega-(exp(floor((log(-omega)+5)/omegaStep)*omegaStep-5)))
	  /((exp(ceil((log(-omega)+5)/omegaStep)*omegaStep-5))
	    -exp(floor((log(-omega)+5)/omegaStep)*omegaStep-5));
      else omegaFrac = 0.;
    }

  // interpolate the logs linearly for large omegas 
  // since there the spectrum is dropping exp in q
  if (omega > 20.) 
    {
      qFrac = (log(q)-(floor((log(q)+5.)/qStep)*qStep-5.))/qStep;
    }
  // linear interpolation in q for small omegas
  else 
    {
      if (exp(ceil((log(q)+5.)/qStep)*qStep-5.)
	  !=exp(floor((log(q)+5.)/qStep)*qStep-5.))
	qFrac = (q-(exp(floor((log(q)+5.)/qStep)*qStep-5.)))
	  /((exp(ceil((log(q)+5.)/qStep)*qStep-5.))
	    -exp(floor((log(q)+5.)/qStep)*qStep-5.));
      else qFrac = 0.;
    }

  if (which_kind == 5 || which_kind == 7)
    {
      if (position >= 0 && iAlphas<Nalphas && iOmega<Nomega && iQ < Nq )
	rate = dGamma_qq_q->at(position);
      else
	rate = 0.;

      if (iAlphas+1 < Nalphas)
	rateAlphaUp = dGamma_qq_q->at(positionAlphaUp);
      else
	rateAlphaUp = rate;

      if (iOmega+1 < Nomega)
	rateOmegaUp = dGamma_qq_q->at(positionOmegaUp);
      else
	rateOmegaUp = rate;

      if (iQ+1 < Nq)
	rateQUp = dGamma_qq_q->at(positionQUp);
      else
	rateQUp = rate;

      if (iAlphas < Nalphas && iOmega < Nomega)
	rateAlphaUpOmegaUp = dGamma_qq_q->at(positionAlphaUpOmegaUp);
      else
	rateAlphaUpOmegaUp = rate;

      if (iAlphas < Nalphas && iQ < Nq)
	rateAlphaUpQUp = dGamma_qq_q->at(positionAlphaUpQUp);
      else
	rateAlphaUpQUp = rate;

      if (iOmega+1 < Nomega && iQ+1 < Nq)
	rateOmegaUpQUp = dGamma_qq_q->at(positionOmegaUpQUp);
      else
	rateOmegaUpQUp = rate;

      if (iAlphas < Nalphas && iOmega < Nomega && iQ < Nq)
	rateAlphaUpOmegaUpQUp = dGamma_qq_q->at(positionAlphaUpOmegaUpQUp);
      else
	rateAlphaUpOmegaUpQUp = rate;

      // used for extrapolation when the data points are too far apart
      if (omega > 20.)
	{ 
	  if (iQ+2 < Nq )
	    rate2QUp = dGamma_qq_q->at(position2QUp);
	  else
	    rate2QUp = rateQUp;

	  if (iAlphas < Nalphas && iQ+2 < Nq )
	    rateAlphaUp2QUp = dGamma_qq_q->at(positionAlphaUpQUp+1);
	  else
	    rateAlphaUp2QUp = rateAlphaUpQUp;

	  if (iOmega < Nomega && iQ+2 < Nq )
	    rateOmegaUp2QUp = dGamma_qq_q->at(positionOmegaUpQUp+1);
	  else
	    rateOmegaUp2QUp = rateOmegaUpQUp;

	  if (iAlphas < Nalphas && iOmega < Nomega && iQ+2 < Nq )
	    rateAlphaUpOmegaUp2QUp = dGamma_qq_q->at(positionAlphaUpOmegaUpQUp+1);
	  else
	    rateAlphaUpOmegaUp2QUp = rateAlphaUpOmegaUpQUp;
	}
    }
  else
    {
      if (position > 0 && iAlphas < Nalphas && iOmega < Nomega && iQ < Nq )
	rate = dGamma_qg_q->at(position);
      else
	rate = 0.;

      if (iAlphas+1 < Nalphas)
	rateAlphaUp = dGamma_qg_q->at(positionAlphaUp);
      else
	rateAlphaUp = rate;

      if (iOmega+1 < Nomega)
	rateOmegaUp = dGamma_qg_q->at(positionOmegaUp);
      else
	rateOmegaUp = rate;

      if (iQ+1 < Nq)
	rateQUp = dGamma_qg_q->at(positionQUp);
      else
	rateQUp = rate;

      if (iAlphas < Nalphas && iOmega < Nomega)
	rateAlphaUpOmegaUp = dGamma_qg_q->at(positionAlphaUpOmegaUp);
      else
	rateAlphaUpOmegaUp = rate;

      if (iAlphas < Nalphas && iQ < Nq)
	rateAlphaUpQUp = dGamma_qg_q->at(positionAlphaUpQUp);
      else
	rateAlphaUpQUp = rate;

      if (iOmega+1 < Nomega && iQ+1 < Nq)
	rateOmegaUpQUp = dGamma_qg_q->at(positionOmegaUpQUp);
      else
	rateOmegaUpQUp = rate;

      if (iAlphas < Nalphas && iOmega < Nomega && iQ < Nq)
	rateAlphaUpOmegaUpQUp = dGamma_qg_q->at(positionAlphaUpOmegaUpQUp);
      else
	rateAlphaUpOmegaUpQUp = rate;

      // used for extrapolation when the data points are too far apart
      if (omega > 20.)
	{ 
	  if (iQ+2 < Nq )
	    rate2QUp = dGamma_qg_q->at(position2QUp);
	  else
	    rate2QUp = rateQUp;

	  if (iAlphas < Nalphas && iQ+2 < Nq )
	    rateAlphaUp2QUp = dGamma_qg_q->at(positionAlphaUpQUp+1);
	  else
	    rateAlphaUp2QUp = rateAlphaUpQUp;

	  if (iOmega < Nomega && iQ+2 < Nq )
	    rateOmegaUp2QUp = dGamma_qg_q->at(positionOmegaUpQUp+1);
	  else
	    rateOmegaUp2QUp = rateOmegaUpQUp;

	  if (iAlphas < Nalphas && iOmega < Nomega && iQ+2 < Nq )
	    rateAlphaUpOmegaUp2QUp = dGamma_qg_q->at(positionAlphaUpOmegaUpQUp+1);
	  else
	    rateAlphaUpOmegaUp2QUp = rateAlphaUpOmegaUpQUp;
	}
    }
  
  if (omega > 0. && omega <= 20.)
    {
      rateOmegaAv = (1.-omegaFrac)*rate + omegaFrac*rateOmegaUp;
      rateAlphaUpOmegaAv = (1.-omegaFrac)*rateAlphaUp + omegaFrac*rateAlphaUpOmegaUp;
      rateQUpOmegaAv = (1.-omegaFrac)*rateQUp + omegaFrac*rateOmegaUpQUp;
      rateAlphaUpQUpOmegaAv = (1.-omegaFrac)*rateAlphaUpQUp + omegaFrac*rateAlphaUpOmegaUpQUp;
    }
  else if (omega > 20.)
    {
      if (rate != 0. && rateOmegaUp != 0.)
	rateOmegaAv = exp((1.-omegaFrac)*log(rate)+omegaFrac*log(rateOmegaUp));
      else if (rate == 0.)
	rateOmegaAv = rateOmegaUp;
      else if (rateOmegaUp == 0.)
	rateOmegaAv = rate;
      else 
	rateOmegaAv = 0.;

      if (rateAlphaUpOmegaUp != 0. && rateAlphaUp != 0.)
	rateAlphaUpOmegaAv = exp((1.-omegaFrac)*log(rateAlphaUp)
				 +omegaFrac*log(rateAlphaUpOmegaUp));
      else if (rateAlphaUp == 0.)
	rateAlphaUpOmegaAv = rateAlphaUpOmegaUp;
      else if (rateAlphaUpOmegaUp == 0.)
	rateAlphaUpOmegaAv = rateAlphaUp;
      else 
	rateAlphaUpOmegaAv = 0.;

      if (rateOmegaUpQUp != 0. && rateQUp != 0.)
	rateQUpOmegaAv = exp((1.-omegaFrac)*log(rateQUp)
			     +omegaFrac*log(rateOmegaUpQUp));
      else if (rateOmegaUpQUp == 0.)
	rateQUpOmegaAv = rateQUp;
      else if (rateQUp == 0.)
	rateQUpOmegaAv = rateOmegaUpQUp;
      else 
	rateQUpOmegaAv = 0.;

      if (rateAlphaUpOmegaUpQUp != 0. && rateAlphaUpQUp != 0.)
	rateAlphaUpQUpOmegaAv = exp((1.-omegaFrac)*log(rateAlphaUpQUp)
				    +omegaFrac*log(rateAlphaUpOmegaUpQUp));
      else if (rateAlphaUpQUp == 0.)
	rateAlphaUpQUpOmegaAv = rateAlphaUpOmegaUpQUp;
      else if (rateAlphaUpOmegaUpQUp == 0.)
	rateAlphaUpQUpOmegaAv = rateAlphaUpQUp;
      else 
	rateAlphaUpQUpOmegaAv = 0.;

      rate2QUpOmegaAv = exp((1.-omegaFrac)*log(rate2QUp)+omegaFrac*log(rateOmegaUp2QUp));
      rateAlphaUp2QUpOmegaAv = exp((1.-omegaFrac)*log(rateAlphaUp2QUp)
				   +omegaFrac*log(rateAlphaUpOmegaUp2QUp));
    }
  else if (omega < 0.)
    {
      rateOmegaAv = (omegaFrac)*rate+(1.-omegaFrac)*rateOmegaUp;
      rateQUpOmegaAv = (omegaFrac)*rateQUp+(1.-omegaFrac)*rateOmegaUpQUp;
      rateAlphaUpOmegaAv = (omegaFrac)*rateAlphaUp+(1.-omegaFrac)*rateAlphaUpOmegaUp;
      rateAlphaUpQUpOmegaAv = (omegaFrac)*rateAlphaUpQUp+(1.-omegaFrac)*rateAlphaUpOmegaUpQUp;
    }      

  // interpolate logs for large omega
  if (omega > 20.)
    {
      if (rateOmegaAv > 0.)
	{
	  rateQAv = exp((1.-qFrac)*log(rateOmegaAv)+qFrac*log(rateQUpOmegaAv));
	}
      else if (rateOmegaAv < 0.)  // use extrapolation
	{
	  slope = (log(rate2QUpOmegaAv)-log(rateQUpOmegaAv))/qStep;
	  rateQAv = exp(log(rateQUpOmegaAv)-slope*((1.-qFrac)*qStep));
	}
      else
	{
	  rateQAv = 0.;
	}

      if (rateAlphaUpOmegaAv > 0.)
	{
	  rateAlphaUpQAv = exp((1.-qFrac)*log(rateAlphaUpOmegaAv) + qFrac*log(rateAlphaUpQUpOmegaAv));
	}
      else if (rateAlphaUpOmegaAv < 0.)  // use extrapolation
	{
	  slopeAlphaUp = (log(rateAlphaUp2QUpOmegaAv)-log(rateAlphaUpQUpOmegaAv))/qStep;
	  rateAlphaUpQAv = exp(log(rateAlphaUpQUpOmegaAv)-slopeAlphaUp*((1.-qFrac)*qStep));
	}
      else
	{
	  rateAlphaUpQAv = 0.;
	}
    }
  // interpolate linearly for small omega
  else
    {
      rateQAv = (1.-qFrac)*rateOmegaAv + qFrac*rateQUpOmegaAv;
      rateAlphaUpQAv = (1.-qFrac)*rateAlphaUpOmegaAv + qFrac*rateAlphaUpQUpOmegaAv;
    }

  result = (1.-alphaFrac)*rateQAv + alphaFrac*rateAlphaUpQAv;
  
  // the absolute normalization doesn't matter when sampling the shape
  // it only matters in "totalRate" etc.
  return result;
}

double LambertW(double z)
{
  double w_new, w_old, ratio, e_old, tol;
  int n;

  tol = 1.0e-14;

  if(z <= -exp(-1.0))
    {
      JSWARN << "LambertW is not defined for z = " << z;
      JSWARN << "z needs to be bigger than " << -exp(-1.0);
      throw std::runtime_error("LambertW small z problem");
    }

  if(z > 3.0)
    {
      w_old = log(z) - log(log(z));
    }
  else {w_old = 1.0;}
 
  w_new = 0.0;
  ratio = 1.0;
  n = 0;
  while(std::abs(ratio) > tol) 
    {
      e_old = exp(w_old);
      ratio = w_old*e_old - z;
      ratio /= ( e_old*(w_old + 1.0) - (w_old+2.0)*(w_old*e_old-z)/(2.0*w_old + 2.0) );
      w_new = w_old - ratio;
      w_old = w_new;
      n++;
      if(n > 99) 
	{
          JSWARN << "LambertW is not converging after 100 iterations.";
          JSWARN << "LambertW: z = " << z;
          JSWARN << "LambertW: w_old = " << w_old;
          JSWARN << "LambertW: w_new = " << w_new;
          JSWARN << "LambertW: ratio = " << ratio;
          throw std::runtime_error("LambertW not conversing");
	}
    }

  return w_new;
}// LambertW by Sangyong Jeon
