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

#include "Matter.h"
#include "JetScapeLogger.h"
#include "JetScapeParticles.h"
#include "Pythia8/Pythia.h"

#include <string>

#include<iostream>

#include "FluidDynamics.h"
#include <GTL/dfs.h>

#define MAGENTA "\033[35m"

using namespace Jetscape;
using namespace std;

const double QS = 1.0 ;

Matter::Matter()
{
  SetId("Matter");
  VERBOSE(8);
}

Matter::~Matter()
{
  VERBOSE(8);
}

void Matter::Init()
{
  JSINFO<<"Intialize Matter ...";

  // Redundant (get this from Base) quick fix here for now
  tinyxml2::XMLElement *eloss= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Eloss" );
  if ( !eloss ) {
    JSWARN << "Couldn't find tag Eloss";
    throw std::runtime_error ("Couldn't find tag Eloss");    
  }
  tinyxml2::XMLElement *matter=eloss->FirstChildElement("Matter");
  if ( !matter ) {
    JSWARN << "Couldn't find tag Eloss -> Matter";
    throw std::runtime_error ("Couldn't find tag Eloss -> Matter");    
  }
 
  if (matter) {   
    string s = matter->FirstChildElement( "name" )->GetText();
    JSDEBUG << s << " to be initializied ...";
 
    in_vac = false;
    brick_med = true;
    recoil_on = false;

    qhat = 0.0;
    Q00 = 1.0; // virtuality separation scale
    qhat0 = 2.0; // GeV^2/fm for gluon at s = 96 fm^-3
    alphas = 0.3; // only useful when qhat0 is a negative number
    hydro_Tc = 0.16;
    brick_length = 4.0;
    vir_factor = 1.0;
      MaxColor = 101;

    double m_qhat=-99.99;
    matter->FirstChildElement("qhat0")->QueryDoubleText(&m_qhat);
    SetQhat(m_qhat);
    //qhat = GetQhat()/fmToGeVinv ;
    //qhat0 = GetQhat()/fmToGeVinv ;
    qhat0 = GetQhat();
    JSDEBUG  << s << " with qhat0 = "<<GetQhat();
 
    int flagInt=-100;
    double inputDouble=-99.99;
      
      if ( !matter->FirstChildElement("matter_on") ) {
          JSWARN << "Couldn't find sub-tag Eloss -> Matter -> matter_on";
          throw std::runtime_error ("Couldn't find sub-tag Eloss -> Matter -> matter_on");
      }
      matter->FirstChildElement("matter_on")->QueryIntText(&flagInt);
      matter_on = flagInt;

    if ( !matter->FirstChildElement("in_vac") ) {
	JSWARN << "Couldn't find sub-tag Eloss -> Matter -> in_vac";
        throw std::runtime_error ("Couldn't find sub-tag Eloss -> Matter -> in_vac");
    }
    matter->FirstChildElement("in_vac")->QueryIntText(&flagInt);
    in_vac = flagInt;

    if ( !matter->FirstChildElement("recoil_on") ) {
	JSWARN << "Couldn't find sub-tag Eloss -> Matter -> recoil_on";
        throw std::runtime_error ("Couldn't find sub-tag Eloss -> Matter -> recoil_on");
    }
    matter->FirstChildElement("recoil_on")->QueryIntText(&flagInt);
    recoil_on = flagInt;

    if ( !matter->FirstChildElement("broadening_on") ) {
	JSWARN << "Couldn't find sub-tag Eloss -> Matter -> broadening_on";
        throw std::runtime_error ("Couldn't find sub-tag Eloss -> Matter -> broadening_on");
    }
    matter->FirstChildElement("broadening_on")->QueryIntText(&flagInt);
    broadening_on = flagInt;


    if ( !matter->FirstChildElement("brick_med") ) {
	JSWARN << "Couldn't find sub-tag Eloss -> Matter -> brick_med";
        throw std::runtime_error ("Couldn't find sub-tag Eloss -> Matter -> brick_med");
    }
    matter->FirstChildElement("brick_med")->QueryIntText(&flagInt);
    brick_med = flagInt;

    if ( !matter->FirstChildElement("Q0") ) {
	JSWARN << "Couldn't find sub-tag Eloss -> Matter -> Q0";
        throw std::runtime_error ("Couldn't find sub-tag Eloss -> Matter -> Q0");
    }
    matter->FirstChildElement("Q0")->QueryDoubleText(&inputDouble);
    Q00 = inputDouble;

    if ( !matter->FirstChildElement("T0") ) {
	JSWARN << "Couldn't find sub-tag Eloss -> Matter -> T0";
        throw std::runtime_error ("Couldn't find sub-tag Eloss -> Matter -> T0");
    }
    matter->FirstChildElement("T0")->QueryDoubleText(&inputDouble);
    T0 = inputDouble;

    if ( !matter->FirstChildElement("alphas") ) {
	JSWARN << "Couldn't find sub-tag Eloss -> Matter -> alphas";
        throw std::runtime_error ("Couldn't find sub-tag Eloss -> Matter -> alphas");
    }
    matter->FirstChildElement("alphas")->QueryDoubleText(&inputDouble);
    alphas = inputDouble;

    if ( !matter->FirstChildElement("hydro_Tc") ) {
	JSWARN << "Couldn't find sub-tag Eloss -> Matter -> hydro_Tc";
        throw std::runtime_error ("Couldn't find sub-tag Eloss -> Matter -> hydro_Tc");
    }
    matter->FirstChildElement("hydro_Tc")->QueryDoubleText(&inputDouble);
    hydro_Tc = inputDouble;

    if ( !matter->FirstChildElement("brick_length") ) {
	JSWARN << "Couldn't find sub-tag Eloss -> Matter -> brick_length";
        throw std::runtime_error ("Couldn't find sub-tag Eloss -> Matter -> brick_length");
    }
    matter->FirstChildElement("brick_length")->QueryDoubleText(&inputDouble);
    brick_length = inputDouble;

    if ( !matter->FirstChildElement("vir_factor") ) {
	JSWARN << "Couldn't find sub-tag Eloss -> Matter -> vir_factor";
        throw std::runtime_error ("Couldn't find sub-tag Eloss -> Matter -> vir_factor");
    }
    matter->FirstChildElement("vir_factor")->QueryDoubleText(&inputDouble);
    vir_factor = inputDouble;

    if(vir_factor<0.0) {
        cout << "Reminder: negative vir_factor is set, initial energy will be used as initial t_max" << endl;
    }

    VERBOSE(7)<< MAGENTA << "MATTER input parameter";
    VERBOSE(7) << "in_vac: " << in_vac << "  brick_med: " << brick_med << "  recoil_on: " << recoil_on;
    VERBOSE(7) << "Q0: " << Q00 << " vir_factor: " << vir_factor << "  qhat0: " << qhat0 << " alphas: " << alphas << " hydro_Tc: " << hydro_Tc << " brick_length: " << brick_length;

  }
  else {
    JSWARN << " : Matter not properly initialized in XML file ...";
    throw std::runtime_error("Matter not properly initialized in XML file ...");
  }

  // Initialize random number distribution
  ZeroOneDistribution = uniform_real_distribution<double> { 0.0, 1.0 };

  //...initialize the random number generator
  srand((unsigned)time(NULL));
  NUM1=-1*rand();
  //    NUM1=-33;
    iEvent = 0;
}


void Matter::WriteTask(weak_ptr<JetScapeWriter> w)
{
  VERBOSE(8);  auto f = w.lock();
  if ( !f ) return;
  f->WriteComment("ElossModule Parton List: "+GetId());
  f->WriteComment("Energy loss to be implemented accordingly ...");
}

void Matter::DoEnergyLoss(double deltaT, double time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut)
{

  double z=0.5;
  double blurb,zeta,tQ2 ;
  int iSplit,pid_a,pid_b;
  unsigned int max_color, min_color, min_anti_color;
  double velocity[4],xStart[4],velocity_jet[4];
  bool photon_brem = false;
    
  VERBOSESHOWER(9)<< MAGENTA << "SentInPartons Signal received : "<<deltaT<<" "<<Q2<<" "<< pIn.size();
    
  VERBOSE(8) << BOLDYELLOW << " ********************************************************** " ;
    
  double rNum;
        
  double delT = deltaT;
  double Time = time*fmToGeVinv;
  double deltaTime = delT*fmToGeVinv;
  double ehat=0;
  double ehat_over_T2 = 10.0;

  std::unique_ptr<FluidCellInfo> check_fluid_info_ptr;


 VERBOSE(7)<< " the time in fm is " << time << " The time in GeV-1 is " << Time ;
 VERBOSE(7) << "pid = " << pIn[0].pid() << " E = " << pIn[0].e() << " px = " << pIn[0].p(1) << " py = " << pIn[0].p(2) << "  pz = " << pIn[0].p(3) << " virtuality = " << pIn[0].t() << " form_time in fm = " << pIn[0].form_time() << " split time = " << pIn[0].form_time() + pIn[0].x_in().t();
  JSDEBUG << " color = " << pIn[0].color() << " anti-color = " << pIn[0].anti_color();
    
    unsigned int ShowerMaxColor = pIn[0].max_color();
    unsigned int CurrentMaxColor;

    if ( pIn[0].max_color() < MaxColor )
    {
        pIn[0].set_max_color(MaxColor);
        
    }
    else
    {
        MaxColor = pIn[0].max_color();
    }
    
    //JSINFO << MAGENTA << " Max color = " << MaxColor;
   //JSDEBUG << " For MATTER, the qhat in GeV^-3 = " << qhat ;
 
    double qhatbrick;
    if(brick_med) qhatbrick=qhat0/3.0;
    qhat = qhat0;

   VERBOSE(8) << " qhat0 = " << qhat0 << " qhat = " << qhat;
    
    
  for (int i=0;i<pIn.size();i++)
  {


      // Reject photons
      if (pIn[i].pid()==photonid)
      {
          JSINFO << BOLDYELLOW << " A photon was RECEIVED with px = " << pIn[i].px() << " from framework and sent back " ;

          pOut.push_back(pIn[i]);
          return;
      }
      
      VERBOSE(2) << BOLDYELLOW << " *  parton formation spacetime point= "<< pIn[i].x_in().t() << "  " << pIn[i].x_in().x() << "  " << pIn[i].x_in().y() << "  " << pIn[i].x_in().z();

      //JSINFO << MAGENTA << " particle rest mass = " << pIn[i].restmass();


      int jet_stat=pIn[i].pstat();  // daughter of recoil will always be recoil
      
      JSDEBUG << MAGENTA << "MATTER -- status: " << pIn[i].pstat() << "  energy: " << pIn[i].e() << " color: " << pIn[i].color() << "  " << pIn[i].anti_color() << "  clock: " << time ;

      velocity[0] = 1.0;
      // Define 3 velocity of the parton
      for(int j=1;j<=3;j++)
      {
          velocity[j] = pIn[i].p(j)/pIn[i].e();
      }
      // velocityMod will be 1 for most partons except maybe heavy quarks, we say "maybe", as for very
      // energetic heavy quarks, the velocity may be very close to 1.
      double velocityMod = std::sqrt(std::pow(velocity[1],2) + std::pow(velocity[2],2) + std::pow(velocity[3],2));

      VERBOSE(2) << BOLDYELLOW << " velocityMod = " << velocityMod ;
      
      if(pIn[i].form_time()<0.0) pIn[i].set_jet_v(velocity); // jet velocity is set only once
      // Notice the assumption that partons passed from hard scattering are on shell.
      // If they had virtuality then this would not be correct. 

      // Define a vector in the direction of the jet originating parton.
      // there is some amount of redundancy here, pIn[i].jet_v() is basically the jet velocity
      // we are defining a local copy in velocity_jet
      velocity_jet[0]=1.0;
      velocity_jet[1]=pIn[i].jet_v().x();
      velocity_jet[2]=pIn[i].jet_v().y();
      velocity_jet[3]=pIn[i].jet_v().z();
    
      // Modulus of the vector pIn[i].jet_v
      double mod_jet_v = std::sqrt( pow(pIn[i].jet_v().x(),2) +  pow(pIn[i].jet_v().y(),2) + pow(pIn[i].jet_v().z(),2) );
      
      for(int j=0;j<=3;j++)
      {
          xStart[j] = pIn[i].x_in().comp(j) ;
      }

      // SC: read in hydro
      initR0 = xStart[0];
      initRx = xStart[1];
      initRy = xStart[2];
      initRz = xStart[3];
      initVx = velocity[1]/velocityMod;
      initVy = velocity[2]/velocityMod;
      initVz = velocity[3]/velocityMod;
      initRdotV = ( initRx*pIn[i].jet_v().x() + initRy*pIn[i].jet_v().y() + initRz*pIn[i].jet_v().z() )/mod_jet_v;
      initVdotV = ( initVx*pIn[i].jet_v().x() + initVy*pIn[i].jet_v().y() + initVz*pIn[i].jet_v().z() )/mod_jet_v;
      // Note: jet_v()/mod_jet_v is a unit 3 vector in the direction of the jet originating parton.
      
      initEner = pIn[i].e(); // initial Energy of parton
      if(!in_vac) length = fillQhatTab();
      if(brick_med) length = brick_length*fmToGeVinv; /// length in GeV-1 will have to changed for hydro
      //if(brick_med) length = 5.0*fmToGeVinv; /// length in GeV-1 will have to changed for hydro

      // SC
      zeta = ( ( xStart[0] + initRdotV )/std::sqrt(2) )*fmToGeVinv;
      
      VERBOSE(8)<< BOLDYELLOW << " zeta = " << zeta;

      double now_R0 = time;
      double now_Rx = initRx+(time-initR0)*initVx;
      double now_Ry = initRy+(time-initR0)*initVy; 
      double now_Rz = initRz+(time-initR0)*initVz;
      double now_temp; 

      if(!in_vac && now_R0>=tStart)
      {
         if(now_R0*now_R0<now_Rz*now_Rz) cout << "Warning 1: " << now_R0 << "  " << now_Rz << endl;
         GetHydroCellSignal(now_R0, now_Rx, now_Ry, now_Rz, check_fluid_info_ptr);
         //VERBOSE(8)<<MAGENTA<<"Temperature from medium = "<<check_fluid_info_ptr->temperature;
         now_temp = check_fluid_info_ptr->temperature;
      }
      else
      {
         now_temp = 0.0;
      }

      int pid = pIn[i].pid();
      
      if (pIn[i].form_time()<0.0) /// A parton without a virtuality or formation time, must set...
      {
          
        //  if (pIn[i].t()<0.0)
        //  {
        //      JSWARN << " parton with a negative virtuality was sent to MATTER and will now have its virtuality reset!, press 1 and return to proceed... ";
        //      cin >> blurb;
        //  }
          
          iSplit = 0;
          if (pIn[i].pid()==gid)
          {
              JSDEBUG << " parton is a gluon ";
              iSplit = 1;
          }
          else
          {
              JSDEBUG << " parton is a quark ";
          }
	  
          //tQ2 = generate_vac_t(pIn[i].pid(), pIn[i].nu(), QS/2.0, pIn[i].e()*pIn[i].e() ,zeta , iSplit);

	  // SC:  
          double pT2 = pIn[i].p(1)*pIn[i].p(1)+pIn[i].p(2)*pIn[i].p(2);
          double max_vir;
          if(vir_factor<0.0) max_vir = pIn[i].e()*pIn[i].e() - pIn[i].restmass()*pIn[i].restmass();
          else max_vir = pT2 * vir_factor;

          if(max_vir<=QS)
          {
              tQ2 = 0.0;
          }
          else
          {
             VERBOSE(8) << BOLDYELLOW << " at x,y,z,t = " <<  pIn[i].x_in().x() << "  " << pIn[i].x_in().y() << "  " << pIn[i].x_in().z() << "  " << pIn[i].x_in().t() ;
             if(abs(pIn[i].pid()) == 4 || abs(pIn[i].pid()) == 5)
             { 
                 tQ2 = generate_vac_t_w_M(pIn[i].pid(), pIn[i].restmass(), pIn[i].nu(), QS/2.0, max_vir, zeta, iSplit);
                 
                 //JSINFO  << BOLDYELLOW << " virtuality calculated as = " << tQ2;
             }
	     else
             {
	         tQ2 = generate_vac_t(pIn[i].pid(), pIn[i].nu(), QS/2.0, max_vir, zeta, iSplit);
             }
    	  }


	  // SC: if matter_on = false, set zero virtuality and MATTER will not do parton shower
          if(matter_on)
          {
              pIn[i].set_t(tQ2); // Also resets momentum!
              //JSINFO << BOLDYELLOW << " virtuality set to " << tQ2 ;
          }
          else pIn[i].set_t(rounding_error);
          //else pIn[i].set_t(0.0);


          pIn[i].set_mean_form_time();
          double ft = generate_L(pIn[i].mean_form_time());
          pIn[i].set_form_time(ft);
          
          unsigned int color=0, anti_color=0;
          std::uniform_int_distribution<short> uni(102,103);
          
          if ( pIn[i].pid()>0 )
          {
              color = 101;
          }
          pIn[i].set_color(color);
          if ( (pIn[i].pid()<0)||(pIn[i].pid()==21) )
          {
              anti_color = uni(*GetMt19937Generator());
          }
          pIn[i].set_anti_color(anti_color);
          
          max_color = color;
          
          if (anti_color > color) max_color = anti_color ;
          
          min_color = color;
          
          min_anti_color = anti_color;
          
          pIn[i].set_max_color(max_color);
          pIn[i].set_min_color(min_color);
          pIn[i].set_min_anti_color(min_anti_color);
          MaxColor = max_color;
          
          // VERBOSE OUTPUT ON INITIAL STATUS OF PARTICLE:
         VERBOSE(8) ;
         VERBOSE(8) << " ***************************************************************************** " ;
         VERBOSE(8) << " ID = " << pIn[i].pid() << " Color = " << pIn[i].color() << " Anti-Color = " << pIn[i].anti_color() ;
         VERBOSE(8) << " E = " << pIn[i].e() << " px = " << pIn[i].px() << " py = " << pIn[i].py() << " pz = " << pIn[i].pz() ;
         VERBOSE(8) << " *  New generated virtuality = " << tQ2 << " Mean formation time = " << pIn[i].mean_form_time();
         VERBOSE(8) << " *  set new formation time to " << pIn[i].form_time() ;
         VERBOSE(8) << " * Maximum allowed virtuality = " << pIn[i].e()*pIn[i].e() - pIn[i].restmass()*pIn[i].restmass() << "   Minimum Virtuality = " << QS;
         VERBOSE(8) << " * Qhat = " << qhat << "  Length in fm = "  << length/5.0 ;
         VERBOSE(8) << " * Jet velocity = " << pIn[i].jet_v().comp(0) << " " << pIn[i].jet_v().comp(1) << "  " << pIn[i].jet_v().comp(2) << "  " << pIn[i].jet_v().comp(3);
         VERBOSE(8) << " * reset location of parton formation = "<< pIn[i].x_in().t() << "  " << pIn[i].x_in().x() << "  " << pIn[i].x_in().y() << "  " << pIn[i].x_in().z();
         VERBOSE(8) << " ***************************************************************************** " ;
         VERBOSE(8) ;
          // end VERBOSE OUTPUT:
 
          
      }

      
      // SC: Q0 can be changed based on different setups
      if(in_vac)
      { // for vaccuum
          qhat = 0.0;
          if(Q00 < 0.0) Q0 = 1.0; // set Q0 = 1 if Q00 < 0
          else Q0 = Q00;
      }
      else
      { // for medium
          double tempEner = initEner;
          qhat = fncQhat(zeta);
          ehat = 0.0;

          if (now_temp>0.0) ehat = 0.0*qhat/4.0/now_temp ;
         VERBOSE(8) << BOLDYELLOW << "at Origin of parton, qhat = " << qhat << " ehat = " << ehat;

// set the minimum Q0 for MATTER using Q^2 = qhat * tau  ELSE use the positive value in XML
          if(Q00 < 0.0)
          { // use dynamical Q0 if Q00 < 0
              if(pid==gid) Q0 = sqrt(sqrt(2.0*tempEner*qhat*sqrt(2.0)));
              else Q0 = sqrt(sqrt(2.0*tempEner*qhat*sqrt(2.0)/Ca*Cf));
              if(Q0 < 1.0) Q0 = 1.0;
              if(zeta > length) Q0 = 1.0;
          }
          else
          {
              Q0 = Q00;
          }
      }

// I dont care what you say, we are not doing pQCD below 1 GeV
      if(Q0<1.0) Q0=1.0;

      //if (pIn[i].t() > QS + rounding_error)
      if (pIn[i].t() > Q0*Q0 + rounding_error || ((!in_vac) && now_temp<=T0 && pIn[i].t() > QS*QS + rounding_error))
      {

          cout << "In matter ...   " << pIn[i].pid() << "  " << pIn[i].t() << endl;

          TakeResponsibilityFor ( pIn[i] ); // Generate error if another module already has responsibility.
          double decayTime = pIn[i].mean_form_time()  ;
	    
          //JSDEBUG << "  deltaT = " << deltaT;
          VERBOSE(8) << " parton origin time = " << pIn[i].x_in().t() << " parton formation time = " << pIn[i].form_time();
          VERBOSE(8) << " parton id " << pIn[i].pid() << " parton virtuality = " << pIn[i].t();
          //JSDEBUG << " parton momentum " << pIn[i].e() << "  " << pIn[i].px() << "  " << pIn[i].py() << "  " << pIn[i].pz();
	    
          double splitTime = pIn[i].form_time() + pIn[i].x_in().t() ;
          
         VERBOSE(8) << " splitTime = " << splitTime;
         VERBOSE(8) << " qhat before splitime loop = " << qhat ;

          if (splitTime<time) // it is time to split and calculate the effect of scattering
          {

             VERBOSE(8) << "SPLIT in MATTER" ;

              // SC: add elastic scattering that generates recoiled and back-reaction partons
              double el_dt=0.1;
              double el_p0[5];
              double el_CR;
              double el_rand;

              if(pid==gid) el_CR=Ca;
              else el_CR=Cf;

              for(int j=1; j<=3; j++) el_p0[j] = pIn[i].p(j);
              el_p0[0]=pIn[i].e();
              el_p0[4]=pIn[i].t();

              for(double el_time=initR0; el_time<time+rounding_error; el_time=el_time+el_dt)
              {

                  if(in_vac) continue;
                  if(!recoil_on) continue;
                  if(el_time<tStart) continue;
                  if(abs(pid)==4 || abs(pid)==5) continue; // recoil not ready for heavy quarks yet

                  double el_rx=initRx+(el_time-initR0)*initVx;
                  double el_ry=initRy+(el_time-initR0)*initVy;
                  double el_rz=initRz+(el_time-initR0)*initVz;

                  double tempLoc,sdLoc,vxLoc,vyLoc,vzLoc,qhatLoc,enerLoc;
                  double betaLoc,gammaLoc,flowFactor;
                  int hydro_ctl;
                  double dt_lrf;

                  double pc0[4]={0.0};
                  double vc0[4]={0.0};
                  double soln_alphas,prob_el;
                  double muD2;
                  double recordE0;
                  
                // Convert hard parton momentum to onshell
                  pc0[1]=el_p0[1];
                  pc0[2]=el_p0[2];
                  pc0[3]=el_p0[3];
                  pc0[0]=sqrt(pc0[1]*pc0[1]+pc0[2]*pc0[2]+pc0[3]*pc0[3]);

                  recordE0=pc0[0];

                  GetHydroCellSignal(el_time, el_rx, el_ry, el_rz, check_fluid_info_ptr);
               	  VERBOSE(8)<<MAGENTA<<"Temperature from medium = "<<check_fluid_info_ptr->temperature;
               	 	
               	  tempLoc = check_fluid_info_ptr->temperature;
               	  sdLoc = check_fluid_info_ptr->entropy_density;
               	  vxLoc = check_fluid_info_ptr->vx;
               	  vyLoc = check_fluid_info_ptr->vy;
               	  vzLoc = check_fluid_info_ptr->vz;

                  vc0[1]=vxLoc;
                  vc0[2]=vyLoc;
                  vc0[3]=vzLoc;

                  hydro_ctl=0;
               	 
               	  if(hydro_ctl == 0 && tempLoc >= hydro_Tc)
                  {
 
                      trans(vc0,pc0);
                      enerLoc=pc0[0];
                      transback(vc0,pc0);

              	      betaLoc = sqrt(vxLoc*vxLoc+vyLoc*vyLoc+vzLoc*vzLoc);
               	      gammaLoc = 1.0/sqrt(1.0-betaLoc*betaLoc);
               	      flowFactor = gammaLoc*(1.0-(initVx*vxLoc+initVy*vyLoc+initVz*vzLoc));
               	 
               	      if(qhat0 < 0.0)
                      { // calculate qhat with alphas
               	          muD2 = 6.0*pi*alphas*tempLoc*tempLoc;
               	          if(enerLoc > 2.0*pi*tempLoc) qhatLoc = Ca*50.4864/pi*pow(alphas,2)*pow(tempLoc,3)*log(5.7*enerLoc*tempLoc/4.0/muD2);
               	          else qhatLoc = Ca*50.4864/pi*pow(alphas,2)*pow(tempLoc,3)*log(5.7*2.0*pi*tempLoc*tempLoc/4.0/muD2);
               	          if(qhatLoc<0.0) qhatLoc=0.0;
               	      }
                      else
                      { // use input qhat
               	          if(brick_med) qhatLoc = qhat0*0.1973;
                          else qhatLoc = qhat0/96.0*sdLoc*0.1973;
               	      }
               	  }
                  else
                  { // outside the QGP medium
                      continue;
                  }

              // Calculating the delta t in the fluid rest frame
                  if(el_time+el_dt<=time) dt_lrf=el_dt*flowFactor;
                  else dt_lrf=(time-el_time)*flowFactor;

              // solve alphas
                  if(qhat0 < 0.0) soln_alphas=alphas;
                  else soln_alphas=solve_alphas(qhatLoc,enerLoc,tempLoc);

               // Calculate the proability of elastic scattering in time delta t in fluid rest frame
                  muD2=6.0*pi*soln_alphas*tempLoc*tempLoc;
                  prob_el=42.0*zeta3*el_CR*soln_alphas*tempLoc/6.0/pi/pi*dt_lrf/0.1973;

                  el_rand = ZeroOneDistribution(*GetMt19937Generator());

		  //cout << "  qhat: " << qhatLoc << "  alphas: " << soln_alphas << "  ener: " << enerLoc << "  prob_el: " << prob_el << "  " << el_rand << endl;

                  if(el_rand<prob_el)
                  { // elastic scattering happens
 
                      //cout << "elastic scattering happens" << endl;
                      int CT=-1;
                      int pid0=-999;
                      int pid2=-999;
                      int pid3=-999;
                      double pc2[4]={0.0}; // final recoil thermal parton
                      double pc3[4]={0.0}; // initial thermal parton
                      double pc4[4]={0.0}; // not used
                      double qt=0.0; // not used
                      unsigned int el_max_color,el_color0,el_anti_color0,el_color2,el_anti_color2,el_color3,el_anti_color3;

                      el_max_color=pIn[i].max_color();
                      el_color0 = pIn[i].color();
                      el_anti_color0 = pIn[i].anti_color();

                      pid0=pid;

                      // deterimine channel
                      //flavor(CT,pid0,pid2,pid3);
                      flavor(CT,pid0,pid2,pid3,el_max_color,el_color0,el_anti_color0,el_color2,el_anti_color2,el_color3,el_anti_color3);

                      //cout << "color: " << el_color0 << "  " << el_anti_color0 << "  " << el_color2 << "  " << el_anti_color2 << "  " << el_color3 << "  " << el_anti_color3 << "  max: " << el_max_color << endl;

                      // do scattering
                      colljet22(CT,tempLoc,muD2,vc0,pc0,pc2,pc3,pc4,qt);

                      if(pc0[0]<pc2[0] && abs(pid0)!=4 && pid0==pid2)
                      { //disable switch for heavy quark, only allow switch for identical particles
                          double p0temp[4]={0.0};
                          for(int k=0;k<=3;k++)
                          {
                              p0temp[k]=pc2[k];
                              pc2[k]=pc0[k];
                              pc0[k]=p0temp[k];
                          }
                      }

                      // push out recoied and back-reaction (negative) partons
                      // need to add color information later!!!
                      double el_vertex[4];
                      el_vertex[0] = el_time;
                      el_vertex[1] = el_rx;
                      el_vertex[2] = el_ry;
                      el_vertex[3] = el_rz;
                 
                      int iout;
                      double ft;

                      pOut.push_back(Parton(0,pid2,1,pc2,el_vertex)); // recoiled
                      iout = pOut.size()-1;
                      pOut[iout].set_jet_v(velocity_jet); // use initial jet velocity
                      pOut[iout].set_mean_form_time();
                      ft = 10000.0; /// a really large formation time.
                      pOut[iout].set_form_time(ft);
                      ////pOut[iout].set_color(el_color2);
                      ////pOut[iout].set_anti_color(el_anti_color2);
                      ////pOut[iout].set_max_color(max_color);
                      ////pOut[iout].set_min_color(pIn[i].min_color());
                      ////pOut[iout].set_min_anti_color(pIn[i].min_anti_color());
                      ////// comment out realistic color above, assume colorless for recoiled and back-reaction parton
                      ////// for the convenience of color string fragmentation in Pythia
                      pOut[iout].set_color(0);
                      pOut[iout].set_anti_color(0);
                      pOut[iout].set_max_color(pIn[i].max_color());
                      pOut[iout].set_min_color(pIn[i].min_color());
                      pOut[iout].set_min_anti_color(pIn[i].min_anti_color());
              
                      pOut.push_back(Parton(0,pid3,-1,pc3,el_vertex)); // back reaction
                      iout = pOut.size()-1;
                      pOut[iout].set_jet_v(velocity_jet);
                      pOut[iout].set_mean_form_time();
                      ft = 10000.0; /// STILL a really large formation time.
                      pOut[iout].set_form_time(ft);
                      ////pOut[iout].set_color(el_color3);
                      ////pOut[iout].set_anti_color(el_anti_color3);
                      ////pOut[iout].set_max_color(max_color);
                      ////pOut[iout].set_min_color(pIn[i].min_color());
                      ////pOut[iout].set_min_anti_color(pIn[i].min_anti_color());
                      pOut[iout].set_color(0);
                      pOut[iout].set_anti_color(0);
                      pOut[iout].set_max_color(pIn[i].max_color());
                      pOut[iout].set_min_color(pIn[i].min_color());
                      pOut[iout].set_min_anti_color(pIn[i].min_anti_color());

                      ////// assumption: pIn doesn't change color in elastic scattering
                      ////pIn[i].set_color(el_color0);
                      ////pIn[i].set_anti_color(el_anti_color0);
                      ////pIn[i].set_max_color(el_max_color);
 
                  } // end if scattering

                  // E1'+E2 = E3'+E4 (all massless in scattering)
		  // Now I want E1+E2 = E3+E4 where 1 and 3 have mass
		  // -> E1-E1' = E3-E3' -> E3 = E1-E1'+E3'
		  // m3^2 = E3^2-E3'^2
                  for(int j=1; j<=3; j++) el_p0[j] = pc0[j];
		  
                  el_p0[0]=el_p0[0]-recordE0+pc0[0];
                  el_p0[4]=el_p0[0]*el_p0[0]-pc0[0]*pc0[0];
                  if(el_p0[4]<0) cout << "complain negative virt" << endl;

              } // end time loop for elastic scattering


       	      // do split
              double t_used = pIn[i].t();
              //if (t_used<QS)  t_used = QS; // SC: not necessary
              double tau_form = 2.0*pIn[i].nu()/t_used;
              double z_low = QS/t_used/2.0 ;
              double z_hi = 1.0 - z_low;
              
              VERBOSE(8) << " zeta = " << zeta ;

              if (pid==gid)
              { // gluon
                  double val1 = P_z_gg_int(z_low, z_hi, zeta, t_used, tau_form,pIn[i].nu() );
                  double val2 = nf*P_z_qq_int(z_low, z_hi, zeta, t_used, tau_form,pIn[i].nu() );
                    
                  if ( val1<0.0 || val2<0.0 )
                  {
                      cerr << " minus log of sudakov negative val1 , val2 = " << val1 << "  " << val2 << endl;
                      throw std::runtime_error("minus log of sudakov negative");
                      // cin >> blurb ;
                  }
                  
                  double ratio = val1/(val1+val2);
                  double r = ZeroOneDistribution(*GetMt19937Generator());
                  if (r>ratio)
                  { // qqbar
		  
                      double r2 = ZeroOneDistribution(*GetMt19937Generator());
		  
                      // assign flavors
                      if (r2>0.6666)
                      {
                          pid_a = uid ;
                          pid_b = -1*uid ;
                      }
                      else if ( r2 > 0.3333 )
                      {
                          pid_a = did ;
                          pid_b = -1*did ;

                      }
                      else
                      {
                          pid_a = sid;
                          pid_b = -1*sid;
                      }
                      iSplit = 2;
                  }
                  else
                  { // gg
                      pid_a = gid ;
                      pid_b = gid ;
                      iSplit = 1;
                  }
              }
              else
              { // we had a quark
                  
                  double ProbGluon = 1.0 - sudakov_Pqg( QS/2, pIn[i].t(), zeta, pIn[i].nu() ) ;
                  double ProbPhoton = 1.0 - sudakov_Pqp( QS/2, pIn[i].t(), zeta, pIn[i].nu() ) ;
                  
                  double val = ProbGluon/(ProbGluon + ProbPhoton) ;
                  
//                  JSINFO << MAGENTA << " val = " << val ;
                  
                  double r2 = ZeroOneDistribution(*GetMt19937Generator());
                  
                  if (r2<=val)
                  { // quark decay to quark and gluon
                      pid_a = pid ;
                      pid_b = gid ;
                      iSplit = 0;
                  }
                  else
                  { // quark decay to quark and photon
                      pid_a = pid;
                      pid_b = photonid;
                      iSplit = 3;
                      photon_brem = true;
                  }
                  
                  
              }
              int ifcounter = 0;
              double l_perp2 = -1.0; // SC: initialization
              // daughter virtualities
              double tQd1 = QS;
              double tQd2 = QS;

              double new_parent_p0 = el_p0[0];
              double new_parent_px = el_p0[1];
              double new_parent_py = el_p0[2];
              double new_parent_pz = el_p0[3];
              double new_parent_t = el_p0[4];
              double new_parent_pl = ( new_parent_px*pIn[i].jet_v().x() + new_parent_py*pIn[i].jet_v().y() + new_parent_pz*pIn[i].jet_v().z() )/mod_jet_v;
              if (new_parent_pl<0.0)
              {
                  JSINFO << BOLDYELLOW << " parton traversing opposite to jet direction ... Just letting you know ! " ;
                 // cin >> blurb ;
              }
              //JSINFO << BOLDYELLOW << " old virtuality = " << pIn[i].t() << " new virtuality = " << new_parent_t ;
              //double new_parent_pl = sqrt(pow(new_parent_px,2)+pow(new_parent_py,2)+pow(new_parent_pz,2));
	      //double new_parent_vx = new_parent_px/new_parent_pl;
	      //double new_parent_vy = new_parent_py/new_parent_pl;
	      //double new_parent_vz = new_parent_pz/new_parent_pl;
              double new_parent_nu = (new_parent_p0+new_parent_pl)/sqrt(2.0);
              
              //set color of daughters here
              unsigned int d1_col, d1_acol, d2_col, d2_acol, color, anti_color;
              //std::uniform_int_distribution<short> uni(101,103);
              //color = pIn[i].color();
              
              if (iSplit<3) // not photon radiation, generate new colors
              {
                  max_color = pIn[i].max_color();
                  //if (pIn[i].anti_color()>maxcolor) color = pIn[i].anti_color();
                  JSDEBUG << " old max color = " << max_color;
                  max_color=++MaxColor;
                  color = max_color;
                  anti_color = max_color;
                  pIn[i].set_max_color(max_color);
                  JSDEBUG << " new color = " << color;
              }
              
              if (iSplit==1)///< gluon splits into two gluons
              {
                  d1_col = pIn[i].color();
                  d2_col = color;
                  d1_acol = anti_color;
                  d2_acol = pIn[i].anti_color();
              }
              else if (iSplit==0) ///< (anti-)quark splits into (anti-)quark + gluon
              {
                  if (pIn[i].pid()>0)
                  {
                      d1_col = color;
                      d1_acol = 0;
                      d2_col = pIn[i].color();
                      d2_acol = anti_color;
                  }
                  else
                  {
                      d1_col = 0; /// < gluon splits into quark anti-quark
                      d1_acol = anti_color;
                      d2_col = color;
                      d2_acol = pIn[i].anti_color();
                  }
              }
              else if (iSplit==2)
              {
                  d1_col = pIn[i].color();
                  d1_acol = 0;
                  d2_acol = pIn[i].anti_color();
                  d2_col = 0;
              }
              else if (iSplit==3) // radiating a photon has col = acol = 0, all color remains in quark
              {
                  d1_col = pIn[i].color();
                  d1_acol = pIn[i].anti_color();
                  d2_col = 0;
                  d2_acol = 0;
              }
              else
              {
                 throw std::runtime_error("error in iSplit");
              }


              while ((l_perp2<=Lambda_QCD*Lambda_QCD)&&(ifcounter<100))
              {

                  if(abs(pid) == 4 || abs(pid) == 5)
                  {
                      z = generate_vac_z_w_M(pid, pIn[i].restmass(), QS/2.0, pIn[i].t(), zeta, pIn[i].nu(), iSplit);
                  }
                  else
                  {
                      z = generate_vac_z(pid, QS/2.0, pIn[i].t(), zeta, pIn[i].nu(), iSplit);
                  }
                  //JSDEBUG << " generated z = " << z;
        
                  int iSplit_a = 0;
                  if (pid_a==gid) iSplit_a = 1;

                  // use pIn information to sample z above, but use new_parent to calculate daughter partons below            
                  if (z*z*new_parent_t>QS)
                  {
                      if(abs(pid_a) == 4 || abs(pid_a) == 5)
                      {
                          double M = pIn[i].restmass();
                          if (QS*(1.0 + std::sqrt(1.0 + 4.0*M*M/QS))/2.0 < z*z*new_parent_t)
                          {
                              tQd1 = generate_vac_t_w_M(pid_a, pIn[i].restmass(), z*new_parent_nu, QS/2.0, z*z*new_parent_t, zeta+std::sqrt(2)*pIn[i].form_time()*fmToGeVinv, iSplit_a);
                          }
                          else
                          {
                              tQd1 = z*z*new_parent_t;
                          }
                      }
                      else
                      {
                           tQd1 = generate_vac_t(pid_a, z*new_parent_nu, QS/2.0, z*z*new_parent_t, zeta+std::sqrt(2)*pIn[i].form_time()*fmToGeVinv, iSplit_a);
                      }
                  }
                  else
                  { // SC
                      tQd1 = z*z*new_parent_t;
                  }
        

                  int iSplit_b = 0;
                  if (pid_b==gid) iSplit_b = 1;
        
                  if ( ( (1.0-z)*(1.0-z)*new_parent_t>QS )&&(iSplit<3) )
                  {
                      if(abs(pid_b) == 4 || abs(pid_b) == 5)
                      {
                          double M = pIn[i].restmass();
                          
                          if (QS*(1.0 + std::sqrt(1.0 + 4.0*M*M/QS))/2.0 < (1.0-z)*(1.0-z)*new_parent_t)
                          {
                              tQd2 = generate_vac_t_w_M(pid_a, pIn[i].restmass(), (1.0-z)*new_parent_nu, QS/2.0, (1.0-z)*(1.0-z)*new_parent_t, zeta+std::sqrt(2)*pIn[i].form_time()*fmToGeVinv, iSplit_a);
                          }
                          else
                          {
                              tQd1 = (1.0-z)*(1.0-z)*new_parent_t;
                          }
                          
                      }
                      else
                      {
                           tQd2 = generate_vac_t(pid_b, (1.0-z)*new_parent_nu, QS/2.0, (1.0-z)*(1.0-z)*new_parent_t, zeta+std::sqrt(2)*pIn[i].form_time()*fmToGeVinv,iSplit_b);
                      }
                  }
                  else
                  { // SC
                      tQd2 = (1.0-z)*(1.0-z)*new_parent_t;
                  }

                  
                  if (iSplit == 3)
                  {
                    tQd2 = rounding_error; // forcing the photon to have no virtuality
                  }
                  
                  l_perp2 = new_parent_t*z*(1.0 - z) - tQd2*z - tQd1*(1.0-z) ; ///< the transverse momentum squared

		
                  if(abs(pid) == 4 || abs(pid) == 5)
                  { 
                      l_perp2 = new_parent_t*z*(1.0 - z) - tQd2*z - tQd1*(1.0-z) - std::pow((1.0-z)*pIn[i].restmass(),2); ///< the transverse momentum squared
                  }
                  else
                  {
                      l_perp2 = new_parent_t*z*(1.0 - z) - tQd2*z - tQd1*(1.0-z) ; ///< the transverse momentum squared
                  }

                  ifcounter++;
              }
              
              if (l_perp2<=Lambda_QCD*Lambda_QCD) l_perp2 = Lambda_QCD*Lambda_QCD; ///< test if negative
              double l_perp = std::sqrt(l_perp2); ///< the momentum transverse to the parent parton direction
              //JSINFO << BOLDYELLOW << " after ifcounter = " << ifcounter << " l_perp2 = " << l_perp2 ;
              //JSINFO << BOLDYELLOW << " z = " << z << " tQd1 = " << tQd1 << " tQd2 = " << tQd2 ;
              
              
              // axis of split
              double angle = generate_angle();

              // KK: changed to x,y,z
              //double parent_perp = std::sqrt( pow(pIn[i].px(),2) + pow(pIn[i].py(),2) + pow(pIn[i].pz(),2) - pow(pIn[i].pl(),2) );
              //double mod_jet_v = std::sqrt( pow(pIn[i].jet_v().x(),2) +  pow(pIn[i].jet_v().y(),2) + pow(pIn[i].jet_v().z(),2) ) ;
              double c_t = pIn[i].jet_v().z()/mod_jet_v;  // c_t is cos(theta) for the jet_velocity unit vector
              double s_t = std::sqrt( 1.0 - c_t*c_t) ;  // s_t is sin(theta)
                  
              double s_p = pIn[i].jet_v().y()/std::sqrt( pow( pIn[i].jet_v().x() , 2 ) + pow( pIn[i].jet_v().y(), 2 ) ) ;
              double c_p = pIn[i].jet_v().x()/std::sqrt( pow( pIn[i].jet_v().x() , 2 ) + pow( pIn[i].jet_v().y(), 2 ) ) ;

              //JSINFO << BOLDYELLOW << " cos 0 = " << c_t << " sin 0 = " << s_t << " cos p = " << c_p << " sin p = " << s_p ;
              
              //double c_t = new_parent_vz;
              //double s_t = std::sqrt( 1.0 - c_t*c_t) ;
              //    
              //double s_p = new_parent_vy/std::sqrt( pow(new_parent_vx,2) + pow(new_parent_vy,2) ) ;
              //double c_p = new_parent_vx/std::sqrt( pow(new_parent_vx,2) + pow(new_parent_vy,2) ) ;
                  
              // First daughter
              double k_perp1[4];
              k_perp1[0] = 0.0;
              k_perp1[1] = z*(new_parent_px-new_parent_pl*s_t*c_p) + l_perp*std::cos(angle)*c_t*c_p - l_perp*std::sin(angle)*s_p ;
              k_perp1[2] = z*(new_parent_py-new_parent_pl*s_t*s_p) + l_perp*std::cos(angle)*c_t*s_p + l_perp*std::sin(angle)*c_p ;
              k_perp1[3] = z*(new_parent_pz-new_parent_pl*c_t) - l_perp*std::cos(angle)*s_t ;
              double k_perp1_2 = pow(k_perp1[1],2)+pow(k_perp1[2],2)+pow(k_perp1[3],2);
              
              double M=0.0;
              
              if ((pid_a == 4)||(pid_a == 5)) M = pIn[i].restmass();
              
              double energy = ( z*new_parent_nu + (tQd1 + k_perp1_2 + M*M)/(2.0*z*new_parent_nu) )/std::sqrt(2.0) ;
              double plong =  ( z*new_parent_nu - (tQd1 + k_perp1_2 + M*M)/(2.0*z*new_parent_nu) )/std::sqrt(2.0) ;
                  if (energy < 0.0)
                  {
                      JSWARN << " Energy negative after rotation, press 1 and return to continue " ;
                      cin >> blurb ;
                  }
              
              double newp[4];
              newp[0] = energy;
              newp[1] = plong*s_t*c_p + k_perp1[1];
              newp[2] = plong*s_t*s_p + k_perp1[2];
              newp[3] = plong*c_t + k_perp1[3];
                  
              double newx[4];
              //newx[0] = time + deltaT;
              newx[0] = time;
              for (int j=1;j<=3;j++)
              {
                  //newx[j] = pIn[i].x_in().comp(j) + (time + deltaT - pIn[i].x_in().comp(0))*velocity[j]/velocityMod;
                  newx[j] = pIn[i].x_in().comp(j) + (time - pIn[i].x_in().comp(0))*velocity[j]/velocityMod;
              }
                  
              pOut.push_back(Parton(0,pid_a,jet_stat,newp,newx));
              int iout = pOut.size()-1;
                  
              pOut[iout].set_jet_v(velocity_jet); // use initial jet velocity
              pOut[iout].set_mean_form_time();
              double ft = generate_L (pOut[iout].mean_form_time());
              pOut[iout].set_form_time(ft);
              pOut[iout].set_color(d1_col);
              pOut[iout].set_anti_color(d1_acol);
              pOut[iout].set_max_color(max_color);
              pOut[iout].set_min_color(pIn[i].min_color());
              pOut[iout].set_min_anti_color(pIn[i].min_anti_color());
              
              
              // Second daughter
              double k_perp2[4];
              k_perp2[0] = 0.0;
              k_perp2[1] = (1.0-z)*(new_parent_px-new_parent_pl*s_t*c_p) - l_perp*std::cos(angle)*c_t*c_p + l_perp*std::sin(angle)*s_p ;
              k_perp2[2] = (1.0-z)*(new_parent_py-new_parent_pl*s_t*s_p) - l_perp*std::cos(angle)*c_t*s_p - l_perp*std::sin(angle)*c_p ;
              k_perp2[3] = (1.0-z)*(new_parent_pz-new_parent_pl*c_t) + l_perp*std::cos(angle)*s_t ;
              double k_perp2_2 = pow(k_perp2[1],2)+pow(k_perp2[2],2)+pow(k_perp2[3],2);
              
              M=0.0;
              if ((pid_b == 4)||(pid_b == 5)) M = pIn[i].restmass();
              
              energy = ( (1.0-z)*new_parent_nu + (tQd2 + k_perp2_2+M*M)/( 2.0*(1.0-z)*new_parent_nu ) )/std::sqrt(2.0) ;
              plong =  ( (1.0-z)*new_parent_nu - (tQd2 + k_perp2_2+M*M)/( 2.0*(1.0-z)*new_parent_nu ) )/std::sqrt(2.0) ;

              //parent_perp = std::sqrt( pow(pIn[i].p(1),2) + pow(pIn[i].p(2),2) + pow(pIn[i].p(3),2) - pow(pIn[i].pl(),2) );
              //mod_jet_v = std::sqrt( pow(pIn[i].jet_v().x(),2) +  pow(pIn[i].jet_v().y(),2) + pow(pIn[i].jet_v().z(),2) ) ;
    
              if (energy < 0.0)
              {
                  JSWARN << " Energy of 2nd daughter negative after rotation, press 1 and return to continue " ;
                  cin >> blurb ;
              }
              
              newp[0] = energy;
              newp[1] = plong*s_t*c_p + k_perp2[1];
              newp[2] = plong*s_t*s_p + k_perp2[2];
              newp[3] = plong*c_t + k_perp2[3];
              
              //newx[0] = time + deltaT;
              newx[0] = time;
              for (int j=1;j<=3;j++)
              {
                  //newx[j] = pIn[i].x_in().comp(j) + (time + deltaT - pIn[i].x_in().comp(0))*velocity[j]/velocityMod;
                  newx[j] = pIn[i].x_in().comp(j) + (time - pIn[i].x_in().comp(0))*velocity[j]/velocityMod;
              }
	      
              if (iSplit<3) // not a photon
              {
                  pOut.push_back(Parton(0,pid_b,jet_stat,newp,newx));
                  iout = pOut.size()-1;
                  pOut[iout].set_jet_v(velocity_jet); // use initial jet velocity
                  pOut[iout].set_mean_form_time();
                  ft = generate_L (pOut[iout].mean_form_time());
                  pOut[iout].set_form_time(ft);
                  pOut[iout].set_color(d2_col);
                  pOut[iout].set_anti_color(d2_acol);
                  pOut[iout].set_max_color(max_color);
                  pOut[iout].set_min_color(pIn[i].min_color());
                  pOut[iout].set_min_anti_color(pIn[i].min_anti_color());
              }
              else // is a photon
              {
                  pOut.push_back(Photon(0,pid_b,jet_stat,newp,newx));
                  iout = pOut.size()-1;
                  pOut[iout].set_jet_v(velocity_jet); // use initial jet velocity
                  pOut[iout].set_mean_form_time();
                  ft = generate_L (pOut[iout].mean_form_time());
                  pOut[iout].set_form_time(1.0/rounding_error);
                  pOut[iout].set_color(0);
                  pOut[iout].set_anti_color(0);
                  pOut[iout].set_max_color(max_color);
                  pOut[iout].set_min_color(pIn[i].min_color());
                  pOut[iout].set_min_anti_color(pIn[i].min_anti_color());
                  
                  JSINFO << BOLDYELLOW << " A photon was made with px = " << pOut[iout].px() << " and sent to the framework " ;
                  
                  
              }
              
          }
          else
          { // not time to split yet broadening it

              if(broadening_on) 
              {

                  double now_zeta = ( ( time + initRdotV + (time - initR0) )/std::sqrt(2) )*fmToGeVinv;
                  qhat = fncQhat(now_zeta);
                  if (now_temp>0.1)
                  {
                      ehat = 0.0*qhat/4.0/now_temp;
                  }
                  else
                  {
                      ehat = 0.0*qhat/4.0/0.3;
                  }

                 VERBOSE(8) <<BOLDRED<< " between splits broadening qhat = " << qhat << " ehat = " << ehat << " and delT = " << delT ;
                 VERBOSE(8) <<BOLDBLUE<< " zeta at formation = " << zeta << " zeta now = " << now_zeta ;
                  
                  if ((!recoil_on)&&(qhat>0.0))
                  {
                      double kt = 0;
                      if (pIn[i].pid() == 21) // particle is a gluon
                      {
                          kt = generate_kt(qhat*1.414/0.197, delT);
                      }
                      else if ( (pIn[i].pid() < 6)&&( pIn[i].pid() > -6  ) ) // particle is a quark
                      {
                          kt = generate_kt(qhat*1.414/0.197*Cf/Ca, delT);  // scale down q-hat by Cf/Ca
                      }
                      else // a photon, or something else that does not have color
                      {
                          kt = 0;
                      }
                      
                      JSDEBUG << " kt generated = "  << kt << " for qhat = " << qhat*1.414/0.197 << " and delT = " << delT ;
                      
                      double ktx, kty, ktz;
                      ktx=kty=ktz=0.0;
                      double vx = initVx;
                      double vy = initVy;
                      double vz = initVz;

                      
                      bool solved = false;
                      int trials = 0;
                      while ((!solved)&&(trials<1000))
                      {
                          
                          if ((abs(vy)>approx)||(abs(vz)>approx))
                          {
                              ktx = kt*(1 - 2*ZeroOneDistribution(*GetMt19937Generator()) );
                      
                              JSDEBUG << " vx = " << vx << " vy = " << vy << " vz = " << vz ;
                      
                              double rad = sqrt( 4*ktx*ktx*vx*vx*vy*vy - 4*(vy*vy + vz*vz)*(ktx*ktx*(vx*vx+vz*vz) - kt*kt*vz*vz)  );
                          
                              double sol1 = (-2*ktx*vx*vy + rad)/2/(vy*vy + vz*vz);
                              double sol2 = (-2*ktx*vx*vy - rad)/2/(vy*vy + vz*vz);
                          
                              kty = sol1;
                          
                              if ((ktx*ktx + sol1*sol1) > kt*kt  )  kty = sol2;
                      
                              if ((ktx*ktx + kty*kty) < kt*kt  )
                              {
                                  ktz = sqrt( kt*kt - ktx*ktx - kty*kty ) ;
                              
                                  double sign = ZeroOneDistribution(*GetMt19937Generator());
                                  
                                  if (sign>0.5) ktz = -1*ktz;
                              
                                  if (vz!=0) ktz = (-1*ktx*vx - kty*vy )/vz ;
                              
                                  solved = true;
                              }
                              
                          }
                          else
                          {
                              ktx = 0;
                              kty = kt*(1 - 2*ZeroOneDistribution(*GetMt19937Generator()) );
                              double sign = ZeroOneDistribution(*GetMt19937Generator());
                              if (sign>0.5) kty = -1*kty;
                              ktz = sqrt(kt*kt-kty*kty);
                              sign = ZeroOneDistribution(*GetMt19937Generator());
                              if (sign>0.5) ktz = -1*ktz;
                        
                          }
                          trials++;
                      }
                      
                     VERBOSE(8) << MAGENTA << " ktx = " << ktx << " kty = " << kty << " ktz = " << ktz ;

                      double px = pIn[i].px();
                      double py = pIn[i].py();
                      double pz = pIn[i].pz();
                      double energy = pIn[i].e();
                      
                      double p = sqrt(px*px + py*py + pz*pz);

                     VERBOSE(8) << BOLDBLUE << " p before b & d, E = " << energy << " pz = " << pz << " px = " << px << " py = " << py ;


                      
                      px += ktx;
                      py += kty;
                      pz += ktz;
                      
                      double np = sqrt(px*px + py*py + pz*pz);
                      JSDEBUG << " p = " << p << " np = " << np;
 
                      double correction = 0.0;
                      
                      if (np*np > p*p) correction = sqrt(np*np - p*p);
                      
                      px -= correction*vx;
                      py -= correction*vy;
                      pz -= correction*vz;

//                      double nnp = sqrt(px*px + py*py + pz*pz);
                      
/*                      if (abs(nnp-p)>abs(np-p))
                      {

                         VERBOSE(8) << MAGENTA << " negative condition invoked ! " ;

                          px += 2*(np-p)*vx;
                          py += 2*(np-p)*vy;
                          pz += 2*(np-p)*vz;
                      }
*/
                      np = sqrt(px*px + py*py + pz*pz);
                      JSDEBUG << "FINAL p = " << p << " np = " << np;
                      
                      pOut.push_back(pIn[i]);
                      int iout = pOut.size()-1;
                      
                      pOut[iout].reset_p(px,py,pz);
                      
                      double drag = ehat*delT;
                      
                     VERBOSE(8) << MAGENTA << " drag = " << drag << " temperature = " <<  now_temp ;
                      
                      if ((np > drag)&&(energy >drag)&&( energy>sqrt(px*px + py*py + pz*pz) ) )
                      {
                          px -= drag*vx;
                          py -= drag*vy;
                          pz -= drag*vz;
                          energy -= drag;
                          pOut[iout].reset_momentum( px , py , pz , energy );
                      }
                      
                     VERBOSE(8) << BOLDYELLOW << " p after b & d, E = " << energy << " pz = " << pz << " px = " << px << " py = " << py ;


                  }
              
	      } // end if(broadening_on) 
              //pOut.push_back(pIn[i]);
          }
      }
      else
      { // virtuality too low lets broaden it

          if(broadening_on)
          {

              double now_zeta = ( ( time + initRdotV + (time - initR0) )/std::sqrt(2) )*fmToGeVinv;
              //double now_zeta = ( ( time + initRdotV )/std::sqrt(2) )*fmToGeVinv;
              qhat = fncQhat(now_zeta);
              if (now_temp>0.1)
              {
                  ehat = 0.0*qhat/4.0/now_temp;
              }
              else
              {
                  ehat = 0.0*qhat/4.0/0.3;
              }

             VERBOSE(8) <<BOLDRED<< " after splits broadening qhat = " << qhat << " ehat = " << ehat << " and delT = " << delT ;
             VERBOSE(8) <<BOLDBLUE<< " zeta at formation = " << zeta << " zeta now = " << now_zeta ;
              
              //JSINFO << " broadening qhat = " << qhat << " and delT = " << delT ;
              
              if ((!recoil_on)&&(qhat>0.0))
              {
                  double kt = 0.0 ;
                  
                  if (pIn[i].pid() == 21)
                  {
                      kt = generate_kt(qhat*1.414/0.197, delT);
                  }
                  else
                  {
                      kt = generate_kt(qhat*1.414/0.197*Cf/Ca, delT);
                  }
                  
                  JSDEBUG << " kt generated = "  << kt << " for qhat = " << qhat*1.414/0.197 << " and delT = " << delT ;
                  
                  double ktx, kty, ktz;
                  ktx=kty=ktz=0;
                  double vx = initVx;
                  double vy = initVy;
                  double vz = initVz;

                  bool solved = false;
                  
                  int trials = 0;
                  while ((!solved)&&(trials<1000))
                  {
                      if ((abs(vy)>approx)||(abs(vz)>approx))
                      {
                      
                          ktx = kt*(1 - 2*ZeroOneDistribution(*GetMt19937Generator()) );
                      
                          JSDEBUG << " vx = " << vx << " vy = " << vy << " vz = " << vz ;
                          
                          double rad = sqrt( 4*ktx*ktx*vx*vx*vy*vy - 4*(vy*vy + vz*vz)*(ktx*ktx*(vx*vx+vz*vz) - kt*kt*vz*vz)  );
                      
                          double sol1 = (-2*ktx*vx*vy + rad)/2/(vy*vy + vz*vz);
                          double sol2 = (-2*ktx*vx*vy - rad)/2/(vy*vy + vz*vz);
                      
                          kty = sol1;
                      
                          if ((ktx*ktx + sol1*sol1) > kt*kt  )  kty = sol2;
                      
                          if ((ktx*ktx + kty*kty) < kt*kt  )
                          {
                              ktz = sqrt( kt*kt - ktx*ktx - kty*kty ) ;
                          
                              double sign = ZeroOneDistribution(*GetMt19937Generator());
                              if (sign>0.5) ktz = -1*ktz;
                          
                          
                              if (vz!=0) ktz = (-1*ktx*vx - kty*vy )/vz ;
                          
                              solved = true;
                          }
                      }
                      else
                      {
                          ktx = 0;
                          kty = kt*(1 - 2*ZeroOneDistribution(*GetMt19937Generator()) );
                          double sign = ZeroOneDistribution(*GetMt19937Generator());
                          if (sign>0.5) kty = -1*kty;
                          ktz = sqrt(kt*kt-kty*kty);
                          sign = ZeroOneDistribution(*GetMt19937Generator());
                          if (sign>0.5) ktz = -1*ktz;
                          
                      }
                      trials++;
                  }
                  
                 VERBOSE(8) << MAGENTA << " ktx = " << ktx << " kty = " << kty << " ktz = " << ktz ;

                  double px = pIn[i].px();
                  double py = pIn[i].py();
                  double pz = pIn[i].pz();
                  double energy = pIn[i].e();

                  
                  double p = sqrt(px*px + py*py + pz*pz);
                  
                 VERBOSE(8) << BOLDBLUE << " p before b & d, E = " << energy << " pz = " << pz << " px = " << px << " py = " << py ;
 
                  px += ktx;
                  py += kty;
                  pz += ktz;
                  
                  double np = sqrt(px*px + py*py + pz*pz);
                  JSDEBUG << " p = " << p << " np = " << np;
                  double correction = 0.0;
                  if (np*np > p*p) correction = sqrt(np*np - p*p);
                  
                  px -= correction*vx;
                  py -= correction*vy;
                  pz -= correction*vz;
                  
                //  double nnp = sqrt(px*px + py*py + pz*pz);
                  
/*                  if (abs(nnp-p)>abs(np-p))
                  {
                      px += 2*(np-p)*vx;
                      py += 2*(np-p)*vy;
                      pz += 2*(np-p)*vz;
                  }
*/
                  np = sqrt(px*px + py*py + pz*pz);
                  JSDEBUG << "FINAL p = " << p << " nnnp = " << np;

                  pOut.push_back(pIn[i]);
                  int iout = pOut.size()-1;
                  
                  pOut[iout].reset_p(px,py,pz);
                  np = sqrt(px*px + py*py + pz*pz );
                  double drag = ehat*delT;
           
                 VERBOSE(8) << MAGENTA << " drag = " << drag << " temperature = " <<  now_temp ;

                  if ((np > drag)&&(energy >drag)&&( energy>sqrt(px*px + py*py + pz*pz) ) )
                  {
                      px -= drag*vx;
                      py -= drag*vy;
                      pz -= drag*vz;
                      energy -= drag;
                      pOut[iout].reset_momentum( px, py , pz , energy );
                  }
                
                 VERBOSE(8) << BOLDYELLOW << " p after b & d, E = " << energy << " pz = " << pz << " px = " << px << " py = " << py ;

              }
              
              
              //pOut.push_back(pIn[i]);
          }
      }
      
  } // particle loop
      
}

double Matter::generate_kt(double local_qhat, double dzeta)
{
    double width = local_qhat*dzeta;
    
    double x,r;
    
    r = ZeroOneDistribution(*GetMt19937Generator());
    x = -log(1.0-r);
    if (x<0.0)  throw std::runtime_error(" k_t^2 < 0.0 ");
    double kt = sqrt(x * local_qhat * dzeta  );
    
    return(kt);
    
}
    
    


double Matter::generate_angle()
{
  double ang, r, blurb;
    
  // r = double(random())/ (maxN );
  r = ZeroOneDistribution(*GetMt19937Generator());
  //    r = mtrand1();
    
  //    cout << " r = " << r << endl;
  //    cin >> blurb ;
    
  ang = r*2.0*pi ;
    
  return(ang);
    
}


double Matter::generate_vac_t(int p_id, double nu, double t0, double t, double loc_a, int is)
{
  double r,z,ratio,diff,scale,t_low, t_hi, t_mid, numer, denom, test ;
    
    
  // r = double(random())/ (maxN );
  r = ZeroOneDistribution(*GetMt19937Generator());
  //        r = mtrand1();
    
    
  if ((r>=1.0)||(r<=0.0))
    {
      throw std::runtime_error("error in random number in t *GetMt19937Generator()");
    }
    
  ratio = 1.0 ;
    
  diff = (ratio - r)/r;
    
  t_low = 2.0*t0;
    
  t_hi = t;
    
  //    cout << " in gen_vac_t : t_low , t_hi = " << t_low << "  " << t_hi << endl;
  //    cin >> test ;
    
  if (p_id==gid)
  {
      numer = sudakov_Pgg(t0,t,loc_a,nu)*std::pow(sudakov_Pqq(t0,t,loc_a,nu),nf);
        
      if ((is!=1)&&(is!=2)) // there is almost no use of is = 2, the `is' in this function is redundant, just for consistency
      {
          throw std::runtime_error(" error in isp ");
      }
  }
  else
  {
      if ( (is!=0)&&(is!=3) ) // there is almost no use of is = 3, the `is' in this function is redundant, just for consistency
      {
          throw std::runtime_error("error in isp in quark split");
      }
      numer = sudakov_Pqg(t0,t,loc_a,nu)*sudakov_Pqp(t0,t,loc_a,nu);
  }
    
  t_mid = t_low;
    
  if (numer>r)
    {
      // cout << " numer > r, i.e. ; " << numer << " > " << r << endl ;        
      return(t_mid) ;
    }
    
  //t_mid = (t_low+t_hi)/2.0 ;
    
    
  scale = t0;
    
  //   cout << " s_approx, s_error = " << s_approx << "  " << s_error << endl;
    
  do {

      t_mid = (t_low + t_hi)/2.0;

      if (p_id==gid)
      {
          denom = sudakov_Pgg(t0, t_mid, loc_a, nu)*std::pow(sudakov_Pqq(t0, t_mid, loc_a, nu),nf);
          if ((is!=1)&&(is!=2))
          {
              throw std::runtime_error(" error in isp numerator");
          }
          
      }
      else
      {
          if ( (is!=0)&&(is!=3) )
          {
              throw std::runtime_error(" error in isp in quark split numerator  ");
              
          }
          denom = sudakov_Pqg(t0, t_mid, loc_a, nu)*sudakov_Pqp(t0, t_mid, loc_a, nu);
          
      }
        
      ratio = numer/denom ;
        
      diff = (ratio - r)/r ;
        
      //       cout << "num, den, r = " << numer << " "<< denom << " " << r << " " << endl;
      //       cout << "diff, t_mid = " << diff << " " << t_mid << endl;
      //       cout << " t_low, t_hi = " << t_low << "  " << t_hi << endl;
      //       cin >> test ;
                
      if (diff<0.0)
        {
	  t_low = t_mid ;
	  //t_mid = (t_low + t_hi)/2.0;
        }
      else
        {
	  t_hi = t_mid ;
	  //t_mid = (t_low + t_hi)/2.0;
        }
        
  } while ((abs(diff)>s_approx)&&(abs(t_hi-t_low)/t_hi>s_error));
    
  return(t_mid);
}

double Matter::generate_vac_t_w_M(int p_id, double M, double nu, double t0, double t, double loc_a, int is)
{
  double r,z,ratio,diff,scale,t_low_M0,t_low_MM,t_low_00, t_hi_M0, t_hi_MM, t_hi_00, t_mid_M0, t_mid_MM, t_mid_00, numer, denom, test ;
  double M_charm=1.5;//InternalHelperPythia.particleData.m0(4);
  double M_bottom=4.8;//InternalHelperPythia.particleData.m0(5);  
    
  // r = double(random())/ (maxN );
  r = ZeroOneDistribution(*GetMt19937Generator());
  //        r = mtrand1();
    
    
  if ((r>=1.0)||(r<=0.0))
    {
      throw std::runtime_error("error in random number in t *GetMt19937Generator()");
    }
    
  ratio = 1.0 ;
    
  diff = (ratio - r)/r;
    
//    if (t0<M*M) t0 = M*M;
    
  t_low_M0 = t0*(1.0 + std::sqrt(1.0 + 2.0*M*M/t0));
  t_low_MM = 2.0*(M*M + t0);
  t_low_00 = 2.0*t0;
    

  t_hi_M0 = t;
  t_hi_MM = t;
  t_hi_00 = t;
    
    //JSINFO << MAGENTA << " in gen_vac_t_w_M : t_low , t_hi = " << t_low << "  " << t_hi ;
    //cin >> test ;
    
  if (p_id==gid)
    {
      numer = sudakov_Pgg(t0,t,loc_a,nu)*std::pow(sudakov_Pqq(t0,t,loc_a,nu),3.0)*sudakov_Pqq_w_M_vac_only(M_charm,t0,t,loc_a,nu)*sudakov_Pqq_w_M_vac_only(M_bottom,t0,t,loc_a,nu);
        
      if ((is!=1)&&(is!=2))
        {
	  throw std::runtime_error(" error in isp ");
        }
    }
  else
    {
      if (is!=0)
        {
	  throw std::runtime_error("error in isp in quark split");            
        }
      if (((int) std::abs((double) p_id))==4 || ((int) std::abs((double) p_id))==5)
        {
         numer = sudakov_Pqg_w_M(M,t0,t,loc_a,nu);
        }
      else
        {
         numer = sudakov_Pqg(t0,t,loc_a,nu);
        }
    }
    
  t_mid_M0 = t_low_M0;
  t_mid_MM = t_low_MM;
  t_mid_00 = t_low_00;
    
  if (numer>r)
    {
      // cout << " numer > r, i.e. ; " << numer << " > " << r << endl ;
      return(t_mid_M0) ;
    }
    
  //t_mid = (t_low+t_hi)/2.0 ;
    
    
  scale = t0;
    
  //   cout << " s_approx, s_error = " << s_approx << "  " << s_error << endl;
    
  do {

      t_mid_M0 = (t_low_M0 + t_hi_M0)/2.0;
      t_mid_MM = (t_low_MM + t_hi_MM)/2.0;
      t_mid_00 = (t_low_00 + t_hi_00)/2.0;

      if (p_id==gid)
        {
	  denom = sudakov_Pgg(t0, t_mid_00, loc_a, nu)*std::pow(sudakov_Pqq(t0, t_mid_00, loc_a, nu),3.0)*sudakov_Pqq_w_M_vac_only(M_charm, t0, t_mid_MM, loc_a, nu)*sudakov_Pqq_w_M_vac_only(M_bottom, t0, t_mid_MM, loc_a, nu);
            
	  if ((is!=1)&&(is!=2))
            {
	      throw std::runtime_error(" error in isp numerator");             
            }
        }
      else
        {
	  if (is!=0)
            {
	      throw std::runtime_error(" error in isp in quark split numerator  ");
            }
          if (((int) std::abs((double) p_id))==4 || ((int) std::abs((double) p_id))==5)
            {
             denom = sudakov_Pqg_w_M(M, t0, t_mid_M0, loc_a, nu);
            }
          else
            {
	     denom = sudakov_Pqg(t0, t_mid_00, loc_a, nu);
            }
        }
        
      ratio = numer/denom ;
        
      diff = (ratio - r)/r ;
        
      //       cout << "num, den, r = " << numer << " "<< denom << " " << r << " " << endl;
      //       cout << "diff, t_mid = " << diff << " " << t_mid << endl;
      //       cout << " t_low, t_hi = " << t_low << "  " << t_hi << endl;
      //       cin >> test ;
                
      if (diff<0.0)
        {
	  t_low_M0 = t_mid_M0 ;
	  t_low_MM = t_mid_MM ;
	  t_low_00 = t_mid_00 ;

	  //t_mid = (t_low + t_hi)/2.0;
        }
      else
        {
	  t_hi_M0 = t_mid_M0 ;
	  t_hi_MM = t_mid_MM ;
	  t_hi_00 = t_mid_00 ;

	  //t_mid = (t_low + t_hi)/2.0;
        }
        
  } while ((abs(diff)>s_approx)&&(abs(t_hi_M0-t_low_M0)/t_hi_M0>s_error));
    
  return(t_mid_M0);
}

/*
 
 
 
  New function
 
 
 
*/



double  Matter::generate_vac_z(int p_id, double t0, double t, double loc_b, double nu, int is )
{
  double r,z, ratio,diff,e,numer1, numer2, numer, denom, z_low, z_hi, z_mid, test;

  r = ZeroOneDistribution(*GetMt19937Generator());
    
  if ((r>1)||(r<0))
    {
      throw std::runtime_error(" error in random number in z *GetMt19937Generator()");
    }
    
  ratio = 1.0 ;
    
  diff = (ratio - r)/r;
    
  e = t0/t;
    
  if (e>0.5)
    {
      throw std::runtime_error(" error in epsilon");
    }
    
  z_low = e ;
    
  z_hi = double(1.0) - e ;
    
  if (p_id==gid)
  {
      if (is==1)
      {
          denom = P_z_gg_int(z_low,z_hi, loc_b, t, 2.0*nu/t, nu );
      }
      else
      {
          denom = P_z_qq_int(z_low,z_hi,loc_b,t,2.0*nu/t, nu);
      }
      
  }
  else if ( (p_id!=gid)&&(is==0) )
  {
      denom = P_z_qg_int(z_low,z_hi, loc_b, t, 2.0*nu/t , nu);
  }
  else if ( (p_id!=gid)&&(is==3) )
  {
      denom = P_z_qp_int(z_low,z_hi, loc_b, t, 2.0*nu/t , nu);
  }
  else
  {
      throw std::runtime_error(" I do not understand your combination of particle id and split id in denominator of generate_z ") ;
  }
    
  //z_mid = (z_low + z_hi)/2.0 ;
  
  int itcounter=0;
  // cout << " generate_vac_z called with p_id = " << p_id << " t0 = " << t0 << " t = " << t << " loc_b=" << loc_b<< " nu = " <<  nu << " is = " << is << endl;

  do { // Getting stuck in here for some reason
    if ( itcounter++ > 10000 ) {
    	cout << " in here" << " abs(diff) = " << abs(diff) << "  approx = " << approx << "  r = " << r << "  zmid = " << z_mid << "  denom = " << denom << "  numer = " << numer << "  e = " << e << "   " << numer/denom << endl;
        throw std::runtime_error("Stuck in endless loop") ;
    }

    z_mid = (z_low + z_hi)/2.0 ;

    if (p_id==gid)
    {
      if (is==1)
      {
          numer = P_z_gg_int(e, z_mid, loc_b, t, 2.0*nu/t , nu );
      }
      else
      {
          numer = P_z_qq_int(e, z_mid, loc_b, t, 2.0*nu/t , nu);
      }
    }
    else if ( (p_id!=gid)&&(is==0) )
    {
      numer = P_z_qg_int(e, z_mid, loc_b, t, 2.0*nu/t , nu );
    }
    else if ( (p_id!=gid)&&(is==3) )
    {
      numer = P_z_qp_int(e, z_mid, loc_b, t, 2.0*nu/t , nu);
    }
    else
    {
        throw std::runtime_error(" I do not understand your combination of particle id and split id in numerator of generate_z ") ;
    }
    ratio = numer/denom ;      
    diff = (ratio - r)/r ;
        
    if (diff>0.0)
    {
        z_hi = z_mid;
        //z_mid = (z_low + z_hi)/2.0;
    }
    else
    {
        z_low = z_mid;
        //z_mid = (z_low + z_hi)/2.0 ;
    }
        
  } while ( (abs(diff)>approx)&&(abs(z_hi-z_low)/z_hi>s_error) );
    
  return(z_mid);
}	


double  Matter::generate_vac_z_w_M(int p_id, double M, double t0, double t, double loc_b, double nu, int is )
{
  double r,z, ratio,diff,e,numer1, numer2, numer, denom, z_low_00,z_low_M0,z_low_MM, z_hi_00, z_hi_M0, z_hi_MM, z_mid_00,z_mid_M0,z_mid_MM, test;

  r = ZeroOneDistribution(*GetMt19937Generator());
    
//    if (t0<M*M) t0 = M*M;
    
  if ((r>1)||(r<0))
    {
      throw std::runtime_error(" error in random number in z *GetMt19937Generator()");
    }
    
  ratio = 1.0 ;
    
  diff = (ratio - r)/r;
    
  e = t0/t;

  double M2 = M*M;
    
  if (e>0.5)
    {
      throw std::runtime_error(" error in epsilon");
    }
    
  if (((int) std::abs((double) p_id))==4 || ((int) std::abs((double) p_id))==5)
    {
      z_low_M0 = e + M2/(t+M2);
    }
  else
    {
      z_low_00 = e;
      z_low_MM = e+M2/t;
    }

  z_hi_00 = double(1.0) - e;
  z_hi_M0 = double(1.0) - e;
  z_hi_MM = double(1.0) - e  - M2/t;
    
  if (p_id==gid)
    {
      if (is==1)
        {
	  denom = P_z_gg_int(z_low_00,z_hi_00, loc_b, t, 2.0*nu/t, nu );
        }
      else
        {
          if (is==3) //THIS IS FOR GLUON-> HEAVY Q-QBAR
            {
             denom = P_z_qq_int_w_M_vac_only(M,z_low_MM,z_hi_MM,loc_b,t,2.0*nu/t, nu);
            }
          else
            {
	     denom = P_z_qq_int(z_low_00,z_hi_00,loc_b,t,2.0*nu/t, nu);
            }
        }
    }
  else
    {
      if (((int) std::abs((double) p_id))==4 || ((int) std::abs((double) p_id))==5)
        {
         denom = P_z_qg_int_w_M(M, z_low_M0, z_hi_M0, loc_b, t, 2.0*nu/t , nu);
        }
      else
        {
         denom = P_z_qg_int(z_low_00,z_hi_00, loc_b, t, 2.0*nu/t , nu);
        }
    }
    
    
  //z_mid = (z_low + z_hi)/2.0 ;
  
  int itcounter=0;
  //JSINFO << BOLDYELLOW << " generate_vac_z called with p_id = " << p_id << " t0 = " << t0 << " t = " << t << " loc_b=" << loc_b<< " nu = " <<  nu << " is = " << is ;

  do
  { // Getting stuck in here for some reason

    if ( itcounter++ > 10000 )
    {
    	cout << " in here" << " abs(diff) = " << abs(diff) << "  approx = " << approx << "  r = " << r << "  zmid = " << z_mid_M0 << "  denom = " << denom << "  numer = " << numer << "  e = " << e << "   " << numer/denom << endl;
        throw std::runtime_error("Stuck in endless loop") ;
        
    }


    z_mid_00 = (z_low_00 + z_hi_00)/2.0 ;
    z_mid_M0 = (z_low_M0 + z_hi_M0)/2.0 ;
    z_mid_MM = (z_low_MM + z_hi_MM)/2.0 ;

    if (p_id==gid)
    {
      if (is==1)
      {
          numer = P_z_gg_int(z_low_00, z_mid_00, loc_b, t, 2.0*nu/t , nu );
      }
      else
      {
        if (is==3)
        {
          numer = P_z_qq_int_w_M_vac_only(M, z_low_MM, z_mid_MM, loc_b, t, 2.0*nu/t , nu); 
        }
        else
        {
          numer = P_z_qq_int(z_low_00, z_mid_00, loc_b, t, 2.0*nu/t , nu);
        }
      }
    }
    else
    {
 
      if (((int) std::abs((double) p_id))==4 || ((int) std::abs((double) p_id))==5)
      {
        numer = P_z_qg_int_w_M(M, z_low_M0, z_mid_M0, loc_b, t, 2.0*nu/t , nu );
      }
      else
      {
        numer = P_z_qg_int(z_low_M0, z_mid_M0, loc_b, t, 2.0*nu/t , nu );
      }
    }
      ratio = numer/denom ;
      diff = (ratio - r)/r ;
    // cout << "num, den, r, diff = " << numer << " "<< denom << " " << r << " " << endl;	
    // cout << " diff, z_mid = " << diff << " " << z_mid << endl ;		
    //		cin >> test ;
      
        
    if (diff>0.0)
    {
        z_hi_M0 = z_mid_M0;
        //z_mid = (z_low + z_hi)/2.0;
    }
    else
    {
        z_low_M0 = z_mid_M0;
        //z_mid = (z_low + z_hi)/2.0 ;
    }
        
  } while ((abs(diff)>approx)&&(abs(z_hi_M0-z_low_M0)/z_hi_M0>s_error) );
    
  return(z_mid_M0);
}




double Matter::generate_L(double form_time)
{
  double r, x_low, x_high , x , diff , span, val, arg, norm ;
    
  // r = double(random())/ (maxN );
  r = ZeroOneDistribution(*GetMt19937Generator());
  //    r = mtrand1();
    
  if ((r>1)||(r<0))
    {
      throw std::runtime_error(" error in random number in z *GetMt19937Generator()");
    }
    
  x_low = 0;
    
  x_high = 8.0*form_time;
  // the value of x_high is slightly arbitrary, the erf function is more or less zero at this distance.
  // picking 10*form_time will not lead to any different results
    
  x = (x_low + x_high)/2.0;
    
  span = (x_high - x_low)/x_high ;
    
  arg = x/form_time/std::sqrt(pi) ;
    
  val = std::erf(arg);
    
  diff = std::abs(val - r) ;
    
  while( (diff>approx)&&(span>error) )
    {
      if ((val-r)>0.0)
        {
            x_high = x ;
        }
      else
        {
            x_low = x ;
        }
        
      x = (x_low + x_high)/2.0;
        
      arg = x/form_time/std::sqrt(pi) ;
        
      val = std::erf(arg);
        
      diff = std::abs(val - r) ;
        
      span = (x_high - x_low)/x_high ;
        
    }
    
  //	cout << " random number for dist = " << r << " distance generated = " << x << endl;
    
    
  return(x);
    

}


double Matter::sudakov_Pgg(double g0, double g1, double loc_c, double E)
{
  double sud,g;
  int blurb;
    
  sud = 1.0 ;
    
  if (g1<2.0*g0)
    {
      cerr << " warning: the lower limit of the sudakov > 1/2 upper limit, returning 1 " << endl;
      cerr << " in sudakov_P glue glue, g0, g1 = " << g0 << "  " << g1 << endl;
      throw std::runtime_error(" warning: the lower limit of the sudakov > 1/2 upper limit, returning 1");
        
      return(sud) ;
        
    }
  g = 2.0*g0;
    
  if (g1>g)
    {
        
      sud = exp( -1.0*(Ca/2.0/pi)*sud_val_GG(g0,g,g1,loc_c, E) );
        
    }
  return(sud);
    
}


double Matter::sud_val_GG(double h0, double h1, double h2, double loc_d, double E1)
{
  double val, h , intg, hL, hR, diff, intg_L, intg_R, t_form, span;
    
  val = 0.0;
    
  h = (h1+h2)/2.0 ;
    
  span = (h2 - h1)/h2;
    
  t_form = 2.0*E1/h;
    
  val = alpha_s(h)*sud_z_GG(h0,h, loc_d, t_form,E1);
    
  intg = val*(h2-h1);
    
  hL = (h1 + h)/2.0 ;
    
  t_form = 2.0*E1/hL;
    
  intg_L = alpha_s(hL)*sud_z_GG(h0,hL,loc_d,t_form,E1)*(h - h1) ;
    
  hR = (h + h2)/2.0 ;
    
  t_form = 2.0*E1/hR;
    
  intg_R = alpha_s(hR)*sud_z_GG(h0,hR,loc_d,t_form,E1)*(h2 - h) ;
    
  diff = std::abs( (intg_L + intg_R - intg)/intg ) ;
    
  //	cout << " iline, gap, diff = " << i_line << " " << h2 << " " << h1 << "  " << diff << endl ;
  //	cout << " intg, Left , right = " << intg << " " << intg_L << "  " << intg_R << endl;
    
  if ( (diff>approx)||(span>error) )
    {
      intg = sud_val_GG(h0,h1,h,loc_d,E1) + sud_val_GG(h0,h,h2,loc_d,E1);
    }
    
  //	cout << " returning with intg = " << intg << endl;
    
  return(intg);
    
}


double Matter::sud_z_GG(double cg, double cg1, double loc_e , double l_fac, double E2)
{
    
  double t2,t3,t7,t11,t12,t15,t21,t25, q2, q3, q8, q12, qL, tau, res, z_min, limit_factor, lz, uz, mz, m_fac;
  double t_q1,t_q3,t_q4,t_q6,t_q8,t_q9,t_q12,q_q1,q_q4,q_q6,q_q9,q_q11;
    
  z_min = std::sqrt(2)*E_minimum/E2;
    
  if (cg1<2.0*cg)
    {
        
      //        cout << " returning with cg, cg1 = " << cg << "   " <<  cg1 << "    " << E_minimum << "  " << E2 << endl ;
      return(0.0);
    };
    
  t2 = std::pow(cg1, 2);
  t3 = t2 * cg1;
  t7 = std::log(cg);
  t11 = std::abs(cg - cg1);
  t12 = std::log(t11);
  t15 = std::pow(cg, 2);
  t21 = t2 * t2;
  t25 = -(5.0 * t3 - 12.0 * cg * t2 + 6.0 * t7 * t3 - 6.0 * t12 * t3 - 4.0 * t15 * cg + 6.0 * t15 * cg1) / t21 / 3.0;
    
  res = t25;
    
  limit_factor = 2.0*std::sqrt(2.0)*cg1/E2/0.1 ;
    
  if (limit_factor<0.0) {
    cerr << " error in z limit factor for medium calculation in sud-z-gg = " << limit_factor << endl ;
    throw std::runtime_error("error in z limit factor for medium calculation in sud-z-gg");
  }
    
  q2 = 1.0 / cg1;
  q3 = cg * q2;
  q8 = 1.0 - q3;
  q12 = (2.0 - 4.0 * q3 + 2.0 / cg * cg1 - 2.0 / q8) * q2;
    
  if (q12<0.0)
    {
      cerr << "ERROR: medium contribution negative in sud_z_GG: q12 = " << q12 << endl;
      cerr << "cg, cg1 = " << cg << "  " << cg1 << endl;
      cerr << " t25 = " << t25 << endl;
      throw std::runtime_error("ERROR: medium contribution negative in sud_z_GG");
    }
    
  tau = l_fac ;
    
  if ((length - loc_e) < tau) tau = (length - loc_e);
    
  if (loc_e > length) tau = 0.0 ;
    
  m_fac = 1.0;
    
  // SC
  //qL = m_fac*qhat*0.6*tau*profile(loc_e+tau) ;
  if(tau<rounding_error) {
      qL = 0.0;
  } else {
      qhat = fncAvrQhat(loc_e,tau);
      qL = qhat*0.6*tau;
  }

    
    
  res = t25 + 2.0*qL*q12/cg1 ;
  //        }
  //        else{
  //            cout << " z trap for medium enabled in sud-val-z " << endl ;
  //        }
  //    }
    
  return(res);
    
}

double Matter::P_z_gg_int(double cg, double cg1, double loc_e, double cg3, double l_fac,double E2)
{
    
  double t3,t4,t5,t10,t11,t12,t15, t9, qL, tau, res, limit_factor, lz, uz, m_fac;
    
  if ((cg< cg1/(2.0*E2*E2/cg1+1.0) )) cg = cg1/( 2.0*E2*E2/cg1 + 1.0 );
    
  t3 = std::log((1.0 - cg1));
  t4 = std::log(cg1);
  t5 = std::pow(cg1,2);
  t10 = std::log((1.0 - cg));
  t11 = std::log(cg);
  t12 = std::pow(cg,2);
  t15 = -(2.0 * cg1) - t3 + t4 - (2.0/3.0) * t5 * cg1 + t5 + (2.0 * cg) + t10 - t11 + (2.0/3.0) * t12 * cg - t12;
    
  res = t15 ;
    
  limit_factor = 2.0*std::sqrt(2.0)*cg3/E2/0.1 ;
    
    
  if (limit_factor<0.0) {
    cerr << " error in z limit factor for medium calculation = " << limit_factor << endl ;
    throw std::runtime_error(" error in z limit factor for medium calculation");
  }
    
    
  tau = l_fac;
    
  if ((length - loc_e) < tau) tau = (length - loc_e);
    
  //             if ((length - loc_e) < tau) tau = 0;
    
  if (loc_e > length) tau = 0.0 ;
    
  m_fac = 1.0;
    
  //            if ((qhat*tau<1.0)&&(in_vac==false)) m_fac = 1.0/qhat/tau;
    
  // SC
  //qL = m_fac*qhat*0.6*tau*profile(loc_e + tau) ;
  if(tau<rounding_error) {
      qL = 0.0;
  } else {
      qhat = fncAvrQhat(loc_e,tau);
      qL = qhat*0.6*tau;
  }

  t9 = 2.0*cg1 - 1.0 / (-1.0 + cg1) - 1.0 / cg1 - 2.0 * cg + 1.0 / (-1.0 + cg) + 1.0 / cg;
    
  if (t9<0.0)
    {
      cerr << "ERROR: medium contribution negative in P_z_gg_int : t9 = " << t9 << endl;
        
      cerr << " cg, cg1 = " << cg << "  " << cg1 << endl;
      throw std::runtime_error("ERROR: medium contribution negative in P_z_gg_int");
    }
    
  res = t15 + 2.0*t9*qL/cg3 ;
        
  return(res);
    
}




double Matter::sudakov_Pqq(double q0, double q1, double loc_c, double E)
{
  double sud,q;
    
  sud = 1.0 ;
    
    
    
  if (q1<2.0*q0)
    //	if (g1<g0)
    {
      JSWARN << " warning: the lower limit of the sudakov > 1/2 upper limit, returning 1 ";
      JSWARN << " in sudakov_Pquark quark, q0, q1 = " << q0 << "  " << q1;
      return(sud) ;
    }
  q = 2.0*q0;
    
  //	g = g0 ;
    
  sud = exp( -1.0*(Tf/2.0/pi)*sud_val_QQ(q0,q,q1,loc_c, E) );
    
  return(sud);
    
}

double Matter::sudakov_Pqq_w_M_vac_only(double M, double q0, double q1, double loc_c, double E)
{
  double sud,q;
    
  sud = 1.0;
       
  if (q1<2.0*(q0+M*M))
    {
      JSWARN << " warning: the lower limit of the sudakov > 1/2 upper limit, returning 1 ";
      JSWARN << " in sudakov_Pquark quark, q0, q1 = " << q0 << "  " << q1;
      return(sud) ;
    }
  else
    {
      q = 2.0*(q0+M*M);
      sud = exp( -1.0*(Tf/2.0/pi)*sud_val_QQ_w_M_vac_only(M,q0,q,q1,loc_c, E) );
      return(sud);
    }
}




double Matter::sud_val_QQ(double h0, double h1, double h2, double loc_d, double E1)
{
  double val, h , intg, hL, hR, diff, intg_L, intg_R, t_form, span;
    
  h = (h1+h2)/2.0 ;
    
  span = (h2 - h1)/h2;
    
  t_form = 2.0*E1/h;
    
  val = alpha_s(h)*sud_z_QQ(h0,h, loc_d, t_form,E1);
    
  intg = val*(h2-h1);
    
  hL = (h1 + h)/2.0 ;
    
  t_form = 2.0*E1/hL;
    
  intg_L = alpha_s(hL)*sud_z_QQ(h0,hL,loc_d,t_form,E1)*(h - h1) ;
    
  hR = (h + h2)/2.0 ;
    
  t_form = 2.0*E1/hR;
    
  intg_R = alpha_s(hR)*sud_z_QQ(h0,hR,loc_d,t_form,E1)*(h2 - h) ;
    
  diff = std::abs( (intg_L + intg_R - intg)/intg ) ;
    
  //	cout << " iline, gap, diff = " << i_line << " " << h2 << " " << h1 << "  " << diff << endl ;
  //	cout << " intg, Left , right = " << intg << " " << intg_L << "  " << intg_R << endl;
    
  if ( (diff>approx)||(span>error) )
    {
      intg = sud_val_QQ(h0,h1,h,loc_d,E1) + sud_val_QQ(h0,h,h2,loc_d,E1);
    }
    
  //	cout << " returning with intg = " << intg << endl;
    
  return(intg);
    
}

double Matter::sud_val_QQ_w_M_vac_only(double M, double h0, double h1, double h2, double loc_d, double E1)
{
  double val, h , intg, hL, hR, diff, intg_L, intg_R, t_form, span;
    
  h = (h1+h2)/2.0 ;
    
  span = (h2 - h1)/h2;
    
  t_form = 2.0*E1/h;
    
  val = alpha_s(h)*sud_z_QQ(h0,h, loc_d, t_form,E1);
    
  intg = val*(h2-h1);
    
  hL = (h1 + h)/2.0 ;
    
  t_form = 2.0*E1/hL;
    
  intg_L = alpha_s(hL)*sud_z_QQ_w_M_vac_only(M,h0,hL,loc_d,t_form,E1)*(h - h1) ;
    
  hR = (h + h2)/2.0 ;
    
  t_form = 2.0*E1/hR;
    
  intg_R = alpha_s(hR)*sud_z_QQ_w_M_vac_only(M,h0,hR,loc_d,t_form,E1)*(h2 - h) ;
    
  diff = std::abs( (intg_L + intg_R - intg)/intg ) ;
    
  //	cout << " iline, gap, diff = " << i_line << " " << h2 << " " << h1 << "  " << diff << endl ;
  //	cout << " intg, Left , right = " << intg << " " << intg_L << "  " << intg_R << endl;
    
  if ( (diff>approx)||(span>error) )
    {
      intg = sud_val_QQ_w_M_vac_only(M,h0,h1,h,loc_d,E1) + sud_val_QQ_w_M_vac_only(M,h0,h,h2,loc_d,E1);
    }
    
  //	cout << " returning with intg = " << intg << endl;
    
  return(intg);
    
}




double Matter::sud_z_QQ(double cg, double cg1, double loc_e , double l_fac, double E2)
{
    
  double t2,t4,t5,t7,t9,t14, q2, q3, q5, q6, q8, q15, qL, tau, res, z_min;
    
    
  z_min = std::sqrt(2)*E_minimum/E2;
    
  //    if (cg<cg1*z_min) cg = cg1*z_min;
    
    
    
  //    if ((cg< cg1/(2.0*E2*E2/cg1+1.0) )) cg = cg1/( 2.0*E2*E2/cg1 + 1.0 );
    
  if (cg1<2.0*cg)
    {
        
      //        cout << " returning with cg, cg1 = " << cg << "   " <<  cg1 << "    " << E_minimum << "  " << E2 << endl ;
      return(0.0);
    };
    
    
    
  t2 = 1.0 / cg1;
  t4 = 1.0 - cg * t2;
  t5 = t4 * t4;
  t7 = std::pow( cg, 2.0);
  t9 = std::pow(cg1, 2.0);
  t14 = ((t5 * t4) - t7 * cg / t9 / cg1) * t2 / 3.0;
    
  //	return(t25);
    
  q2 = 1.0 / cg1;
  q3 = (cg * q2);
  q5 = 1.0 - q3;
  q6 = std::log(std::abs(q5));
  q8 = std::log(q3);
  /*	q10 = std::log(q3);
	q12 = std::log(q5);*/
  q15 = (-1.0 + (2.0 * q3) + q6 - q8) * q2;
    
    
    
  if (q15<0.0)
    {
      cerr << "ERROR: medium contribution negative in sud_z_QQ: q15 = " << q15 << endl;
      cerr << "cg, cg1 = " << cg << "  " << cg1 << endl;
      cerr << " t14 = " << t14 << endl;
      throw std::runtime_error("ERROR: medium contribution negative in sud_z_QQ");
    }
    
  tau = l_fac ;
    
  if ((length - loc_e) < tau) tau = (length - loc_e);
    
  if (loc_e > length) tau = 0.0 ;

  // SC  
  //qL = qhat*0.6*tau*profile(loc_e + tau) ;
  if(tau<rounding_error) {
      qL = 0.0;
  } else {
      qhat = fncAvrQhat(loc_e,tau);
      qL = qhat*0.6*tau;
  }

  res = t14 + 2.0*qL*q15/cg1 ;
    
  return(res);
    
    
}

double Matter::sud_z_QQ_w_M_vac_only(double M, double cg, double cg1, double loc_e , double l_fac, double E2)
{
    
  double q2, q3, q5, q6, q8, q15, qL, tau, res, z_min;
    
    
  z_min = std::sqrt(2)*E_minimum/E2;
    
  //    if (cg<cg1*z_min) cg = cg1*z_min;
    
    
    
  //    if ((cg< cg1/(2.0*E2*E2/cg1+1.0) )) cg = cg1/( 2.0*E2*E2/cg1 + 1.0 );
    
  if (cg1<2.0*(cg+M*M))
    {
        
      //        cout << " returning with cg, cg1 = " << cg << "   " <<  cg1 << "    " << E_minimum << "  " << E2 << endl ;
      return(0.0);
    };
     
  double t1 = M * M;
  double t2 = t1 + cg;
  double t3 = 1.0/cg1;
  double t5 = -t2*t3 + 1.0;
  double t6 = t5*t5;
  double t8 = t2*t2;
  double t10 = pow(cg1, 2.0);
  double t15 = 2.0/3.0*(t6*t5-t8*t2/t10/cg1)*t3;    

    
  q2 = 1.0 / cg1;
  q3 = (cg * q2);
  q5 = 1.0 - q3;
  q6 = std::log(std::abs(q5));
  q8 = std::log(q3);
  /*	q10 = std::log(q3);
	q12 = std::log(q5);*/
  q15 = (-1.0 + (2.0 * q3) + q6 - q8) * q2;
    
    
    
  if (q15<0.0)
    {
      cerr << "ERROR: medium contribution negative in sud_z_QQ: q15 = " << q15 << endl;
      cerr << "cg, cg1 = " << cg << "  " << cg1 << endl;
      cerr << " t15 = " << t15 << endl;
      throw std::runtime_error("ERROR: medium contribution negative in sud_z_QQ");
    }
    
  tau = l_fac ;
    
  if ((length - loc_e) < tau) tau = (length - loc_e);
    
  if (loc_e > length) tau = 0.0 ;

  // SC  
  //qL = qhat*0.6*tau*profile(loc_e + tau) ;
  if(tau<rounding_error) {
      qL = 0.0;
  } else {
      qhat = fncAvrQhat(loc_e,tau);
      qL = qhat*0.6*tau;
  }

  res = t15 + 2.0*qL*q15/cg1 ;
    
  return(res);
    
    
}


double Matter::P_z_qq_int(double cg, double cg1, double loc_e, double cg3, double l_fac, double E2)
{
  double t_q1,t_q3,t_q4,t_q6,t_q8,t_q9,t_q12,q_q1,q_q4,q_q6,q_q9,q_q11,qL,tau,res;
    
  if ((cg< cg1/(2.0*E2*E2/cg1+1.0) )) cg = cg1/( 2.0*E2*E2/cg1 + 1.0 );
    
  t_q1 = std::pow(cg1, 2);
  t_q3 = 1.0 - cg1;
  t_q4 = t_q3 * t_q3;
  t_q6 = std::pow(cg, 2);
  t_q8 = 1.0 - cg;
  t_q9 = t_q8 * t_q8;
  t_q12 = t_q1 * cg1 / 6.0 - t_q4 * t_q3 / 6.0 - t_q6 * cg / 6.0 + t_q9 * t_q8 / 6.0;
    
  tau = l_fac;
    
  if ((length - loc_e) < tau) tau = (length - loc_e);
    
  if (loc_e > length) tau = 0.0 ;
 
  // SC  
  //qL = qhat*0.6*tau*profile(loc_e + tau) ;
  if(tau<rounding_error) {
      qL = 0.0;
  } else {
      qhat = fncAvrQhat(loc_e,tau);
      qL = qhat*0.6*tau;
  }

   
  q_q1 = std::log(cg1);
  q_q4 = std::log(1.0 - cg1);
  q_q6 = std::log(cg);
  q_q9 = std::log(1.0 - cg);
  q_q11 = -cg1 + q_q1 / 2.0 - q_q4 / 2.0 + cg - q_q6 / 2.0 + q_q9 / 2.0;
    
  if (q_q11<0.0)
    {
      cerr << "ERROR: medium contribution negative in P_z_gg_int : q_q11 = " << q_q11 << endl;
      throw std::runtime_error("ERROR: medium contribution negative in P_z_gg_int");
    }
    
  res = t_q12*Tf/Ca + 2.0*qL*q_q11/cg3*(Tf*Cf/Ca/Ca);
    
  return(res);
    
}
//
//

// Sudakov for a quark to radiate a quark + photon
double Matter::sudakov_Pqp(double g0, double g1, double loc_c, double E)
{
    double sud,g;
    int blurb;
    
    sud = 1.0 ;
    
    if (g1<2.0*g0)
    {
        JSWARN << " warning: the lower limit of the sudakov > 1/2 upper limit, returning 1 ";
        JSWARN << " in sudakov_Pquark Photon, g0, g1 = " << g0 << "  " << g1;
        return(sud) ;
    }
    g = 2.0*g0;
    
    double logsud =sud_val_QP(g0,g,g1, loc_c, E );
    
    sud = exp( (-1.0/2.0/pi)*logsud);
    
    return(sud);
    
}

double Matter::sud_val_QP(double h0, double h1, double h2, double loc_d, double E1)
{
    double val, h , intg, hL, hR, diff, intg_L, intg_R, t_form, span;
    int blurb;
    
    double alphaEM = 1.0/137.0;
    
    val = 0.0;
    
    h = (h1+h2)/2.0 ;
    
    span = (h2-h1)/h2;
    
    t_form = 2.0*E1/h;
    
    val = alphaEM*sud_z_QP(h0,h, loc_d, t_form,E1);
    
    intg = val*(h2-h1);
    
    hL = (h1 + h)/2.0 ;
    
    t_form = 2.0*E1/hL;
    
    intg_L = alphaEM*sud_z_QP(h0,hL,loc_d,t_form,E1)*(h - h1) ;
    
    hR = (h + h2)/2.0 ;
    
    t_form = 2.0*E1/hR;
    
    intg_R = alphaEM*sud_z_QP(h0,hR,loc_d,t_form,E1)*(h2 - h) ;
    
    diff = std::abs( (intg_L + intg_R - intg)/intg ) ;
    
    if ( (diff>approx )||(span>error ))
    {
        intg = sud_val_QP(h0,h1,h,loc_d, E1) + sud_val_QP(h0,h,h2,loc_d, E1);
    }
    
    
    return(intg);
    
}

double Matter::P_z_qq_int_w_M_vac_only(double M, double cg, double cg1, double loc_e, double cg3, double l_fac, double E2)
{
  double t_q1,t_q3,t_q4,t_q6,t_q8,t_q9,t_q12,q_q1,q_q4,q_q6,q_q9,q_q11,qL,tau,res;
    
  if ((cg< cg1/(2.0*E2*E2/cg1+1.0) )) cg = cg1/( 2.0*E2*E2/cg1 + 1.0 );
    
  t_q1 = std::pow(cg1, 2);
  t_q3 = 1.0 - cg1;
  t_q4 = t_q3 * t_q3;
  t_q6 = std::pow(cg, 2);
  t_q8 = 1.0 - cg;
  t_q9 = t_q8 * t_q8;
  t_q12 = t_q1 * cg1 / 6.0 - t_q4 * t_q3 / 6.0 - t_q6 * cg / 6.0 + t_q9 * t_q8 / 6.0;
    
  tau = l_fac;
    
  if ((length - loc_e) < tau) tau = (length - loc_e);
    
  if (loc_e > length) tau = 0.0 ;
 
  // SC  
  //qL = qhat*0.6*tau*profile(loc_e + tau) ;
  if(tau<rounding_error) {
      qL = 0.0;
  } else {
      qhat = fncAvrQhat(loc_e,tau);
      qL = qhat*0.6*tau;
  }

   
  q_q1 = std::log(cg1);
  q_q4 = std::log(1.0 - cg1);
  q_q6 = std::log(cg);
  q_q9 = std::log(1.0 - cg);
  q_q11 = -cg1 + q_q1 / 2.0 - q_q4 / 2.0 + cg - q_q6 / 2.0 + q_q9 / 2.0;
    
  if (q_q11<0.0)
    {
      cerr << "ERROR: medium contribution negative in P_z_gg_int : q_q11 = " << q_q11 << endl;
      throw std::runtime_error("ERROR: medium contribution negative in P_z_gg_int");
    }
    
  res = t_q12*Tf/Ca + 2.0*qL*q_q11/cg3*(Tf*Cf/Ca/Ca);
    
  return(res);
    
}



double Matter::sud_z_QP(double cg, double cg1, double loc_e, double l_fac,double E2)
{
    
    double t2,t6,t10,t11,t17, q2, q3, q4, q5,q6,q10,q14, qL, tau, res, z_min;
    int blurb;
    
    z_min = std::sqrt(2)*E_minimum/E2;
    
    if (cg1<2.0*cg)
    {
        return(0.0);
    };
    
    t2 = std::pow(cg1, 2);
    t6 = std::log(cg);
    t10 = std::abs(cg - cg1);
    t11 = std::log(t10);
    t17 = -1.0 / t2 * (3.0 * cg1 - 6.0 * cg + 4.0 * t6 * cg1 - 4.0 * t11 * cg1) / 2.0;
    
    //    return(t17);
    
    q14=0.0;
    
    tau = l_fac ;
    
    if ((length - loc_e) < tau) tau = (length - loc_e);
    
    if (loc_e > length) tau = 0.0 ;
    
    // SC
    //qL = qhat*0.6*tau*profile(loc_e + tau) ;
    if(tau<rounding_error)
    {
        qL = 0.0;
    }
    else
    {
        qhat = fncAvrQhat(loc_e,tau)*Cf/Ca; //for photon production, only the quark scatters
        qL = qhat*0.6*tau;
    }
    
    //JSINFO << BOLDRED << " qhat L = " << qL << " location = " << loc_e << " tau = " << tau << " length = " << length;
    
    res = t17 + 2.0*qL*q14/cg1 ;
    
    //   cout << " t0 , t , res = " << cg << "  "  << cg1 << "   " << res << endl ;
    
    
    if (q14<0.0)
    {
        cerr << "ERROR: medium contribution negative in sud_z_QG : q14 = " << q14 << endl;
        throw std::runtime_error("ERROR: medium contribution negative in sud_z_QG");
    }
    
    return(res);
    
}
double Matter::P_z_qp_int(double cg, double cg1, double loc_e, double cg3, double l_fac, double E2 )
{
    
    double t2, t5, t7, t10, t12, q2, q6, q10, tau, qL, res ;
    
    
    if ((cg< cg1/(2.0*E2*E2/cg1+1.0) )) cg = cg1/( 2.0*E2*E2/cg1 + 1.0 );
    
    t2 = std::pow(cg1, 2);
    t5 = std::log(1.0 - cg1);
    t7 = std::pow(cg, 2);
    t10 = std::log(1.0 - cg);
    t12 = -cg1 - t2 / 2.0 - 2.0 * t5 + cg + t7 / 2.0 + 2.0 * t10;
    
    //    return(t12);
    
    q10 = 0.0;
    tau = l_fac;
    
    if ((length - loc_e) < tau) tau = (length - loc_e);
    
    if (loc_e > length) tau = 0.0 ;
    
    // SC
    //qL = qhat*0.6*tau*profile(loc_e + tau) ;
    if(tau<rounding_error) {
        qL = 0.0;
    } else {
        qhat = fncAvrQhat(loc_e,tau);
        qL = qhat*0.6*tau;
    }
    
    
    res = t12 + 2.0*qL*q10/cg3 ;
    
    return(res);

    
}





//
//
// Sudakov for a quark to radiate a quark + gluon.
double Matter::sudakov_Pqg(double g0, double g1, double loc_c, double E)
{
  double sud,g;
  int blurb;
    
  sud = 1.0 ;
    
  if (g1<2.0*g0)
    {
      JSWARN << " warning: the lower limit of the sudakov > 1/2 upper limit, returning 1 ";
      JSWARN << " in sudakov_Pquark gluon, g0, g1 = " << g0 << "  " << g1;
      return(sud) ;
    }
  g = 2.0*g0;
    
  sud = exp( -1.0*(Cf/2.0/pi)*sud_val_QG(g0,g,g1, loc_c, E ) );
    
  return(sud);
    
}

double Matter::sudakov_Pqg_w_M(double M, double g0, double g1, double loc_c, double E)
{
  double sud,g;
  int blurb;
    
  sud = 1.0 ;
    
    if (g1<g0*( 1.0 + std::sqrt( 1.0 + 2.0*M*M/g0 ) ) )
    {
      JSWARN << " warning: Not enough separation between upper and lower limits of Sudakov to have resolvable radiation ";
      JSWARN << " in sudakov_Pquark gluon, g0*( 1.0 + std::sqrt( 1.0 + 2.0*M*M/g0 ) ), g1 = " << g0*( 1.0 + std::sqrt( 1.0 + 2.0*M*M/g0 ) ) << "  " << g1;
      return(sud) ;
    }
  g = g0*( 1.0 + std::sqrt( 1.0 + 2.0*M*M/g0 ) );
    
  sud = exp( -1.0*(Cf/2.0/pi)*sud_val_QG_w_M(M,g0,g,g1,loc_c, E ) );
    
  return(sud);
    
}


double Matter::sud_val_QG(double h0, double h1, double h2, double loc_d, double E1)
{
  double val, h , intg, hL, hR, diff, intg_L, intg_R, t_form, span;
  int blurb;
    
    
  val = 0.0;
    
  h = (h1+h2)/2.0 ;
    
  span = (h2-h1)/h2;
    
  t_form = 2.0*E1/h;
    
  val = alpha_s(h)*sud_z_QG(h0,h, loc_d, t_form,E1);
    
  intg = val*(h2-h1);
    
  hL = (h1 + h)/2.0 ;
    
  t_form = 2.0*E1/hL;
    
  intg_L = alpha_s(hL)*sud_z_QG(h0,hL,loc_d,t_form,E1)*(h - h1) ;
    
  hR = (h + h2)/2.0 ;
    
  t_form = 2.0*E1/hR;
    
  intg_R = alpha_s(hR)*sud_z_QG(h0,hR,loc_d,t_form,E1)*(h2 - h) ;
    
  diff = std::abs( (intg_L + intg_R - intg)/intg ) ;
    
  if ( (diff>approx )||(span>error ))
    {
      intg = sud_val_QG(h0,h1,h,loc_d, E1) + sud_val_QG(h0,h,h2,loc_d, E1);
    }
    
  return(intg);
    
}

double Matter::sud_val_QG_w_M(double M, double h0, double h1, double h2, double loc_d, double E1)
{
  double val, h , intg, hL, hR, diff, intg_L, intg_R, t_form, span;
  int blurb;
    
    
  val = 0.0;
    
  h = (h1+h2)/2.0 ;
    
  span = (h2-h1)/h2;
    
  t_form = 2.0*E1/h;
    
  val = alpha_s(h)*sud_z_QG_w_M(M,h0,h,loc_d, t_form,E1);
    
  intg = val*(h2-h1);
    
  hL = (h1 + h)/2.0 ;
    
  t_form = 2.0*E1/hL;
    
  intg_L = alpha_s(hL)*sud_z_QG_w_M(M,h0,hL,loc_d,t_form,E1)*(h - h1) ;
    
  hR = (h + h2)/2.0 ;
    
  t_form = 2.0*E1/hR;
    
  intg_R = alpha_s(hR)*sud_z_QG_w_M(M,h0,hR,loc_d,t_form,E1)*(h2 - h) ;
    
  diff = std::abs( (intg_L + intg_R - intg)/intg ) ;
    
  //	cout << " iline, gap, diff = " << i_line << " " << h2 << " " << h1 << "  " << diff << endl ;
  //	cout << " intg, Left , right = " << intg << " " << intg_L << "  " << intg_R << endl;
    
  if ( (diff>approx )||(span>error ))
    {
      intg = sud_val_QG_w_M(M,h0,h1,h,loc_d, E1) + sud_val_QG_w_M(M,h0,h,h2,loc_d, E1);
    }
    
  //    cout << " returning with intg = " << intg << endl;
    
  return(intg);
    
}


double Matter::sud_z_QG(double cg, double cg1, double loc_e, double l_fac,double E2)
{
    
  double t2,t6,t10,t11,t17, q2, q3, q4, q5,q6,q10,q14, qL, tau, res, z_min;
  int blurb;
    
  z_min = std::sqrt(2)*E_minimum/E2;
    
  if (cg1<2.0*cg)
    {
              return(0.0);
    };
    
  t2 = std::pow(cg1, 2);
  t6 = std::log(cg);
  t10 = std::abs(cg - cg1);
  t11 = std::log(t10);
  t17 = -1.0 / t2 * (3.0 * cg1 - 6.0 * cg + 4.0 * t6 * cg1 - 4.0 * t11 * cg1) / 2.0;
    
  //	return(t17);
    
  q2 = 1.0/ cg1;
  q3 = cg * q2;
  q4 = q3 - 1.0;
  q5 = std::abs(q4);
  q6 = std::log(q5);
  q10 = std::log(q3);
  q14 = (q6 + 2.0/cg*cg1 - q10 + 2.0/q4) * q2;
    
    
  tau = l_fac ;
    
  if ((length - loc_e) < tau) tau = (length - loc_e);
    
  if (loc_e > length) tau = 0.0 ;

  // SC  
  //qL = qhat*0.6*tau*profile(loc_e + tau) ;
  if(tau<rounding_error)
  {
      qL = 0.0;
  }
  else
  {
      qhat = fncAvrQhat(loc_e,tau);
      if (qhat*sqrt(2)>0.6)
      {
        // JSINFO << BOLDYELLOW << " length = " << length << " loc = " << loc_e << " tau = " << tau ;
        //JSINFO << BOLDYELLOW << " parton formed at x = " << initRx << " y = " << initRy << " z = " << initRz << " t = " << initR0 ;
        // JSINFO << BOLDYELLOW << " mean qhat for sudakov in GeV^2/fm = " << qhat*5*sqrt(2) ;
      }
      qL = qhat*0.6*tau;
  }

    //JSINFO << BOLDRED << " qhat L = " << qL << " location = " << loc_e << " tau = " << tau << " length = " << length;
    
  res = t17 + 2.0*qL*q14/cg1 ;
    
  //   cout << " t0 , t , res = " << cg << "  "  << cg1 << "   " << res << endl ;
    
    
  if (q14<0.0)
    {
      cerr << "ERROR: medium contribution negative in sud_z_QG : q14 = " << q14 << endl;
      throw std::runtime_error("ERROR: medium contribution negative in sud_z_QG");
    }
    
  return(res);
    
    
    
}

double Matter::sud_z_QG_w_M(double M, double cg, double cg1, double loc_e, double l_fac,double E2)//(t0,t,loc)
{
    
  double  qL, tau, res, z_min;
  int blurb;
    
  z_min = std::sqrt(2)*E_minimum/E2;
    
  if (cg1<2.0*cg + M*M/(1.0 + M*M/cg1 ))
    {
        
      JSINFO << MAGENTA << " returning with cg, cg1 = " << cg << "   " <<  cg1 << "    " << E_minimum << "  " << E2 ;
      return(M*M);
    };
    
  double t1  = 1.0/cg1;
  double t2  = t1*cg;
  double t4  = std::pow(1.0 - t2, 2.0);
  double t7  = std::log(t2);
  double t9  = M*M;
  double t10 = t1*t9;
  double t13 = 1.0/(t10+1.0)*t10;
  double t15 = std::pow(t2 + t13, 2.0);
  double t18 = std::log(1.0 - t2 - t13);
  double t21 = t1*(-t4/2.0 - 1.0 + 2.0*t2 - 2.0*t7 + t15/2.0 + t13 + 2.0*t18);

  double q1  = M*M;
  double q2  = 1.0/cg1;
  double q3  = q2*q1;
  double q5  = 4.0*q3+1.0;
  double q6  = q2*cg;
  double q7  = std::log(q6);
  double q9  = q1*q1;
  double q10 = std::pow(cg1, 2.0);
  double q12 = 1.0/q10*q9;
  double q14 = q12-2.0*q3+1.0/2.0;
  double q15 = 1.0-q6;
  double q16 = std::log(q15);
  double q18 = q3+1.0;
  double q28 = q15*q15;
  double q33 = 1.0/q18*q3;
  double q34 = 1.0-q6-q33;
  double q35 = std::log(q34);
  double q37 = q6+q33;
  double q38 = std::log(q37);
  double q48 = q37*q37;
  double q52 = q7*q5+q16*q14+q15*q18/2.0+2.0/cg*cg1+3.0/2.0*q2/q15*q1-1.0/q28*q12/2.0-q35*q5-q38*q14-q37*q18/2.0-2.0/q34-3.0/2.0*q2/q37*q1+1.0/q48*q12/2.0;
  double q53 = q2*q52;
    
  tau = l_fac ;
    
  if ((length - loc_e) < tau) tau = (length - loc_e);
    
  if (loc_e > length) tau = 0.0 ;

  // SC  
  //qL = qhat*0.6*tau*profile(loc_e + tau) ;
  if(tau<rounding_error) {
      qL = 0.0;
  } else {
      qhat = fncAvrQhat(loc_e,tau);
      if (qhat*sqrt(2)>0.6)
      {
        // JSINFO << BOLDYELLOW << " length = " << length << " loc = " << loc_e << " tau = " << tau ;
        //JSINFO << BOLDYELLOW << " parton formed at x = " << initRx << " y = " << initRy << " z = " << initRz << " t = " << initR0 ;
        // JSINFO << BOLDYELLOW << " mean qhat for sudakov in GeV^2/fm = " << qhat*5*sqrt(2) ;
      }
      qL = qhat*0.6*2.0*tau;
  }

  double e1  = M*M;
  double e2  = 1.0/cg1;
  double e3  = e2*e1;
  double e4  = e2*cg;
  double e5  = std::log(e4);
  double e8  = std::log(1.0-e4);
  double e13 = 1.0/(e3+1.0)*e3;
  double e15 = std::log(1.0-e4-e13);
  double e18 = std::log(e4+e13);
  double e22 = e2*(-(2.0*e5-e8+1.0-e4)*e3+(2.0*e15-e18+e4+e13)*e3);

  double eL;

  if(tau<rounding_error) {
      eL = 0.0;
  } else {
      ehat = 0.0;//fncAvrEhat(loc_e,tau);
      if (ehat*sqrt(2)>0.6)
      {
        // JSINFO << BOLDYELLOW << " length = " << length << " loc = " << loc_e << " tau = " << tau ;
        //JSINFO << BOLDYELLOW << " parton formed at x = " << initRx << " y = " << initRy << " z = " << initRz << " t = " << initR0 ;
        // JSINFO << BOLDYELLOW << " mean qhat for sudakov in GeV^2/fm = " << qhat*5*sqrt(2) ;
      }
      eL = ehat*0.6*4.0;
  }

  double f1  = M*M;
  double f2  = 1.0/cg1;
  double f3  = f2*f1;
  double f4  = f2*cg;
  double f5  = std::log(f4);
  double f8  = f1*f1;
  double f9  = std::pow(cg1, 2.0);
  double f11 = 1.0/f9*f8;
  double f14 = 13.0/4.0*f11-15.0/4.0*f3+1.0/2.0;
  double f15 = 1.0-f4;
  double f16 = std::log(f15);
  double f24 = f15*f15;
  double f32 = 1.0/(f3+1.0)*f3;
  double f33 = 1.0-f4-f32;
  double f34 = std::log(f33);
  double f37 = f4+f32;
  double f38 = std::log(f37);
  double f45 = f37*f37;
  double f52 = f2*((15.0/2.0*f5*f3+f16*f14+1.0/cg*cg1+15.0/4.0*f2/f15*f1-13.0/8.0/f24*f11)*f3-(15.0/2.0*f34*f3+f38*f14+1.0/f33+15.0/4.0*f2/f37*f1-13.0/8.0/f45*f11)*f3);
  double e2L;

  if(tau<rounding_error) {
      e2L = 0.0;
  } else {
      e2hat = qhat/2.0;//fncAvrE2hat(loc_e,tau);
      if (e2hat*sqrt(2)>0.6)
      {
        // JSINFO << BOLDYELLOW << " length = " << length << " loc = " << loc_e << " tau = " << tau ;
        //JSINFO << BOLDYELLOW << " parton formed at x = " << initRx << " y = " << initRy << " z = " << initRz << " t = " << initR0 ;
        // JSINFO << BOLDYELLOW << " mean qhat for sudakov in GeV^2/fm = " << qhat*5*sqrt(2) ;
      }
      e2L = e2hat*0.6*8.0/(tau*cg1);
  }


    //JSINFO << BOLDRED << " qhat L = " << qL << " location = " << loc_e << " tau = " << tau << " length = " << length;
    
    res = t21 + qL*q53/cg1 + eL*e22/cg1+e2L*f52/cg1;
    
  //   cout << " t0 , t , res = " << cg << "  "  << cg1 << "   " << res << endl ;
    
    
  if (res<0.0)
    {
      cerr << "ERROR: medium contribution negative in sud_z_QG : res = " << res << endl;
        
      throw std::runtime_error("ERROR: medium contribution negative in sud_z_QG");
    }
    
  return(res);
    
    
    
}


double Matter::P_z_qg_int(double cg, double cg1, double loc_e, double cg3, double l_fac, double E2 )
{
    
  double t2, t5, t7, t10, t12, q2, q6, q10, tau, qL, res ;
    
    
  if ((cg< cg1/(2.0*E2*E2/cg1+1.0) )) cg = cg1/( 2.0*E2*E2/cg1 + 1.0 );
    
  t2 = std::pow(cg1, 2);
  t5 = std::log(1.0 - cg1);
  t7 = std::pow(cg, 2);
  t10 = std::log(1.0 - cg);
  t12 = -cg1 - t2 / 2.0 - 2.0 * t5 + cg + t7 / 2.0 + 2.0 * t10;
    
  //	return(t12);
    
  q2 = std::log(cg1);
  q6 = std::log(cg);
  q10 = q2 - 2.0 / (cg1 - 1.0) - q6 + 2.0 / (cg - 1.0);
    
  tau = l_fac;
    
  if ((length - loc_e) < tau) tau = (length - loc_e);
    
  if (loc_e > length) tau = 0.0 ;

  // SC
  //qL = qhat*0.6*tau*profile(loc_e + tau) ;
  if(tau<rounding_error) {
      qL = 0.0;
  } else {
      qhat = fncAvrQhat(loc_e,tau);
      qL = qhat*0.6*tau;
  }

   
  res = t12 + 2.0*qL*q10/cg3 ;
    
  return(res);
    
    
    
}

double Matter::P_z_qg_int_w_M(double M, double cg, double cg1, double loc_e, double cg3, double l_fac, double E2 )
{
    
  double t2, t5, t7, t10, t12, tau, qL, res ;
    
    
    if (std::abs(cg-cg1)<rounding_error) return(cg) ;
  //if ((cg< cg1/(2.0*E2*E2/cg1+1.0) )) cg = cg1/( 2.0*E2*E2/cg1 + 1.0 );
    
  t2 = std::pow(cg1, 2);
  t5 = std::log(1.0 - cg1);
  t7 = std::pow(cg, 2);
  t10 = std::log(1.0 - cg);
  t12 = -cg1 - t2 / 2.0 - 2.0 * t5 + cg + t7 / 2.0 + 2.0 * t10;
    
  //	return(t12);

  double q1  = M*M;
  double q2  = 1.0/cg3;
  double q3  = q2*q1;
  double q5  = 4.0*q3+1.0;
  double q6  = 1.0-cg1;
  double q7  = std::log(q6);
  double q9  = q1*q1;
  double q10 = cg3*cg3;
  double q12 = 1.0/q10*q9;
  double q14 = q12-2.0*q3+1.0/2.0;
  double q15 = std::log(cg1);
  double q17 = q3+1.0;
  double q26 = std::pow(cg1, 2.0);
  double q30 = 1.0 - cg;
  double q31 = std::log(q30);
  double q33 = std::log(cg);
  double q43 = std::pow(cg, 2.0);
  double q47 = q7*q5+q15*q14+cg1*q17/2.0+2.0/q6+3.0/2.0*q2/cg1*q1-1.0/q26*q12/2.0-q31*q5-q33*q14-cg*q17/2.0-2.0/q30-3.0/2.0*q2/cg*q1+1.0/q43*q12/2.0;

  tau = l_fac;
    
  if ((length - loc_e) < tau) tau = (length - loc_e);
    
  if (loc_e > length) tau = 0.0 ;

  // SC
  //qL = qhat*0.6*tau*profile(loc_e + tau) ;
  if(tau<rounding_error) {
      qL = 0.0;
  } else {
      qhat = fncAvrQhat(loc_e,tau);
      qL = qhat*0.6*2.0*tau;
  }

  double e1  = M*M;
  double e3  = 1.0/cg3*e1;
  double e4  = 1.0-cg1;
  double e5  = std::log(e4);
  double e10 = 1.0-cg;
  double e11 = std::log(e10);
  double e17 = 2.0*(e5+1.0/e4+cg1/2.0)*e3-2.0*(e11+1.0/e10+cg/2.0)*e3;

  double eL;

  if(tau<rounding_error) {
      eL = 0.0;
  } else {
      ehat = 0.0;//fncAvrEhat(loc_e,tau);
      eL = ehat*0.6*4.0;
  }

  double f1  = M*M;
  double f2  = 1.0/cg3;
  double f3  = f2*f1;
  double f4  = 13.0*f3;
  double f6  = f1*f1;
  double f7  = f6*(f4+15.0);
  double f8  = 1.0-cg1;
  double f9  = std::log(f8);
  double f10 = cg3*cg3;
  double f11 = 1.0/f10;
  double f15 = std::log(cg1);
  double f23 = f1*(39.0/4.0*f11*f6+15.0/2.0*f3+1.0);
  double f28 = f6*(f4+15.0/2.0);
  double f29 = f8*f8;
  double f37 = 1.0/f10/cg3*f6*f1;
  double f42 = 1.0-cg;
  double f43 = std::log(f42);
  double f47 = std::log(cg);
  double f54 = f42*f42;
  double f63 = f11*f9*f7/4.0+f2*f1*f15/2.0+1.0/f8*f2*f23-1.0/f29*f11*f28/2.0+13.0/6.0/f29/f8*f37-f11*f43*f7/4.0-f2*f1*f47/2.0-1.0/f42*f2*f23+1.0/f54*f11*f28/2.0-13.0/6.0/f54/f42*f37;
    
  double e2L;

  if(tau<rounding_error) {
      e2L = 0.0;
  } else {
      e2hat = qhat/2.0;//fncAvrE2hat(loc_e,tau);
      e2L = e2hat*0.6*8.0/(tau*cg3);
  }


   
  res = t12 + qL*q47/cg3 + eL*e17/cg3 + e2L*f63/cg3;
    
  return(res);
    
    
    
}


double Matter::alpha_s(double q2)
{
  double a,L2,q24, c_nf;
    
  L2 = std::pow(Lambda_QCD,2);
    
  q24 = q2/4.0;
    
  c_nf = nf ;
    
  if (q24>4.0) {
    c_nf = 4;
  }
    
  if (q24>64.0) {
    c_nf = 5;
  }
    
    
    
  if (q24>L2)
    {
      a = 12.0*pi/(11.0*Nc-2.0*c_nf)/std::log(q24/L2) ;
    }
  else
    {        
      JSWARN << " alpha too large ";
      a=0.6;        
    }

  return(a);
}

double Matter::profile(double zeta)
{
  double prof;
    
  /*  Modify the next set of lines to get a particular profile in brick test or further modify the profile for hydro to introduce coherence effects */
    
  prof = 1.0;
    
  return(prof);
    
}


////////////////////////////////////////////////////////////////////////////////////////

double Matter::fillQhatTab() {


    double xLoc,yLoc,zLoc,tLoc;
    double vxLoc,vyLoc,vzLoc,gammaLoc,betaLoc;
    double edLoc,sdLoc;
    double tempLoc;
    double flowFactor,qhatLoc;
    int hydro_ctl;
    double lastLength=initR0;

    double tStep = 0.1;

    std::unique_ptr<FluidCellInfo> check_fluid_info_ptr;

    for(int i = 0; i < dimQhatTab; i++)
    {
        tLoc = tStep*i;

        //if(tLoc<initR0-tStep) { // potential problem of making t^2<z^2
        if(tLoc<initR0 || tLoc<tStart)
        {
            qhatTab1D[i] = 0.0; 
            continue;	    
        }

	xLoc = initRx+(tLoc-initR0)*initVx; 
        yLoc = initRy+(tLoc-initR0)*initVy;
        zLoc = initRz+(tLoc-initR0)*initVz;

//        if(bulkFlag == 1) { // read OSU hydro
//            readhydroinfoshanshan_(&tLoc,&xLoc,&yLoc,&zLoc,&edLoc,&sdLoc,&tempLoc,&vxLoc,&vyLoc,&vzLoc,&hydro_ctl);
//        } else if(bulkFlag == 2) { // read CCNU hydro
//            hydroinfoccnu_(&tLoc, &xLoc, &yLoc, &zLoc, &tempLoc, &vxLoc, &vyLoc, &vzLoc, &hydro_ctl);
//        } else if(bulkFlag == 0) { // static medium
//            vxLoc = 0.0;
//            vyLoc = 0.0;
//            vzLoc = 0.0;
//            hydro_ctl = 0;
//            tempLoc = T;
//        }

        GetHydroCellSignal(tLoc, xLoc, yLoc, zLoc, check_fluid_info_ptr);
        VERBOSE(8)<< MAGENTA<<"Temperature from medium = "<<check_fluid_info_ptr->temperature;
	
        tempLoc = check_fluid_info_ptr->temperature;
        sdLoc = check_fluid_info_ptr->entropy_density;
        vxLoc = check_fluid_info_ptr->vx;
        vyLoc = check_fluid_info_ptr->vy;
        vzLoc = check_fluid_info_ptr->vz;

        hydro_ctl=0;

        if(hydro_ctl == 0 && tempLoc >= hydro_Tc)
        {
            lastLength = tLoc;
            betaLoc = sqrt(vxLoc*vxLoc+vyLoc*vyLoc+vzLoc*vzLoc);
            gammaLoc = 1.0/sqrt(1.0-betaLoc*betaLoc);
            flowFactor = gammaLoc*(1.0-(initVx*vxLoc+initVy*vyLoc+initVz*vzLoc));

            if(qhat0 < 0.0)
            {
                // calculate qhat with alphas
                double muD2 = 6.0*pi*alphas*tempLoc*tempLoc;
                // if(initEner > pi*tempLoc) qhatLoc = Ca*alphas*muD2*tempLoc*log(6.0*initEner*tempLoc/muD2);
                // else qhatLoc = Ca*alphas*muD2*tempLoc*log(6.0*pi*tempLoc*tempLoc/muD2);
                // fitted formula from https://arxiv.org/pdf/1503.03313.pdf
                if(initEner > pi*tempLoc) qhatLoc = Ca*50.4864/pi*pow(alphas,2)*pow(tempLoc,3)*log(5.7*initEner*tempLoc/4.0/muD2);
                else qhatLoc = Ca*50.4864/pi*pow(alphas,2)*pow(tempLoc,3)*log(5.7*pi*tempLoc*tempLoc/4.0/muD2);
                qhatLoc = qhatLoc*flowFactor;
                if(qhatLoc<0.0) qhatLoc=0.0;
            }
            else
            { // use input qhat
                if(brick_med)
                {
                    qhatLoc = qhat0*0.1973*flowFactor;
                }
                else
                {
                    qhatLoc = qhat0/96.0*sdLoc*0.1973*flowFactor;  // qhat0 at s = 96fm^-3
                }
            }
        //    cout << "check qhat --  ener, T, qhat: " << initEner << "  " << tempLoc << "  " << qhatLoc << endl;
        }
        else
        { // outside the QGP medium
            qhatLoc = 0.0;
        }

        qhatTab1D[i] = qhatLoc/sqrt(2.0); // store qhat value in light cone coordinate         

    }

    for(int i = 0; i < dimQhatTab; i++)
    { // dim of loc

        double totValue = 0.0;

        for(int j = 0; i+j < dimQhatTab; j++)
        { // dim of tau_f

            totValue = totValue+qhatTab1D[i+j];             
            qhatTab2D[i][j] = totValue/(j+1);

        }
    }

    //return(lastLength*sqrt(2.0)*5.0); // light cone + GeV unit
    return( (2.0*lastLength+initRdotV-initR0)/sqrt(2.0)*5.0 ); // light cone + GeV unit
  
}

//////////////////////////////////////////////////////////////////////////////////////

double Matter::fncQhat(double zeta)
{
    if(in_vac) return(0.0);

    double tStep = 0.1;
    //int indexZeta = (int)(zeta/sqrt(2.0)/5.0/tStep+0.5); // zeta was in 1/GeV and light cone coordinate
    int indexZeta = (int)((sqrt(2.0)*zeta/5.0-initRdotV+initR0)/2.0/tStep+0.5); // zeta was in 1/GeV and light cone coordinate

    //if(indexZeta >= dimQhatTab) indexZeta = dimQhatTab-1;
    if(indexZeta >= dimQhatTab) return(0);
    
    double avrQhat = qhatTab1D[indexZeta]; 
    return(avrQhat);
}


//////////////////////////////////////////////////////////////////////////////////////

double Matter::fncAvrQhat(double zeta, double tau) {

    if(in_vac) return(0.0);

    double tStep = 0.1;
    //int indexZeta = (int)(zeta/sqrt(2.0)/5.0/tStep+0.5); // zeta was in 1/GeV and light cone coordinate
    int indexZeta = (int)((sqrt(2.0)*zeta/5.0-initRdotV+initR0)/2.0/tStep+0.5); // zeta was in 1/GeV and light cone coordinate
    int indexTau = (int)(tau/sqrt(2.0)/5.0/tStep+0.5); // tau was in 1/GeV and light cone coordinate

   // if(indexZeta >= dimQhatTab) indexZeta = dimQhatTab-1;
    if(indexZeta >= dimQhatTab) return(0);
    if(indexTau >= dimQhatTab) indexTau = dimQhatTab-1;
    
    double avrQhat = qhatTab2D[indexZeta][indexTau]; 
    return(avrQhat);

}

//////////////////////////////////////////////////////////////////////////////////////

void Matter::flavor(int &CT,int &KATT0,int &KATT2,int &KATT3, unsigned int &max_color, unsigned int &color0, unsigned int &anti_color0, unsigned int &color2, unsigned int &anti_color2, unsigned int &color3, unsigned int &anti_color3){  

  int vb[7]={0};
  int b=0;
  int KATT00=KATT0;
  unsigned int backup_color0=color0;
  unsigned int backup_anti_color0=anti_color0;

  vb[1]=1; 
  vb[2]=2; 
  vb[3]=3; 
  vb[4]=-1; 
  vb[5]=-2; 
  vb[6]=-3;
	  
  if(KATT00==21) {   //.....for gluon
    double R1=16.0; // gg->gg DOF_g
    double R2=0.0;  // gg->qqbar don't consider this channel in MATTER
    double R3=6.0*6*4/9; // gq->gq or gqbar->gqbar flavor*DOF_q*factor 
    double R0=R1+R3;
	   	
    double a=ran0(&NUM1);
	   
    if(a<=R1/R0) { // gg->gg
      CT=1;
      KATT3=21;
      KATT2=21;
      //KATT0=KATT0;
      max_color++;
      color0=backup_color0;
      anti_color0=max_color;
      max_color++;
      color2=anti_color0;
      anti_color2=max_color;
      color3=backup_anti_color0;
      anti_color3=max_color;
    } else { // gq->gq
      CT=3;
      b=floor(ran0(&NUM1)*6+1);
      if(b==7) b=6;
      KATT3=vb[b];
      KATT2=KATT3;
      //KATT0=KATT0;
      if(KATT3>0) { // gq->gq
          max_color++;
          color0=backup_color0;
	  anti_color0=max_color;
	  color2=max_color;
	  anti_color2=0;
          color3=backup_anti_color0;
	  anti_color3=0;
      } else { // gqbar->gqbar	 
          max_color++;
	  color0=max_color;
          anti_color0=backup_anti_color0;
	  color2=0;
	  anti_color2=max_color;
	  color3=0;
          anti_color3=backup_color0;
      }     
    }
  } else { //.....for quark and antiquark (light)
    double R3=16.0; // qg->qg DOF_g
    double R4=4.0*6*4/9; // qq'->qq' scatter with other species
    double R5=1.0*6*4/9; // qq->qq scatter with itself
    double R6=0.0; // qqbar->q'qbar' to other final state species, don't consider in MATTER
    double R7=1.0*6*4/9; // qqbar->qqbar scatter with its anti-particle
    double R8=0.0; // qqbar->gg don't consider in MATTER
    double R00=R3+R4+R5+R7;

    double a=ran0(&NUM1);
    if(a<=R3/R00) { // qg->qg
	CT=13;
	KATT3=21;
	KATT2=21; 
	//KATT0=KATT0;
        if(KATT00>0) { // qg->qg
            max_color++;
            color0=max_color;
            anti_color0=0;
	    max_color++;
            color2=max_color;
            anti_color2=color0;
            color3=max_color;
            anti_color3=backup_color0;
        } else { // qbarg->qbarg	 
            max_color++;
            color0=0;
            anti_color0=max_color;
	    max_color++;
            color2=anti_color0;
            anti_color2=max_color;
            color3=backup_anti_color0;
            anti_color3=max_color;
        }     
    } else if(a<=(R3+R4)/R00) { // qq'->qq'
	CT=4;
	do {
      	    b=floor(ran0(&NUM1)*6+1);
	    if(b==7) b=6;
	    KATT3=vb[b];
	} while (KATT3==KATT0 || KATT3==-KATT0); 
	KATT2=KATT3;
	//KATT0=KATT0;
	if(KATT00>0 && KATT2>0) { // qq->qq
            max_color++;
	    color0=max_color;
	    anti_color0=0;
	    color2=backup_color0;
	    anti_color2=0;
	    color3=max_color;
	    anti_color3=0;
	} else if(KATT00>0 && KATT2<0) { //qqbar->qqbar
	    max_color++;
            color0=max_color;
            anti_color0=0;
            color2=0;
            anti_color2=max_color;
            color3=0;
            anti_color3=backup_color0;
        } else if(KATT00<0 && KATT2>0) { //qbarq->qbarq
            max_color++;
            color0=0;
	    anti_color0=max_color;
	    color2=max_color;
	    anti_color2=0;
	    color3=backup_anti_color0;
	    anti_color3=0;
        } else { //qbarqbar->qbarqbar
            max_color++;
	    color0=0;
	    anti_color0=max_color;
	    color2=0;
	    anti_color2=backup_anti_color0;
	    color3=0;
	    anti_color3=max_color;
        }
    } else if(a<=(R3+R4+R5)/R00) { // scatter with itself
	CT=5;
	KATT3=KATT0;
	KATT2=KATT0;
	//KATT0=KATT0;
	if(KATT00>0) { // qq->qq
            max_color++;
	    color0=max_color;
	    anti_color0=0;
	    color2=backup_color0;
	    anti_color2=0;
	    color3=max_color;
	    anti_color3=0;
         } else { //qbarqbar->qbarqbar
            max_color++;
	    color0=0;
	    anti_color0=max_color;
	    color2=0;
	    anti_color2=backup_anti_color0;
	    color3=0;
	    anti_color3=max_color;
         }
    } else { // scatter with its anti-particle
	CT=7;
	KATT3=-KATT0; 
	KATT2=KATT3;
	//KATT0=KATT0;
	if(KATT00>0) { //qqbar->qqbar
	    max_color++;
            color0=max_color;
            anti_color0=0;
            color2=0;
            anti_color2=max_color;
            color3=0;
            anti_color3=backup_color0;
        } else { //qbarq->qbarq
            max_color++;
            color0=0;
	    anti_color0=max_color;
	    color2=max_color;
	    anti_color2=0;
	    color3=backup_anti_color0;
	    anti_color3=0;
        }
    }	   
  }	
}

void Matter::colljet22(int CT,double temp,double qhat0ud,double v0[4],double p0[4],double p2[4],double p3[4],double p4[4],double &qt){
  //
  //    p0 initial jet momentum, output to final momentum
  //    p2 final thermal momentum,p3 initial termal energy
  //
  //    amss=sqrt(abs(p0(4)**2-p0(1)**2-p0(2)**2-p0(3)**2))
  //
  //************************************************************
  p4[1]=p0[1];
  p4[2]=p0[2];
  p4[3]=p0[3];
  p4[0]=p0[0];	  	  
  //************************************************************	

  //    transform to local comoving frame of the fluid
  //  cout << endl;
  //  cout << "flow  "<< v0[1] << " " << v0[2] << " " << v0[3] << " "<<" Elab " << p0[0] << endl;
  
  trans(v0,p0);
  //  cout << p0[0] << " " << sqrt(qhat0ud) << endl;
  
  //  cout << sqrt(pow(p0[1],2)+pow(p0[2],2)+pow(p0[3],2)) << " " << p0[1] << " " << p0[2] << " " << p0[3] << endl;

  //************************************************************
  trans(v0,p4);
  //************************************************************


  //    sample the medium parton thermal momentum in the comoving frame


  double xw;
  double razim;
  double rcos;
  double rsin;
	  
  double ss;  
  double tmin;
  double tmid;
  double tmax;
	  
  double rant;
  double tt;
	  
  double uu;	  
  double ff=0.0;
  double rank;
	  
  double mmax;
  double msq=0.0;
	  
  double f1;
  double f2;
  
  double p0ex[4]={0.0};
  double vc[4]={0.0};

  //  Initial 4-momentum of jet
  //
  //************************************************************
  p4[1]=p0[1];
  p4[2]=p0[2];
  p4[3]=p0[3];
  p4[0]=p0[0];	  	  
  //************************************************************	  

  int ic=0;

  do
    {
      do{	  
	xw=15.0*ran0(&NUM1);
	razim=2.0*pi*ran0(&NUM1);
	rcos=1.0-2.0*ran0(&NUM1);
	rsin=sqrt(1.0-rcos*rcos);
	//
	p2[0]=xw*temp;
	p2[3]=p2[0]*rcos;
	p2[1]=p2[0]*rsin*cos(razim);
	p2[2]=p2[0]*rsin*sin(razim);
	  
	f1=pow(xw,3 )/(exp(xw)-1)/1.4215;
	f2=pow(xw,3)/(exp(xw)+1)/1.2845;
	//
	//    cms energy
	//
	ss=2.0*(p0[0]*p2[0]-p0[1]*p2[1]-p0[2]*p2[2]-p0[3]*p2[3]);
	  
	//	if(ss.lt.2.d0*qhat0ud) goto 14

	tmin=qhat0ud;
	tmid=ss/2.0;
	tmax=ss-qhat0ud;
	  
	//    use (s^2+u^2)/(t+qhat0ud)^2 as scattering cross section in the
	//
	rant=ran0(&NUM1);
	tt=rant*ss;	  

	//		ic+=1;
	//		cout << p0[0] << "  " << p2[0] <<  endl;
	//		cout << tt << "  " << ss <<  "" << qhat0ud <<endl;
	//		cout << ic << endl;

      }while((tt<qhat0ud) || (tt>(ss-qhat0ud))); 

      uu=ss-tt;	  
	  
      if(CT==1)
	{	
	  ff=f1;	  
	  mmax=4.0/pow(ss,2)*(3.0-tmin*(ss-tmin)/pow(ss,2)+(ss-tmin)*ss/pow(tmin,2)+tmin*ss/pow((ss-tmin),2));
	  msq=pow((1.0/p0[0]/p2[0]/2.0),2)*(3.0-tt*uu/pow(ss,2)+uu*ss/pow(tt,2)+tt*ss/pow(uu,2))/mmax;
	}

	
      if(CT==2)
	{
	  ff=f1;	  
	  mmax=4.0/pow(ss,2)*(4.0/9.0*(pow(tmin,2)+pow((ss-tmin),2))/tmin/(ss-tmin)-(pow(tmin,2)+pow((ss-tmin),2))/pow(ss,2));
	  msq=pow((1.0/p0[0]/p2[0]/2.0),2)*(4.0/9.0*(pow(tt,2)+pow(uu,2))/tt/uu-(pow(tt,2)+pow(uu,2))/pow(ss,2))/(mmax+4.0);
	}

	
      if(CT==3)
	{
	  ff=f2;	  
	  if(((pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/ss/(ss-tmin))>((pow(ss,2)+pow((ss-tmax),2)/pow(tmax,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmax),2))/ss/(ss-tmax))))
	    {
	      mmax=4.0/pow(ss,2)*((pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/ss/(ss-tmin));
	    }
	  else
	    {	  
	      mmax=4.0/pow(ss,2)*((pow(ss,2)+pow((ss-tmax),2))/pow(tmax,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmax),2))/ss/(ss-tmax));
	    }
	  //
	  msq=pow((1.0/p0[0]/p2[0]/2.0),2)*((pow(ss,2)+pow(uu,2))/pow(tt,2)+4.0/9.0*(pow(ss,2)+pow(uu,2))/ss/uu)/mmax;
	}

	
      if(CT==13)
	{
	  ff=f1;
	  
	  if(((pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/ss/(ss-tmin))>((pow(ss,2)+pow((ss-tmax),2)/pow(tmax,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmax),2))/ss/(ss-tmax))))
	    {
	      mmax=4.0/pow(ss,2)*((pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/ss/(ss-tmin));
	    }
	  else 
	    {
	      mmax=4.0/pow(ss,2)*((pow(ss,2)+pow((ss-tmax),2))/pow(tmax,2)+4.0/9.0*(pow(ss,2)+pow((ss-tmax),2))/ss/(ss-tmax));
	    }
	  //
	  msq=pow((1.0/p0[0]/p2[0]/2.0),2)*((pow(ss,2)+pow(uu,2))/pow(tt,2)+4.0/9.0*(pow(ss,2)+pow(uu,2))/ss/uu)/mmax;
	}

	
      if(CT==4)
	{
	  ff=f2;	  
	  mmax=4.0/pow(ss,2)*(4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2));
	  msq=pow((1.0/p0[0]/p2[0]/2.0),2)*(4.0/9.0*(pow(ss,2)+pow(uu,2))/pow(tt,2))/mmax;
	}

	
      if(CT==5)
	{
	  ff=f2;	  
	  mmax=4.0/pow(ss,2)*(4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+(pow(ss,2)+pow(tmin,2))/pow((ss-tmin),2)-2.0/3.0*pow(ss,2)/tmin/(ss-tmin));
	  msq=pow((1.0/p0[0]/p2[0]/2.0),2)*(4.0/9.0*((pow(ss,2)+pow(uu,2))/pow(tt,2)+(pow(ss,2)+pow(tt,2))/pow(uu,2)-2.0/3.0*pow(ss,2)/tt/uu))/mmax;	  
	}

	
      if(CT==6)
	{
	  ff=f2;	  
	  mmax=4.0/pow(ss,2)*(4.0/9.0*(pow(tmin,2)+pow((ss-tmin),2))/pow(ss,2));
	  msq=pow((1.0/p0[0]/p2[0]/2.0),2)*(4.0/9.0*(pow(tt,2)+pow(uu,2))/pow(ss,2))/(mmax+0.5);
	}

	
      if(CT==7)
	{
	  ff=f2;	  
	  mmax=4.0/pow(ss,2)*(4.0/9.0*(pow(ss,2)+pow((ss-tmin),2))/pow(tmin,2)+(pow(tmin,2)+pow((ss-tmin),2))/pow(ss,2)+2.0/3.0*pow((ss-tmin),2)/ss/tmin);
	  msq=(pow((1.0/p0[0]/p2[0]/2.0),2)*(4.0/9.0*(((pow(ss,2)+pow(uu,2))/pow(tt,2))+(pow(tt,2)+pow(uu,2))/pow(ss,2)+2.0/3.0*pow(uu,2)/ss/tt)))/mmax;	  
	}

	
      if(CT==8)
	{
	  ff=f2;	 
	  mmax=4.0/pow(ss,2)*(4.0/9.0*(pow(tmin,2)+pow((ss-tmin),2))/tmin/(ss-tmin)-(pow(tmin,2)+pow((ss-tmin),2))/pow(ss,2));
	  msq=pow((1.0/p0[0]/p2[0]/2.0),2)*(4.0/9.0*(pow(tt,2)+pow(uu,2))/tt/uu-(pow(tt,2)+pow(uu,2))/pow(ss,2))/(mmax+4.0);	  
	}
	
      rank=ran0(&NUM1);
    }while(rank>(msq*ff));
	
  //
  p3[1]=p2[1]; 
  p3[2]=p2[2];
  p3[3]=p2[3];
  p3[0]=p2[0];

  //    velocity of the center-of-mass
  //
  vc[1]=(p0[1]+p2[1])/(p0[0]+p2[0]);
  vc[2]=(p0[2]+p2[2])/(p0[0]+p2[0]);
  vc[3]=(p0[3]+p2[3])/(p0[0]+p2[0]);
  //
  //    transform into the cms frame
  //
  trans(vc,p0);
  trans(vc,p2);
  //
  //    cm momentum
  //
  double pcm=p2[0];
  //
  //    sample transverse momentum transfer with respect to jet momentum
  //    in cm frame
  //
  double ranp=2.0*pi*ran0(&NUM1);
  //
  //    transverse momentum transfer
  //
  qt=sqrt(pow(pcm,2)-pow((tt/2.0/pcm-pcm),2));
  double qx=qt*cos(ranp);
  double qy=qt*sin(ranp);

  //
  //    longitudinal momentum transfer
  //
  double qpar=tt/2.0/pcm;
  //
  //    qt is perpendicular to pcm, need to rotate back to the cm frame
  //
  double upt=sqrt(p2[1]*p2[1]+p2[2]*p2[2])/p2[0];
  double upx=p2[1]/p2[0];
  double upy=p2[2]/p2[0];
  double upz=p2[3]/p2[0];
  //
  //    momentum after collision in cm frame
  //
  p2[1]=p2[1]-qpar*upx;
  p2[2]=p2[2]-qpar*upy;
  if(upt!=0.0) 
    {
      p2[1]=p2[1]+(upz*upx*qy+upy*qx)/upt;
      p2[2]=p2[2]+(upz*upy*qy-upx*qx)/upt;
    }
  p2[3]=p2[3]-qpar*upz-upt*qy;

  p0[1]=-p2[1];
  p0[2]=-p2[2];
  p0[3]=-p2[3];
  //
  //    transform from cm back to the comoving frame
  //
  transback(vc,p2);
  transback(vc,p0);

  //************************************************************
  //
  //     calculate qt in the rest frame of medium
  //
  //  if(p0[4]>p2[4])
  //	{
  rotate(p4[1],p4[2],p4[3],p0,1);
  qt=sqrt(pow(p0[1],2)+pow(p0[2],2));
  rotate(p4[1],p4[2],p4[3],p0,-1);
  //	}
  //  else
  //	{
  //	  rotate(p4[1],p4[2],p4[3],p2,1);
  //	  qt=sqrt(pow(p2[1],2)+pow(p2[2],2));
  //	  rotate(p4[1],p4[2],p4[3],p2,-1);
  //	}
  //************************************************************	  
	  
	  
	  
  //
  //    transform from comoving frame to the lab frame
  //
  transback(v0,p2);
  transback(v0,p0);
  transback(v0,p3);
	  
  //************************************************************
  transback(v0,p4);
  //************************************************************	  
	  
}

void Matter::trans(double v[4],double p[4]){	
  double vv=sqrt(v[1]*v[1]+v[2]*v[2]+v[3]*v[3]);	
  double ga=1.0/sqrt(1.0-vv*vv);
  double ppar=p[1]*v[1]+p[2]*v[2]+p[3]*v[3];
  double gavv=(ppar*ga/(1.0+ga)-p[0])*ga;
  p[0]=ga*(p[0]-ppar);
  p[1]=p[1]+v[1]*gavv;
  p[2]=p[2]+v[2]*gavv;
  p[3]=p[3]+v[3]*gavv;
}

void Matter::transback(double v[4],double p[4]){
  double vv=sqrt(v[1]*v[1]+v[2]*v[2]+v[3]*v[3]);
  double ga=1.0/sqrt(1.0-vv*vv); 
  double ppar=p[1]*v[1]+p[2]*v[2]+p[3]*v[3];
  double gavv=(-ppar*ga/(1.0+ga)-p[0])*ga;
  p[0]=ga*(p[0]+ppar);
  p[1]=p[1]-v[1]*gavv;
  p[2]=p[2]-v[2]*gavv;
  p[3]=p[3]-v[3]*gavv;
}

void Matter::rotate(double px,double py,double pz,double pr[4],int icc){
  //     input:  (px,py,pz), (wx,wy,wz), argument (i)
  //     output: new (wx,wy,wz)
  //     if i=1, turn (wx,wy,wz) in the direction (px,py,pz)=>(0,0,E)
  //     if i=-1, turn (wx,wy,wz) in the direction (0,0,E)=>(px,py,pz)

	   
  double wx,wy,wz,E,pt,w,cosa,sina,cosb,sinb;
  double wx1,wy1,wz1;	   	   

  wx=pr[1];
  wy=pr[2];
  wz=pr[3];

  E=sqrt(px*px+py*py+pz*pz);
  pt=sqrt(px*px+py*py);

  w=sqrt(wx*wx+wy*wy+wz*wz);

  if(pt==0)
    {
      cosa=1;
      sina=0;
    } 
  else
    {
      cosa=px/pt;
      sina=py/pt;
    }

  cosb=pz/E;
  sinb=pt/E;

  if(icc==1)
    { 
      wx1=wx*cosb*cosa+wy*cosb*sina-wz*sinb;
      wy1=-wx*sina+wy*cosa;
      wz1=wx*sinb*cosa+wy*sinb*sina+wz*cosb;
    }  

  else
    {
      wx1=wx*cosa*cosb-wy*sina+wz*cosa*sinb;
      wy1=wx*sina*cosb+wy*cosa+wz*sina*sinb;
      wz1=-wx*sinb+wz*cosb;
    }

  wx=wx1;
  wy=wy1;
  wz=wz1;

  pr[1]=wx;
  pr[2]=wy;
  pr[3]=wz;      

  //  pr[0]=sqrt(pr[1]*pr[1]+pr[2]*pr[2]+pr[3]*pr[3]);

}


float Matter::ran0(long *idum)

{
  const int    IM1=2147483563;
  const int    IM2=2147483399;
  const double AM=(1.0/IM1);
  const int    IMM1=(IM1-1);
  const int    IA1=40014;
  const int    IA2=40692;
  const int    IQ1=53668;
  const int    IQ2=52774;
  const int    IR1=12211;
  const int    IR2=3791;
  const int    NTAB=32;
  const int    NDIV=(1+IMM1/NTAB);
  const double EPS=1.2e-7;
  const double RNMX=(1.0-EPS);

  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0) { 
    if (-(*idum) < 1) *idum=1; 
    else *idum = -(*idum);
    for (j=NTAB+7;j>=0;j--) { 
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1; 
  *idum=IA1*(*idum-k*IQ1)-k*IR1; 
  if (*idum < 0) *idum += IM1; 
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; 
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV; 
  iy=iv[j]-idum2; 
  iv[j] = *idum; 
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; 
  else return temp;
}
	

//////////////////////////////////////////////////////////////////////////////////////

double Matter::solve_alphas(double var_qhat, double var_ener, double var_temp) {

    double preFactor=42.0*Ca*zeta3/pi;

    // reference: qhatLoc = Ca*50.4864/pi*pow(alphas,2)*pow(tempLoc,3)*log(5.7*2.0*pi*tempLoc*tempLoc/4.0/muD2);
    double max_qhat=preFactor*pow(0.5,2)*pow(var_temp,3)*log(5.7*max(var_ener,2.0*pi*var_temp)/24/pi/0.5/var_temp);

    if(max_qhat<var_qhat) {
        JSINFO << "qhat exceeds HTL calculation, use alpha_s = 0.5";
	return(0.5);
    }

    double solution=sqrt(var_qhat/preFactor/pow(var_temp,3)/log(5.7*max(var_ener,2.0*pi*var_temp)/24/pi/0.2/var_temp));
    double fnc_value,fnc_derivative;
    fnc_value=fnc0_alphas(solution,var_qhat,var_ener,var_temp);
    fnc_derivative=fnc0_derivative_alphas(solution,var_qhat,var_ener,var_temp);

    //cout << "initial guess: " << solution << "  " << fnc_value << endl;

    while(fabs(fnc_value/var_qhat)>0.001) {
    
	solution=solution-fnc_value/fnc_derivative;
        fnc_value=fnc0_alphas(solution,var_qhat,var_ener,var_temp);
        fnc_derivative=fnc0_derivative_alphas(solution,var_qhat,var_ener,var_temp);

    }

    if(solution<0.0 || solution>0.5) {
        JSINFO << "unreasonable alpha_s: " << solution << " use alpha_s = 0.5";
	solution=0.5;
    }

    return(solution);
    
}

double Matter::fnc0_alphas(double var_alphas, double var_qhat, double var_ener, double var_temp) {
   
    double preFactor=42.0*Ca*zeta3/pi;
    return(preFactor*var_alphas*var_alphas*pow(var_temp,3)*log(5.7*max(var_ener,2.0*pi*var_temp)/24/pi/var_alphas/var_temp)-var_qhat);

}

double Matter::fnc0_derivative_alphas(double var_alphas, double var_qhat, double var_ener, double var_temp) {
   
    double preFactor=42.0*Ca*zeta3/pi;
    return(preFactor*pow(var_temp,3)*(2.0*var_alphas*log(5.7*max(var_ener,2.0*pi*var_temp)/24/pi/var_alphas/var_temp)-var_alphas));

}
