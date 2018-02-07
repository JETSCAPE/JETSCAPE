// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// Overwritten: Abhijit Majumder (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

// quick and dirty test class implementation for Eloss modules ...
// can be used as a user template ...

#include "ElossModulesTestMatter.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include <string>

//#include "../3rdparty/MATTER/sudakovs.cpp"
//#include "../3rdparty/MATTER/generators.cpp"
//#include "../3rdparty/MATTER/matter_profile.cpp"
//#include "../3rdparty/MATTER/corrector.cpp"


#include "tinyxml2.h"
#include<iostream>

#include "fluid_dynamics.h"

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
  INFO<<"Intialize Matter ...";

  // Redundant (get this from Base) quick fix here for now
  tinyxml2::XMLElement *eloss= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Eloss" );
  if ( !eloss ) {
    WARN << "Couldn't find tag Eloss";
    throw std::runtime_error ("Couldn't find tag Eloss");    
  }
  tinyxml2::XMLElement *matter=eloss->FirstChildElement("Matter");
  if ( !matter ) {
    WARN << "Couldn't find tag Eloss -> Matter";
    throw std::runtime_error ("Couldn't find tag Eloss -> Matter");    
  }
 
  string s = matter->FirstChildElement( "name" )->GetText();
  JSDEBUG << s << " to be initializied ...";
  
  in_vac = false;
  brick_med = true;
  
  qhat = 0.0;
  Q00 = 1.0; // virtuality separation scale
  qhat0 = 2.0; // GeV^2/fm for gluon at s = 96 fm^-3
  alphas = 0.3; // only useful when qhat0 is a negative number
  hydro_Tc = 0.16;
  brick_length = 4.0;
  vir_factor = 1.0;
  
  double m_qhat=-99.99;
  if ( !matter->FirstChildElement("qhat0")) {
    WARN << "Couldn't find sub-tag Eloss -> Matter -> qhat0";
    throw std::runtime_error ("Couldn't find sub-tag Eloss -> Matter -> qhat0");
  }
  matter->FirstChildElement("qhat0")->QueryDoubleText(&m_qhat);
  SetQhat(m_qhat);
  //qhat = GetQhat()/fmToGeVinv ;
  qhat0 = GetQhat()/fmToGeVinv ;
  JSDEBUG  << s << " with qhat0 = "<<GetQhat();
  
  int flagInt=-100;
  double inputDouble=-99.99;
  
  if ( !matter->FirstChildElement("in_vac") ) {
    WARN << "Couldn't find sub-tag Eloss -> Matter -> in_vac";
    throw std::runtime_error ("Couldn't find sub-tag Eloss -> Matter -> in_vac");
  }
  matter->FirstChildElement("in_vac")->QueryIntText(&flagInt);
  in_vac = flagInt;
  
  if ( !matter->FirstChildElement("brick_med") ) {
    WARN << "Couldn't find sub-tag Eloss -> Matter -> brick_med";
    throw std::runtime_error ("Couldn't find sub-tag Eloss -> Matter -> brick_med");
  }
  matter->FirstChildElement("brick_med")->QueryIntText(&flagInt);
  brick_med = flagInt;
  
  if ( !matter->FirstChildElement("Q0") ) {
    WARN << "Couldn't find sub-tag Eloss -> Matter -> Q0";
    throw std::runtime_error ("Couldn't find sub-tag Eloss -> Matter -> Q0");
  }
  matter->FirstChildElement("Q0")->QueryDoubleText(&inputDouble);
  Q00 = inputDouble;
  
  if ( !matter->FirstChildElement("alphas") ) {
    WARN << "Couldn't find sub-tag Eloss -> Matter -> alphas";
    throw std::runtime_error ("Couldn't find sub-tag Eloss -> Matter -> alphas");
  }
  matter->FirstChildElement("alphas")->QueryDoubleText(&inputDouble);
  alphas = inputDouble;
  
  if ( !matter->FirstChildElement("hydro_Tc") ) {
    WARN << "Couldn't find sub-tag Eloss -> Matter -> hydro_Tc";
    throw std::runtime_error ("Couldn't find sub-tag Eloss -> Matter -> hydro_Tc");
  }
  matter->FirstChildElement("hydro_Tc")->QueryDoubleText(&inputDouble);
  hydro_Tc = inputDouble;
  
  if ( !matter->FirstChildElement("brick_length") ) {
    WARN << "Couldn't find sub-tag Eloss -> Matter -> brick_length";
    throw std::runtime_error ("Couldn't find sub-tag Eloss -> Matter -> brick_length");
  }
  matter->FirstChildElement("brick_length")->QueryDoubleText(&inputDouble);
  brick_length = inputDouble;
  
  if ( !matter->FirstChildElement("vir_factor") ) {
    WARN << "Couldn't find sub-tag Eloss -> Matter -> vir_factor";
    throw std::runtime_error ("Couldn't find sub-tag Eloss -> Matter -> vir_factor");
  }
  matter->FirstChildElement("vir_factor")->QueryDoubleText(&inputDouble);
  vir_factor = inputDouble;
  
  if(vir_factor<0.0) {
    cout << "Error: vir_factor < 0, reset to 1.0" << endl;
    vir_factor=1.0;
  }
  
  VERBOSE(7)<< MAGENTA << "MATTER input parameter";
  VERBOSE(7)<< MAGENTA << "Q00:" << Q00;
  INFO << "in_vac: " << in_vac << "  brick_med: " << brick_med;
  INFO << "Q00: " << Q00 << " vir_factor: " << vir_factor << "  qhat0: " << qhat0*fmToGeVinv << " alphas: " << alphas << " hydro_Tc: " << hydro_Tc << " brick_length: " << brick_length;

  // Initialize random number distribution
  ZeroOneDistribution = uniform_real_distribution<double> { 0.0, 1.0 };
}


void Matter::WriteTask(weak_ptr<JetScapeWriter> w)
{
  VERBOSE(8);
  w.lock()->WriteComment("ElossModule Parton List: "+GetId());
  w.lock()->WriteComment("Energy loss to be implemented accordingly ...");
}

void Matter::DoEnergyLoss(double deltaT, double time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut)
{

  double z=0.5;
  double blurb,zeta,tQ2 ;
  int iSplit,pid_a,pid_b;
    unsigned int max_color, min_color, min_anti_color;
  double velocity[4],xStart[4];

  VERBOSESHOWER(8)<< MAGENTA << "SentInPartons Signal received : "<<deltaT<<" "<<Q2<<" "<<&pIn;
      
  // FluidCellInfo* check_fluid_info_ptr = new FluidCellInfo;
  // VERBOSE(8)<< MAGENTA<<"Temperature from Brick (Signal) = "<<check_fluid_info_ptr->temperature;
  // delete check_fluid_info_ptr;
  std::unique_ptr<FluidCellInfo> check_fluid_info_ptr;
  GetHydroCellSignal(1, 1.0, 1.0, 0.0, check_fluid_info_ptr);      

  double rNum;
        
  double delT = deltaT;
  double Time = time*fmToGeVinv;
  double deltaTime = delT*fmToGeVinv;
    
  JSDEBUG << " the time in fm is " << time << " The time in GeV-1 is " << Time ;
  JSDEBUG << "pid = " << pIn[0].pid() << " E = " << pIn[0].e() << " px = " << pIn[0].p(1) << " py = " << pIn[0].p(2) << "  pz = " << pIn[0].p(3) << " virtuality = " << pIn[0].t() << " form_time in fm = " << pIn[0].form_time()/fmToGeVinv ;
  JSDEBUG << " color = " << pIn[0].color() << " anti-color = " << pIn[0].anti_color();
    
  //JSDEBUG << " For MATTER, the qhat in GeV^-3 = " << qhat ;
    
  for (int i=0;i<pIn.size();i++)
  {
      velocity[0] = 1.0;
      for(int j=1;j<=3;j++)
      {
          velocity[j] = pIn[i].p(j)/pIn[i].e();
      }
          
      double velocityMod = std::sqrt(std::pow(velocity[1],2) + std::pow(velocity[2],2) + std::pow(velocity[3],2));
          
      for(int j=0;j<=3;j++)
      {
          xStart[j] = pIn[i].x_in().comp(j) ;
          //+ velocity[j]*delT;
      }

      // SC: read in hydro
      initR0 = xStart[0];
      initRx = xStart[1];
      initRy = xStart[2];
      initRz = xStart[3];
      initVx = velocity[1]/velocityMod;
      initVy = velocity[2]/velocityMod;
      initVz = velocity[3]/velocityMod;
      initRdotV = initRx*initVx + initRy*initVy + initRz*initVz;
      initEner = pIn[i].e();
      if(!in_vac) length = fillQhatTab();
      if(brick_med) length = brick_length*fmToGeVinv; /// length in GeV-1 will have to changed for hydro
      //if(brick_med) length = 5.0*fmToGeVinv; /// length in GeV-1 will have to changed for hydro


      //zeta = ( xStart[0] + std::sqrt( xStart[1]*xStart[1] + xStart[2]*xStart[2] + xStart[3]*xStart[3] )  )/std::sqrt(2);
      // SC
      zeta = ( xStart[0] + initRdotV )/std::sqrt(2)*fmToGeVinv;


      int pid = pIn[i].pid();
          
      if (pIn[i].form_time()<0.0) /// A parton without a virtuality or formation time, must set...
      {
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
	  
	  // SC:  
          double pT2 = pIn[i].p(1)*pIn[i].p(1)+pIn[i].p(2)*pIn[i].p(2);
          //tQ2 = generate_vac_t(pIn[i].pid(), pIn[i].nu(), QS/2.0, pIn[i].e()*pIn[i].e() ,zeta , iSplit);
          tQ2 = generate_vac_t(pIn[i].pid(), pIn[i].nu(), QS/2.0, pT2*vir_factor ,zeta , iSplit);
            	    
          // KK:
          pIn[i].set_jet_v(velocity);
          pIn[i].set_t(tQ2); // Also resets momentum!
          pIn[i].set_mean_form_time();
          double ft = generate_L(pIn[i].mean_form_time()) ;
          pIn[i].set_form_time(ft);
          
          unsigned int color=0, anti_color=0;
          std::uniform_int_distribution<short> uni(102,103);
          
          if ( pIn[i].pid()>0 )
          {
             // color = uni(*get_mt19937_generator());
              color = 101;
          }
          pIn[i].set_color(color);
          if ( (pIn[i].pid()<0)||(pIn[i].pid()==21) )
          {
              anti_color = uni(*get_mt19937_generator());
          }
          pIn[i].set_anti_color(anti_color);
          
          max_color = color;
          
          if (anti_color > color) max_color = anti_color ;
          
          min_color = color;
          
          min_anti_color = anti_color;
          
          pIn[i].set_max_color(max_color);
          pIn[i].set_min_color(min_color);
          pIn[i].set_min_anti_color(min_anti_color);
          
            
          //JSDEBUG:
          JSDEBUG ;
          JSDEBUG << " ***************************************************************************** " ;
          JSDEBUG<< " ID = " << pIn[i].pid() << " Color = " << pIn[i].color() << " Anti-Color = " << pIn[i].anti_color() ;
          JSDEBUG << " E = " << pIn[i].e() << " px = " << pIn[i].px() << " py = " << pIn[i].py() << " pz = " << pIn[i].pz() ;
          JSDEBUG << " *  New generated virtuality = " << tQ2 << " Mean formation time = " << pIn[i].mean_form_time()/fmToGeVinv;
          JSDEBUG << " *  set new formation time to " << pIn[i].form_time()/fmToGeVinv ;
          JSDEBUG << " * Maximum allowed virtuality = " << pIn[i].e()*pIn[i].e() << "   Minimum Virtuality = " << QS;
          JSDEBUG << " * Qhat = " << qhat << "  Length = "  << length ;
          JSDEBUG << " * Jet velocity = " << pIn[i].jet_v().comp(0) << " " << pIn[i].jet_v().comp(1) << "  " << pIn[i].jet_v().comp(2) << "  " << pIn[i].jet_v().comp(3);
          JSDEBUG << " * reset location of parton formation = "<< pIn[i].x_in().t() << "  " << pIn[i].x_in().x() << "  " << pIn[i].x_in().y() << "  " << pIn[i].x_in().z();
          JSDEBUG << " ***************************************************************************** " ;
          JSDEBUG ;
          // end DEBUG:
 
          
      }

      // SC: Q0 can be changed based on different setups
      if(in_vac) { // for vaccuum
	  qhat = 0.0;
	  if(Q00 < 0.0) Q0 = 1.0; // set Q0 = 1 if Q00 < 0
	  else Q0 = Q00;
      } else { // for medium    
	  double tempEner = initEner;
          qhat = fncQhat(zeta);

	  if(Q00 < 0.0) { // use dynamical Q0 if Q00 < 0
	      if(pid==gid) Q0 = sqrt(sqrt(2.0*tempEner*qhat*sqrt(2.0)));
	      else Q0 = sqrt(sqrt(2.0*tempEner*qhat*sqrt(2.0)/Ca*Cf));
	      if(Q0 < 1.0) Q0 = 1.0;
	      if(zeta > length) Q0 = 1.0;
	  } else {
              Q0 = Q00;
          }
      }

      if(Q0<1.0) Q0=1.0;

      //if (pIn[i].t() > QS + rounding_error)
      if (pIn[i].t() > Q0*Q0 + rounding_error)
      { //
          double decayTime = pIn[i].mean_form_time()  ;
	    
          //JSDEBUG << "  deltaT = " << deltaT;
          //JSDEBUG << " parton origin time = " << pIn[i].x_in().t()/fmToGeVinv << " parton formation time = " << pIn[i].form_time()/fmToGeVinv;
          // JSDEBUG << " parton id " << pIn[i].pid() << " parton virtuality = " << pIn[i].t();
          //JSDEBUG << " parton momentum " << pIn[i].e() << "  " << pIn[i].px() << "  " << pIn[i].py() << "  " << pIn[i].pz();
	    
          double splitTime = pIn[i].form_time() + pIn[i].x_in().t() ;
          // JSDEBUG << " splitTime = " << splitTime/fmToGeVinv;
	  
          if (splitTime<Time)
          {
              // do split
              double t_used = pIn[i].t();
              //if (t_used<QS)  t_used = QS; // SC: not necessary
              double tau_form = 2.0*pIn[i].nu()/t_used;
              double z_low = QS/t_used/2.0 ;
              double z_hi = 1.0 - z_low;
              
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
                  double r = ZeroOneDistribution(*get_mt19937_generator());
                  if (r>ratio)
                  { // qqbar
		  
                      double r2 = ZeroOneDistribution(*get_mt19937_generator());
		  
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
                  pid_a = pid ;
                  pid_b = gid ;
                  iSplit = 0;
              }
              int ifcounter = 0;
              double l_perp2 = -1.0; // SC: initialization
              // daughter virtualities
              double tQd1 = QS;
              double tQd2 = QS;
              
              
              //set color of daughters here
              unsigned int d1_col, d1_acol, d2_col, d2_acol, color, anti_color;
              //std::uniform_int_distribution<short> uni(101,103);
              //color = pIn[i].color();
              max_color = pIn[i].max_color();
              //if (pIn[i].anti_color()>maxcolor) color = pIn[i].anti_color();
              JSDEBUG << " old max color = " << max_color;
              max_color++;
              color = max_color;
              anti_color = max_color;
              pIn[i].set_max_color(max_color);
              JSDEBUG << " new color = " << color;
              
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
              else
              {
                 throw std::runtime_error("error in iSplit");
              }

              
              JSDEBUG << " d1_col = " << d1_col << " d1_acol = " << d1_acol << " d2_col = " << d2_col << " d2_acol = " << d2_acol;

              while ((l_perp2<=0.0)&&(ifcounter<100))
              {
                  z = generate_vac_z(pid,QS/2.0,pIn[i].t(),zeta,pIn[i].nu(),iSplit) ;
                  //JSDEBUG << " generated z = " << z;
        
                  int iSplit_a = 0;
                  if (pid_a==gid) iSplit_a = 1;
	
                  if (z*z*pIn[i].t()>QS)
                  {
                      tQd1 = generate_vac_t(pid_a, z*pIn[i].nu(), QS/2.0, z*z*pIn[i].t() ,zeta+std::sqrt(2)*pIn[i].form_time() , iSplit_a);
                  } else { // SC
                      tQd1 = z*z*pIn[i].t();
                  }
        
                  int iSplit_b = 0;
                  if (pid_b==gid) iSplit_b = 1;
        
                  if ((1.0-z)*(1.0-z)*pIn[i].t()>QS)
                  {
                      tQd2 = generate_vac_t(pid_b, (1.0-z)*pIn[i].nu(), QS/2.0, (1.0-z)*(1.0-z)*pIn[i].t() ,zeta+std::sqrt(2)*pIn[i].form_time(),iSplit_b);
                  } else { // SC
                      tQd2 = (1.0-z)*(1.0-z)*pIn[i].t();
                  }
		
                  l_perp2 =  pIn[i].t()*z*(1.0 - z) - tQd2*z - tQd1*(1.0-z) ; ///< the transverse momentum squared
                  ifcounter++;
              }
              
              if (l_perp2<=0.0) l_perp2 = 0.0; ///< test if negative
              double l_perp = std::sqrt(l_perp2); ///< the momentum transverse to the parent parton direction
              
              // axis of split
              double angle = generate_angle();

                  
              // double parent_perp = std::sqrt( pow(pIn[i].p(1),2) + pow(pIn[i].p(2),2) + pow(pIn[i].p(3),2) - pow(pIn[i].pl(),2) );
              // KK: changed to x,y,z
              double parent_perp = std::sqrt( pow(pIn[i].px(),2) + pow(pIn[i].py(),2) + pow(pIn[i].pz(),2) - pow(pIn[i].pl(),2) );
              double mod_jet_v = std::sqrt( pow(pIn[i].jet_v().x(),2) +  pow(pIn[i].jet_v().y(),2) + pow(pIn[i].jet_v().z(),2) ) ;
                  
              double c_t = pIn[i].jet_v().z()/mod_jet_v;
              double s_t = std::sqrt( 1.0 - c_t*c_t) ;
                  
              double s_p = pIn[i].jet_v().y()/std::sqrt( pow( pIn[i].jet_v().x() , 2 ) + pow( pIn[i].jet_v().y(), 2 ) ) ;
              double c_p = pIn[i].jet_v().x()/std::sqrt( pow( pIn[i].jet_v().x() , 2 ) + pow( pIn[i].jet_v().y(), 2 ) ) ;
                  
              // First daughter
              double k_perp1[4];
              k_perp1[0] = 0.0;
              k_perp1[1] = z*(pIn[i].px() - pIn[i].pl()*s_t*c_p) + l_perp*std::cos(angle)*c_t*c_p - l_perp*std::sin(angle)*s_p ;
              k_perp1[2] = z*(pIn[i].py() - pIn[i].pl()*s_t*s_p) + l_perp*std::cos(angle)*c_t*s_p + l_perp*std::sin(angle)*c_p ;
              k_perp1[3] = z*(pIn[i].pz() - pIn[i].pl()*c_t) - l_perp*std::cos(angle)*s_t ;
              double k_perp1_2 = pow(k_perp1[1],2)+pow(k_perp1[2],2)+pow(k_perp1[3],2);
                  
              double energy = ( z*pIn[i].nu() + (tQd1 + k_perp1_2)/(2.0*z*pIn[i].nu() ) )/std::sqrt(2.0) ;
              double plong =  ( z*pIn[i].nu() - (tQd1 + k_perp1_2)/(2.0*z*pIn[i].nu() ) )/std::sqrt(2.0) ;
                  
              double newp[4];
              newp[0] = energy;
              newp[1] = plong*s_t*c_p + k_perp1[1];
              newp[2] = plong*s_t*s_p + k_perp1[2];
              newp[3] = plong*c_t + k_perp1[3];
                  
              double newx[4];
              newx[0] = Time + deltaTime ;
              for (int j=1;j<=3;j++)
              {
                  newx[j] = pIn[i].x_in().comp(j) + (Time + deltaTime - pIn[i].x_in().comp(0) )*velocity[j]/velocityMod;
              }
                  
              pOut.push_back(Parton(0,pid_a,0,newp,newx ));
              int iout = pOut.size()-1 ;
                  
              pOut[iout].set_jet_v(velocity);
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
              k_perp2[1] = (1.0-z)*(pIn[i].px() - pIn[i].pl()*s_t*c_p) - l_perp*std::cos(angle)*c_t*c_p + l_perp*std::sin(angle)*s_p ;
              k_perp2[2] = (1.0-z)*(pIn[i].py() - pIn[i].pl()*s_t*s_p) - l_perp*std::cos(angle)*c_t*s_p - l_perp*std::sin(angle)*c_p ;
              k_perp2[3] = (1.0-z)*(pIn[i].pz() - pIn[i].pl()*c_t) + l_perp*std::cos(angle)*s_t ;
              double k_perp2_2 = pow(k_perp2[1],2)+pow(k_perp2[2],2)+pow(k_perp2[3],2);

              energy = ( (1.0-z)*pIn[i].nu() + (tQd2 + k_perp2_2)/( 2.0*(1.0-z)*pIn[i].nu() ) )/std::sqrt(2.0) ;
              plong =  ( (1.0-z)*pIn[i].nu() - (tQd2 + k_perp2_2)/( 2.0*(1.0-z)*pIn[i].nu() ) )/std::sqrt(2.0) ;

              parent_perp = std::sqrt( pow(pIn[i].p(1),2) + pow(pIn[i].p(2),2) + pow(pIn[i].p(3),2) - pow(pIn[i].pl(),2) );
              mod_jet_v = std::sqrt( pow(pIn[i].jet_v().x(),2) +  pow(pIn[i].jet_v().y(),2) + pow(pIn[i].jet_v().z(),2) ) ;
                  
              newp[0] = energy;
              newp[1] = plong*s_t*c_p + k_perp2[1] ;
              newp[2] = plong*s_t*s_p + k_perp2[2] ;
              newp[3] = plong*c_t + k_perp2[3] ;
              
              newx[0] = Time + deltaTime;
              for (int j=1;j<=3;j++)
              {
                  newx[j] = pIn[i].x_in().comp(j) + (Time + deltaTime - pIn[i].x_in().comp(0) )*velocity[j]/velocityMod;
              }
	      
              pOut.push_back(Parton(0,pid_b,0,newp,newx ));
              iout = pOut.size()-1 ;
              pOut[iout].set_jet_v(velocity);
              pOut[iout].set_mean_form_time();
              ft = generate_L (pOut[iout].mean_form_time());
              pOut[iout].set_form_time(ft);
              pOut[iout].set_color(d2_col);
              pOut[iout].set_anti_color(d2_acol);
              pOut[iout].set_max_color(max_color);
              pOut[iout].set_min_color(pIn[i].min_color());
              pOut[iout].set_min_anti_color(pIn[i].min_anti_color());

              

                  
          }
          else
          { // not time to split yet
              // pOut.push_back(pIn[i]);
          }
      }
      else
      { // virtuality too low
          // pOut.push_back(pIn[i]);
      }
          	  
  } // particle loop
      
}

double Matter::generate_angle()
{
  double ang, r, blurb;
    
  // r = double(random())/ (maxN );
  r = ZeroOneDistribution(*get_mt19937_generator());
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
  r = ZeroOneDistribution(*get_mt19937_generator());
  //        r = mtrand1();
    
    
  if ((r>=1.0)||(r<=0.0))
    {
      throw std::runtime_error("error in random number in t *get_mt19937_generator()");
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
      numer = sudakov_Pqg(t0,t,loc_a,nu);
    }
    
  t_mid = t_low;
    
  if (numer>r)
    {
      // cout << " numer > r, i.e. ; " << numer << " > " << r << endl ;        
      return(t_mid) ;
    }
    
  t_mid = (t_low+t_hi)/2.0 ;
    
    
  scale = t0;
    
  //   cout << " s_approx, s_error = " << s_approx << "  " << s_error << endl;
    
  while ((abs(diff)>s_approx)&&(abs(t_hi-t_low)/t_hi>s_error))
    {
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
	  if (is!=0)
            {
	      throw std::runtime_error(" error in isp in quark split numerator  ");
            }            
	  denom = sudakov_Pqg(t0, t_mid, loc_a, nu);
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
	  t_mid = (t_low + t_hi)/2.0;
        }
      else
        {
	  t_hi = t_mid ;
	  t_mid = (t_low + t_hi)/2.0;
        }
        
    }
    
  return(t_mid);
}

/*
 
 
 
  New function
 
 
 
*/



double  Matter::generate_vac_z(int p_id, double t0, double t, double loc_b, double nu, int is )
{
  double r,z, ratio,diff,e,numer1, numer2, numer, denom, z_low, z_hi, z_mid, test;

  r = ZeroOneDistribution(*get_mt19937_generator());
    
  if ((r>1)||(r<0))
    {
      throw std::runtime_error(" error in random number in z *get_mt19937_generator()");
    }
    
  ratio = 1.0 ;
    
  diff = (ratio - r)/r;
    
  e = t0/t;
    
  if (e>0.5)
    {
      throw std::runtime_error(" error in 	epsilon");
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
  else
    {
      denom = P_z_qg_int(z_low,z_hi, loc_b, t, 2.0*nu/t , nu);
    }
    
    
  z_mid = (z_low + z_hi)/2.0 ;
  
  int itcounter=0;
  // cout << " generate_vac_z called with p_id = " << p_id << " t0 = " << t0 << " t = " << t << " loc_b=" << loc_b<< " nu = " <<  nu << " is = " << is << endl;
  while (abs(diff)>approx) { // Getting stuck in here for some reason
    if ( itcounter++ > 1000 ) throw std::runtime_error("Stuck in endless loop") ;
    // cout << " in here" << " abs(diff) = " << abs(diff) << "  approx = " << approx << endl;
    if (p_id==gid) {
      if (is==1) {
	numer = P_z_gg_int(e, z_mid, loc_b, t, 2.0*nu/t , nu );
      } else {
	numer = P_z_qq_int(e, z_mid, loc_b, t, 2.0*nu/t , nu);
      }
    } else {
      numer = P_z_qg_int(e, z_mid, loc_b, t, 2.0*nu/t , nu );
    }
        
    ratio = numer/denom ;      
    diff = (ratio - r)/r ;     
    // cout << "num, den, r, diff = " << numer << " "<< denom << " " << r << " " << endl;	
    // cout << " diff, z_mid = " << diff << " " << z_mid << endl ;		
    //		cin >> test ;
      
        
    if (diff>0.0)
    {
        z_hi = z_mid;
        z_mid = (z_low + z_hi)/2.0;
    }
    else
    {
        z_low = z_mid;
        z_mid = (z_low + z_hi)/2.0 ;
    }
        
  }
    
  return(z_mid);
}	




double Matter::generate_L(double form_time)
{
  double r, x_low, x_high , x , diff , span, val, arg, norm ;
    
  // r = double(random())/ (maxN );
  r = ZeroOneDistribution(*get_mt19937_generator());
  //    r = mtrand1();
    
  if ((r>1)||(r<0))
    {
      throw std::runtime_error(" error in random number in z *get_mt19937_generator()");
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
      WARN << " warning: the lower limit of the sudakov > 1/2 upper limit, returning 1 ";
      WARN << " in sudakov_Pquark quark, q0, q1 = " << q0 << "  " << q1;
      return(sud) ;
    }
  q = 2.0*q0;
    
  //	g = g0 ;
    
  sud = exp( -1.0*(Tf/2.0/pi)*sud_val_QQ(q0,q,q1,loc_c, E) );
    
  return(sud);
    
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


double Matter::sudakov_Pqg(double g0, double g1, double loc_c, double E)
{
  double sud,g;
  int blurb;
    
  sud = 1.0 ;
    
  if (g1<2.0*g0)
    {
      WARN << " warning: the lower limit of the sudakov > 1/2 upper limit, returning 1 ";
      WARN << " in sudakov_Pquark gluon, g0, g1 = " << g0 << "  " << g1;
      return(sud) ;
    }
  g = 2.0*g0;
    
  sud = exp( -1.0*(Cf/2.0/pi)*sud_val_QG(g0,g,g1, loc_c, E ) );
    
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
    
  //	cout << " iline, gap, diff = " << i_line << " " << h2 << " " << h1 << "  " << diff << endl ;
  //	cout << " intg, Left , right = " << intg << " " << intg_L << "  " << intg_R << endl;
    
  if ( (diff>approx )||(span>error ))
    {
      intg = sud_val_QG(h0,h1,h,loc_d, E1) + sud_val_QG(h0,h,h2,loc_d, E1);
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
        
      //        cout << " returning with cg, cg1 = " << cg << "   " <<  cg1 << "    " << E_minimum << "  " << E2 << endl ;
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
  if(tau<rounding_error) {
      qL = 0.0;
  } else {
      qhat = fncAvrQhat(loc_e,tau);
      qL = qhat*0.6*tau;
  }


    
  res = t17 + 2.0*qL*q14/cg1 ;
    
  //   cout << " t0 , t , res = " << cg << "  "  << cg1 << "   " << res << endl ;
    
    
  if (q14<0.0)
    {
      cerr << "ERROR: medium contribution negative in sud_z_QG : q14 = " << q14 << endl;
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
      WARN << " alpha too large ";
      a=0.6;        
    }

  return(a);
}

double Matter::profile(double zeta)
{
  double prof;
    
  /*  Modify the next set of lines to get a particular profile */
    
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
    double lastLength=0.0;

    double tStep = 0.1;

    std::unique_ptr<FluidCellInfo> check_fluid_info_ptr;

    for(int i = 0; i < dimQhatTab; i++) {
        tLoc = tStep*i;

	if(tLoc<initR0-tStep) {
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
	VERBOSE(7)<< MAGENTA<<"Temperature from Brick (Signal) = "<<check_fluid_info_ptr->temperature;
	
	tempLoc = check_fluid_info_ptr->temperature;
	sdLoc = check_fluid_info_ptr->entropy_density;
	vxLoc = check_fluid_info_ptr->vx;
	vyLoc = check_fluid_info_ptr->vy;
	vzLoc = check_fluid_info_ptr->vz;

	hydro_ctl=0;

        if(hydro_ctl == 0 && tempLoc >= hydro_Tc) { 
            lastLength = tLoc;
            betaLoc = sqrt(vxLoc*vxLoc+vyLoc*vyLoc+vzLoc*vzLoc);
            gammaLoc = 1.0/sqrt(1.0-betaLoc*betaLoc);
            flowFactor = gammaLoc*(1.0-(initVx*vxLoc+initVy*vyLoc+initVz*vzLoc));

            if(qhat0 < 0.0) { // calculate qhat with alphas
                double muD2 = 6.0*pi*alphas*tempLoc*tempLoc;
                // if(initEner > pi*tempLoc) qhatLoc = Ca*alphas*muD2*tempLoc*log(6.0*initEner*tempLoc/muD2);
                // else qhatLoc = Ca*alphas*muD2*tempLoc*log(6.0*pi*tempLoc*tempLoc/muD2);
                // fitted formula from https://arxiv.org/pdf/1503.03313.pdf
                if(initEner > pi*tempLoc) qhatLoc = Ca*50.4864/pi*pow(alphas,2)*pow(tempLoc,3)*log(5.7*initEner*tempLoc/4.0/muD2);
                else qhatLoc = Ca*50.4864/pi*pow(alphas,2)*pow(tempLoc,3)*log(5.7*pi*tempLoc*tempLoc/4.0/muD2);
                qhatLoc = qhatLoc*flowFactor;
                if(qhatLoc<0.0) qhatLoc=0.0;
            } else { // use input qhat
                if(brick_med) qhatLoc = qhat0*0.1973*flowFactor;  
		else qhatLoc = qhat0/96.0*sdLoc*0.1973*flowFactor;  // qhat_over_T3 here is actually qhat0 at s = 96fm^-3
            }
        //    cout << "check qhat --  ener, T, qhat: " << initEner << "  " << tempLoc << "  " << qhatLoc << endl;

        } else { // outside the QGP medium
            qhatLoc = 0.0;
        }

        qhatTab1D[i] = qhatLoc/sqrt(2.0); // store qhat value in light cone coordinate         

    }

    for(int i = 0; i < dimQhatTab; i++) { // dim of loc

        double totValue = 0.0;

        for(int j = 0; i+j < dimQhatTab; j++) { // dim of tau_f

            totValue = totValue+qhatTab1D[i+j];             
            qhatTab2D[i][j] = totValue/(j+1);

        }
    }

    //return(lastLength*sqrt(2.0)*5.0); // light cone + GeV unit
    return( (lastLength+initRdotV)/sqrt(2.0)*5.0 ); // light cone + GeV unit
  
}

//////////////////////////////////////////////////////////////////////////////////////

double Matter::fncQhat(double zeta) {

    if(in_vac) return(0.0);

    double tStep = 0.1;
    //int indexZeta = (int)(zeta/sqrt(2.0)/5.0/tStep+0.5); // zeta was in 1/GeV and light cone coordinate
    int indexZeta = (int)((sqrt(2.0)*zeta-initRdotV)/5.0/tStep+0.5); // zeta was in 1/GeV and light cone coordinate

    if(indexZeta >= dimQhatTab) indexZeta = dimQhatTab-1;      
    
    double avrQhat = qhatTab1D[indexZeta]; 
    return(avrQhat);

}


//////////////////////////////////////////////////////////////////////////////////////

double Matter::fncAvrQhat(double zeta, double tau) {

    if(in_vac) return(0.0);

    double tStep = 0.1;
    //int indexZeta = (int)(zeta/sqrt(2.0)/5.0/tStep+0.5); // zeta was in 1/GeV and light cone coordinate
    int indexZeta = (int)((sqrt(2.0)*zeta-initRdotV)/5.0/tStep+0.5); // zeta was in 1/GeV and light cone coordinate
    int indexTau = (int)(tau/sqrt(2.0)/5.0/tStep+0.5); // tau was in 1/GeV and light cone coordinate

    if(indexZeta >= dimQhatTab) indexZeta = dimQhatTab-1;      
    if(indexTau >= dimQhatTab) indexTau = dimQhatTab-1;      
    
    double avrQhat = qhatTab2D[indexZeta][indexTau]; 
    return(avrQhat);

}

//////////////////////////////////////////////////////////////////////////////////////


