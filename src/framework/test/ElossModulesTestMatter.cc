// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
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
  tinyxml2::XMLElement *matter=eloss->FirstChildElement("Matter");
 
  if (matter) {   
    string s = matter->FirstChildElement( "name" )->GetText();
    DEBUG << s << " to be initializied ...";

    double m_qhat=-99.99;
    matter->FirstChildElement("qhat")->QueryDoubleText(&m_qhat);
    SetQhat(m_qhat);
    qhat = GetQhat()/fmToGeVinv ;
    DEBUG  << s << " with qhat = "<<GetQhat();
      
  }
  else {
    WARN << " : Matter not properly initialized in XML file ...";
    throw std::runtime_error("Matter not properly initialized in XML file ...");
  }

  // Initialize random number distribution
  ZeroOneDistribution = uniform_real_distribution<double> { 0.0, 1.0 };
}


void Matter::WriteTask(weak_ptr<JetScapeWriter> w)
{
  VERBOSE(8);
  w.lock()->WriteComment("ElossModule Parton List: "+GetId());
  w.lock()->WriteComment("Energy loss to be implemented accordingly ...");
}

// stupid toy branching ....
// think about memory ... use pointers ...
//void Matter::DoEnergyLoss(double deltaT, double Q2, const vector<Parton>& pIn, vector<Parton>& pOut)
void Matter::DoEnergyLoss(double deltaT, double Q2, vector<Parton>& pIn, vector<Parton>& pOut)
{

  //DEBUG:
  //cout<<" -----> "<<*GetShowerInitiatingParton()<<endl;
  double z=0.5;
  double blurb,zeta,tQ2 ;
  int iSplit,pid_a,pid_b;
  double velocity[4],xStart[4];
  //    cout << " pIn size = " << pIn.size() << endl;

  VERBOSESHOWER(8)<< MAGENTA << "SentInPartons Signal received : "<<deltaT<<" "<<Q2<<" "<<&pIn;
      
  FluidCellInfo* check_fluid_info_ptr = new FluidCellInfo;
  GetHydroCellSignal(1, 1.0, 1.0, 0.0, check_fluid_info_ptr);      
  VERBOSE(8)<< MAGENTA<<"Temperature from Brick (Signal) = "<<check_fluid_info_ptr->temperature;            
  delete check_fluid_info_ptr;

  double rNum;
        
  double delT = 0.1;
  double Time = deltaT*fmToGeVinv;
  double deltaTime = delT*fmToGeVinv;
    
  //   cout << " QS = " << QS << " Q2 = " << Q2 << endl;

  for (int i=0;i<pIn.size();i++) {
                   		  
    velocity[0] = 1.0;
    for(int j=1;j<=3;j++){
      velocity[j] = pIn[i].p(j)/pIn[i].e();
    }
          
    double velocityMod = std::sqrt(std::pow(velocity[1],2) + std::pow(velocity[2],2) + std::pow(velocity[3],2));
          
    for(int j=0;j<=3;j++) {
      xStart[j] = pIn[i].x_in().comp(j) ;
      //+ velocity[j]*delT;
    }
	  
    length = 5.0*fmToGeVinv; /// length in GeV-1 will have to changed for hydro
          
    //cout << " Length = " << length << endl;
          
    zeta = ( xStart[0] + std::sqrt( xStart[1]*xStart[1] + xStart[2]*xStart[2] + xStart[3]*xStart[3] )  )/std::sqrt(2);
    int pid = pIn[i].pid();
          
    if (pIn[i].form_time()<0.0) /// A parton without a virtuality or formation time, must set...
      {
	iSplit = 0;
	if (pIn[i].pid()==gid)
	  {
	    DEBUG << " parton is a gluon ";
	    iSplit = 1;
	  }
	else {
	  DEBUG << " parton is a quark ";
	}
	    
	tQ2 = generate_vac_t(pIn[i].pid(), pIn[i].nu(), QS/2.0, pIn[i].e()*pIn[i].e() ,zeta , iSplit);
            	    
	// KK:
	pIn[i].set_jet_v(velocity);              
	pIn[i].set_t(tQ2); // Also resets momentum!
	pIn[i].set_mean_form_time();
	double ft = generate_L (pIn[i].mean_form_time() ) ;              
	pIn[i].set_form_time(ft);
            
	//DEBUG:
	DEBUG ;
	DEBUG << " ***************************************************************************** " ;              
	DEBUG << " *  New generated virtuality = " << tQ2 << " Mean formation time = " << pIn[i].mean_form_time()/fmToGeVinv;              
	DEBUG << " *  set new formation time to " << pIn[i].form_time()/fmToGeVinv ;
	DEBUG << " * Maximum allowed virtuality = " << pIn[i].e()*pIn[i].e() << "   Minimum Virtuality = " << QS;
	DEBUG ;
	DEBUG << " * Jet velocity = " << pIn[i].jet_v().comp(0) << " " << pIn[i].jet_v().comp(1) << "  " << pIn[i].jet_v().comp(2) << "  " << pIn[i].jet_v().comp(3);
	DEBUG << " * reset location of parton formation = "<< pIn[i].x_in().t() << "  " << pIn[i].x_in().x() << "  " << pIn[i].x_in().y() << "  " << pIn[i].x_in().z() ;
	DEBUG << " ***************************************************************************** " ;
	DEBUG ;
	// end DEBUG:              
      }
          
          
          
    if (pIn[i].t() > QS + rounding_error) { // 
      double decayTime = pIn[i].mean_form_time()  ;
	    
      //cout << "  deltaT = " << deltaT <<endl;
      //cout << " parton origin time = " << pIn[i].x_in().t()/fmToGeVinv << " parton formation time = " << pIn[i].form_time()/fmToGeVinv << endl ;          
      // cout << " parton id " << pIn[i].pid() << " parton virtuality = " << pIn[i].t() << endl;          
      // cout << " parton momentum " << pIn[i].p_in().t() << "  " << pIn[i].p_in().x() << "  " << pIn[i].p_in().y() << "  " << pIn[i].p_in().z() << endl ;
	    
      double splitTime = pIn[i].form_time() + pIn[i].x_in().t() ;          
      //cout << " splitTime = " << splitTime/fmToGeVinv << endl ;
	    
      if (splitTime<Time)  {
	// do split
	//cout << " doing the split " << endl ;

	double t_used = pIn[i].t();
	if (t_used<QS)  t_used = QS;
              
	double tau_form = 2.0*pIn[i].nu()/t_used;              
	double z_low = QS/t_used/2.0 ;                  
	double z_hi = 1.0 - z_low;
              
	if (pid==gid) { // gluon
		
	  double val1 = P_z_gg_int(z_low, z_hi, zeta, t_used, tau_form,pIn[i].nu() );                
	  double val2 = nf*P_z_qq_int(z_low, z_hi, zeta, t_used, tau_form,pIn[i].nu() );
                    
	  if ( val1<0.0 || val2<0.0 ) {
	    cerr << " minus log of sudakov negative val1 , val2 = " << val1 << "  " << val2 << endl;
	    throw std::runtime_error("minus log of sudakov negative");
	    // cin >> blurb ;
	  }
                  
	  double ratio = val1/(val1+val2);                  
	  double r = ZeroOneDistribution(*get_mt19937_generator());
	  if (r>ratio) { // qqbar
		  
	    double r2 = ZeroOneDistribution(*get_mt19937_generator());
		  
	    // assign flavors
	    if (r2>0.6666){
	      pid_a = uid ;
	      pid_b = -1*uid ;
	    } else if ( r2 > 0.3333 && r2< 0.6666) {
	      pid_a = did ;
	      pid_b = -1*did ;
	    } else {
	      pid_a = sid;
	      pid_b = -1*sid;
	    }
                      
	    iSplit = 2;
	  } else { // gg
	    pid_a = gid ;
	    pid_b = gid ;
	    iSplit = 1;
	  }
	} else { // we had a quark
	  pid_a = pid ;		
	  pid_b = gid ;                  
	  iSplit = 0;
	}

	z = generate_vac_z(pid,QS/2.0,pIn[i].t(),zeta,pIn[i].nu(),iSplit) ;
	//cout << " generated z = " << z << endl;
        
	int iSplit_a = 0;                  
	if (pid_a==gid) iSplit_a = 1;                  
	
	// daughter virtualities
	double tQd1 = QS;                  
	if (z*z*pIn[i].t()>QS) {
	  tQd1 = generate_vac_t(pid_a, z*pIn[i].nu(), QS/2.0, z*z*pIn[i].t() ,zeta+std::sqrt(2)*pIn[i].form_time() , iSplit_a);
	}
        
	int iSplit_b = 0;                  
	if (pid_b==gid) iSplit_b = 1;
        
	double tQd2 = QS;                  
	if ((1.0-z)*(1.0-z)*pIn[i].t()>QS) {
	  tQd2 = generate_vac_t(pid_b, (1.0-z)*pIn[i].nu(), QS/2.0, (1.0-z)*(1.0-z)*pIn[i].t() ,zeta+std::sqrt(2)*pIn[i].form_time() , iSplit_b);
	}
	
	// opening angle
	double angle = generate_angle();                  
	
	double l_perp2 =  pIn[i].t()*z*(1.0 - z) - tQd2*z - tQd1*(1.0-z) ; ///< the transverse momentum squared
	// cout << " l_perp2 = " << l_perp2 << endl ;                   
	if (l_perp2<0.0) l_perp2 = 0.0; ///< test if negative
                  
	double l_perp = std::sqrt(l_perp2); ///< the momentum transverse to the parent parton direction
                  
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
                  
	// cout << " E, plong of d1 , E^2 - plong^2 - k_perp1^2 = " << energy << " " << plong << "  " << energy*energy - plong*plong - k_perp1_2 << endl;


	double newp[4];
	newp[0] = energy;
	newp[1] = plong*s_t*c_p + k_perp1[1];
	newp[2] = plong*s_t*s_p + k_perp1[2];
	newp[3] = plong*c_t + k_perp1[3];
                  
	// cout << " d1 momentum " << newp[0] << "  " << newp[1] << "  " << newp[2] << "  " << newp[3] << endl ;
	// cout << " d1 mass^2 = " << pow(newp[0],2) - pow(newp[1],2) - pow(newp[2],2) - pow(newp[3],2) << endl;

	double newx[4];                  
	newx[0] = Time + deltaTime ;                  
	for (int j=1;j<=3;j++) {
	  newx[j] = pIn[i].x_in().comp(j) + (Time + deltaTime - pIn[i].x_in().comp(0) )*velocity[j]/velocityMod;
	}
                  
	pOut.push_back(Parton(0,pid_a,0,newp,newx ));
	int iout = pOut.size()-1 ;
                  
	//cout << "  created a new parton from split with iout = " << iout << endl;
	pOut[iout].set_jet_v(velocity);
	// pOut[iout].set_t(tQd1); // KK: Not necessary, and in fact wrong
	pOut[iout].set_mean_form_time();
	double ft = generate_L (pOut[iout].mean_form_time());
	pOut[iout].set_form_time(ft);
		  
	// Second daughter
	double k_perp2[4];
	k_perp2[0] = 0.0;                  
	k_perp2[1] = (1.0-z)*(pIn[i].px() - pIn[i].pl()*s_t*c_p) - l_perp*std::cos(angle)*c_t*c_p + l_perp*std::sin(angle)*s_p ;                  
	k_perp2[2] = (1.0-z)*(pIn[i].py() - pIn[i].pl()*s_t*s_p) - l_perp*std::cos(angle)*c_t*s_p - l_perp*std::sin(angle)*c_p ;                  
	k_perp2[3] = (1.0-z)*(pIn[i].pz() - pIn[i].pl()*c_t) + l_perp*std::cos(angle)*s_t ;
	double k_perp2_2 = pow(k_perp2[1],2)+pow(k_perp2[2],2)+pow(k_perp2[3],2);

	energy = ( (1.0-z)*pIn[i].nu() + (tQd2 + k_perp2_2)/( 2.0*(1.0-z)*pIn[i].nu() ) )/std::sqrt(2.0) ;
	plong =  ( (1.0-z)*pIn[i].nu() - (tQd2 + k_perp2_2)/( 2.0*(1.0-z)*pIn[i].nu() ) )/std::sqrt(2.0) ;

	// cout << " E, plong of d2 , E^2 - plong^2 - k_perp^2 = " << energy << " " << plong << "  " << energy*energy - plong*plong - k_perp2_2 << endl;	  
	parent_perp = std::sqrt( pow(pIn[i].p(1),2) + pow(pIn[i].p(2),2) + pow(pIn[i].p(3),2) - pow(pIn[i].pl(),2) );                  
	mod_jet_v = std::sqrt( pow(pIn[i].jet_v().x(),2) +  pow(pIn[i].jet_v().y(),2) + pow(pIn[i].jet_v().z(),2) ) ;
                  
	newp[0] = energy;
	newp[1] = plong*s_t*c_p + k_perp2[1] ;
	newp[2] = plong*s_t*s_p + k_perp2[2] ;
	newp[3] = plong*c_t + k_perp2[3] ;
              
	// cout << " d2 momentum " << newp[0] << "  " << newp[1] << "  " << newp[2] << "  " << newp[3] << endl ;	      
	// cout << " d2 mass^2 = " << pow(newp[0],2) - pow(newp[1],2) - pow(newp[2],2) - pow(newp[3],2) << endl;

	newx[0] = Time + deltaTime;
	for (int j=1;j<=3;j++){
	  newx[j] = pIn[i].x_in().comp(j) + (Time + deltaTime - pIn[i].x_in().comp(0) )*velocity[j]/velocityMod;
	}
	      
	pOut.push_back(Parton(0,pid_b,0,newp,newx ));
              
	iout = pOut.size()-1 ;                  
	//cout << "  created a new parton from split with iout = " << iout << endl;
	// cout << endl;
	pOut[iout].set_jet_v(velocity);
	// pOut[iout].set_t(tQd2);	// KK: Not necessary, and in fact wrong
	pOut[iout].set_mean_form_time();                  
	ft = generate_L (pOut[iout].mean_form_time());	      
	pOut[iout].set_form_time(ft);

                  
      }  else  { // not time to split yet
	// pOut.push_back(pIn[i]);
      }
    } else { // virtuality too low
      // pOut.push_back(pIn[i]);
    }
          
    //          cout << " exit formation time = " << pIn[0].form_time() << endl;
	  
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
    if ( itcounter++ > 20 ) throw std::runtime_error("Stuck in endless loop") ;
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
      
        
    if (diff>0.0) {
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
    
  qL = m_fac*qhat*0.6*tau*profile(loc_e+tau) ;
    
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
    
  qL = m_fac*qhat*0.6*tau*profile(loc_e + tau) ;
    
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
    
  qL = qhat*0.6*tau*profile(loc_e + tau) ;
    
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
    
  qL = qhat*0.6*tau*profile(loc_e + tau) ;
    
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
    
  qL = qhat*0.6*tau*profile(loc_e + tau) ;
    
    
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
    
  qL = qhat*0.6*tau*profile(loc_e + tau) ;
    
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



// obsolete in the future ...
// ----------------------

Martini::Martini()
{
  SetId("Martini");
  VERBOSE(8);
}

Martini::~Martini()
{
  VERBOSE(8);

}

void Martini::Init()
{
  INFO<<"Intialize Martini ...";

  ZeroOneDistribution = uniform_real_distribution<double> { 0.0, 1.0 };
}

//void Martini::DoEnergyLoss(double deltaT, double Q2, const vector<Parton>& pIn, vector<Parton>& pOut)
void Martini::DoEnergyLoss(double deltaT, double Q2, vector<Parton>& pIn, vector<Parton>& pOut)
{
  // *my_file_ << "I'm Martini! my id= " << get_my_task_number() << "  " << ZeroOneDistribution( *get_mt19937_generator()) << endl;
  // return;

  if (Q2<=QS)
    VERBOSESHOWER(8)<< MAGENTA << "SentInPartons Signal received : "<<deltaT<<" "<<Q2<<" "<<&pIn;
}

