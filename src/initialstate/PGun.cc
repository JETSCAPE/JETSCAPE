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
//Parton Gun Test

#include "PGun.h"

using namespace Jetscape;

PGun::PGun() : HardProcess()
{
  fixed_E0=0;
  parId=1;
  SetId("PGun");
  VERBOSE(8);
}

PGun::~PGun()
{
  VERBOSE(8);
}

void PGun::InitTask()
{
  JSDEBUG<<"Initialize PGun Brick (Test) ...";
  VERBOSE(8);
  tinyxml2::XMLElement *pgun=GetHardXML()->FirstChildElement("PGun");

  if (pgun)
    {
      string s = pgun->FirstChildElement( "name" )->GetText();
      JSDEBUG << s << " to be initilizied ...";
      
      pgun->FirstChildElement("E0")->QueryDoubleText(&fixed_E0);
      pgun->FirstChildElement("parId")->QueryIntText(&parId);

      JSDEBUG << s << " with fixed E0 = "<<fixed_E0;
      INFO<<"Parton Gun with fixed E0 = "<<fixed_E0<< " and Id = "<< parId;
      
    }
  else
    {
      WARN << " : PGun not properly initialized in XML file ...";
      exit(-1);
    }
}
 
void PGun::Exec()
{
  VERBOSE(2)<<"Run Hard Process : "<<GetId()<< " ...";

  double p[4], xLoc[4];

  double E0, pT, rapidity, phi;
  double eta_cut = 1.0;
  double tempRand;
  const double maxN = 1.0*RAND_MAX;
  const double PI = 3.1415926;
  
  double ppx,ppy,ppz,pp0,mass; 

  for (int i=0;i<1;i++)
     {
       if (parId != 21) {
	 tempRand = rand()/maxN;
	 if(tempRand < 0.50) parId = -parId;
       }  

  if(std::abs(parId)==4)
  {
    mass=1.3;
  }        
  else if(std::abs(parId)==5)
  {
    mass=4.2;
  }  
  else
  {
    mass = 0.0;
  }

       E0=fixed_E0;
       pT = sqrt(E0*E0-mass*mass); //max_pT*(rand()/maxN);
     
       
       //phi = 2.0*PI*(rand()/maxN);
       rapidity=0;//2.0*eta_cut*(rand()/maxN)-eta_cut;
       phi = 0.0;
         
         
       p[1] = pT*cos(phi);
       p[2] = pT*sin(phi);
       p[3] = E0*sinh(rapidity);
       p[0] = E0*cosh(rapidity);
       //INFO<<"pgun p: "<<p[0]<<" "<<p[1]<<" "<<p[2]<<" "<<p[3];
  
       // Roll for a starting point
       // See: https://stackoverflow.com/questions/15039688/random-generator-from-vector-with-probability-distribution-in-c
       for (int i=0;i<=3; i++) {
	 xLoc[i] = 0.0;
       };
       
       if (!ini) {
	 INFO << "No initial state module, setting the starting location to 0. ";
       } else {
	 auto num_bin_coll = ini->GetNumOfBinaryCollisions();
	 if ( num_bin_coll.size()==0 ){
	   WARN << "num_of_binary_collisions is empty, setting the starting location to 0. Make sure to add e.g. trento before PythiaGun.";
	 } else {	 
	   std::discrete_distribution<> dist( begin(num_bin_coll),end(num_bin_coll) ); // Create the distribution
	   
	   // Now generate values
	   auto idx = dist( *GetMt19937Generator() );
	   auto coord = ini->CoordFromIdx( idx );
	   xLoc[1] = std::get<0>( coord );
	   xLoc[2] = std::get<1>( coord );
	 }
       }

       AddParton(make_shared<Parton>(0,parId,0,pT,rapidity,phi,p[0],xLoc));
  
       VERBOSEPARTON(7,*GetPartonAt(i))
	 <<" added "
	 <<" at x=" << xLoc[1]
	 <<", y=" << xLoc[2]
	 <<", z=" << xLoc[3];
     }
  
  VERBOSE(8)<<GetNHardPartons();
}
