/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * Modular, task-based framework
 * Intial Design: Joern Putschke, Kolja Kauder (Wayne State University)
 * For the full list of contributors see AUTHORS.

 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 * or via email to bugs.jetscape.org@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/
//Parton Gun Test

#include "PGun.h"

using namespace Jetscape;

PGun::PGun() : HardProcess()
{
  fixed_pT=0;
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
      
      pgun->FirstChildElement("pT")->QueryDoubleText(&fixed_pT);

      JSDEBUG << s << " with fixed pT = "<<fixed_pT;
      INFO<<"Parton Gun with fixed pT = "<<fixed_pT;
      
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
  for (int i=0;i<=3; i++) {
     xLoc[i] = 0.0;
   };

  double pT, rapidity, phi;
  double eta_cut = 1.0;
  double tempRand;
  const double maxN = 1.0*RAND_MAX;
  const double PI = 3.1415926;
  
  double parID,ppx,ppy,ppz,pp0,mass; 

  for (int i=0;i<2;i++)
     {
       tempRand = rand()/maxN;
       if(tempRand < 0.25) parID = 21;
       else if(tempRand < 0.50) parID = 1;
       else if(tempRand < 0.75) parID = 2;
       else parID = 3;
       if (parID != 21) {
	 tempRand = rand()/maxN;
	 if(tempRand < 0.50) parID = -parID;
       }            
       
       mass = 0.0;
       pT = fixed_pT; //max_pT*(rand()/maxN);
       
       phi = 2.0*PI*(rand()/maxN);
       rapidity=0;//2.0*eta_cut*(rand()/maxN)-eta_cut;
              
       p[1] = pT*cos(phi);
       p[2] = pT*sin(phi);
       p[3] = sqrt(pT*pT+mass*mass)*sinh(rapidity);
       p[0] = sqrt(pT*pT+mass*mass)*cosh(rapidity);
  
       //AddParton(make_shared<Parton>(0,parID,0,p,xLoc));
       AddParton(make_shared<Parton>(0,parID,0,pT,rapidity,phi,p[0],xLoc));

       // DEBUG: (<< of Parton not working with Logger VERBOSE standard ... Check!
       //JetScapeLogger::Instance()->VerboseParton(6,*GetPartonAt(i))<<__PRETTY_FUNCTION__<<" : ";
       VERBOSEPARTON(6,*GetPartonAt(i));
       //cout<<*GetPartonAt(i)<<endl;
     }
  
  VERBOSE(8)<<GetNHardPartons();
}
