//Parton Gun Test ...
#include "PGun.h"

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
  DEBUG<<"Initialize PGun Brick (Test) ...";
  VERBOSE(8);
  tinyxml2::XMLElement *pgun=GetHardXML()->FirstChildElement("PGun");

  if (pgun)
    {
      string s = pgun->FirstChildElement( "name" )->GetText();
      DEBUG << s << " to be initilizied ...";
      
      pgun->FirstChildElement("pT")->QueryDoubleText(&fixed_pT);

      DEBUG << s << " with fixed pT = "<<fixed_pT;
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
  INFO<<"Run Hard Process : "<<GetId()<< " ...";
  VERBOSE(8)<<"Current Event #"<<GetCurrentEvent();

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
       rapidity=2.0*eta_cut*(rand()/maxN)-eta_cut;
  
       p[0] = pT*cos(phi);
       p[1] = pT*sin(phi);
       p[2] = sqrt(pT*pT+mass*mass)*sinh(rapidity);
       p[3] = sqrt(pT*pT+mass*mass)*cosh(rapidity);
  
       AddParton(make_shared<Parton>(0,parID,0,0,p,xLoc));

       // DEBUG: (<< of Parton not working with Logger ... Check!)
       cout<<*GetPartonAt(i)<<endl;
     }
  
  VERBOSE(8)<<GetNHardPartons();
}
