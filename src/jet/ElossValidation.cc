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

#include "ElossValidation.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include <string>

#include "tinyxml2.h"
#include<iostream>

#include "FluidDynamics.h"

#define MAGENTA "\033[35m"

using namespace Jetscape;
using namespace std;

const double QS = 1.0 ;

ElossValidate::ElossValidate()
{
  SetId("ElossValidate");
  VERBOSE(8);
}

ElossValidate::~ElossValidate()
{
  VERBOSE(8);
}

void ElossValidate::Init()
{
  INFO<<"Intialize ElossValidate ...";

  // Redundant (get this from Base) quick fix here for now
  tinyxml2::XMLElement *eloss= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Eloss" );
  if ( !eloss ) {
    WARN << "Couldn't find tag Eloss";
    throw std::runtime_error ("Couldn't find tag Eloss");    
  }
  tinyxml2::XMLElement *xmle=eloss->FirstChildElement("ElossValidate");
  if ( !xmle ) {
    WARN << "Couldn't find tag Eloss -> ElossValidate";
    throw std::runtime_error ("Couldn't find tag Eloss -> ElossValidate");    
  }

  if ( !xmle->FirstChildElement( "name" ) ){
    WARN << "Couldn't find tag Eloss -> ElossValidate";
    throw std::runtime_error ("Couldn't find tag Eloss -> ElossValidate -> name");    
  }
  string s = xmle->FirstChildElement( "name" )->GetText();
  INFO << s << " to be initializied ...";
  
}


void ElossValidate::WriteTask(weak_ptr<JetScapeWriter> w)
{
  VERBOSE(8);
  auto f = w.lock();
  if ( !f ) return;
  f->WriteComment("ElossModule Parton List: "+GetId());
}

void ElossValidate::DoEnergyLoss(double deltaT, double time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut)
{

  VERBOSESHOWER(8)<< MAGENTA << "SentInPartons Signal received : "<<deltaT<<" "<<Q2<<" "<<&pIn;

  // Check hydro communication
  std::unique_ptr<FluidCellInfo> check_fluid_info_ptr;
  GetHydroCellSignal(1, 1.0, 1.0, 0.0, check_fluid_info_ptr);

  double delT = deltaT;
  double Time = time*fmToGeVinv;
  double deltaTime = delT*fmToGeVinv;
    
  INFO << " the time in fm is " << time << " The time in GeV-1 is " << Time ;
  INFO << " pid = " << pIn[0].pid() << " E = " << pIn[0].e() << " px = " << pIn[0].p(1) << " py = " << pIn[0].p(2) << "  pz = " << pIn[0].p(3) << " virtuality = " << pIn[0].t() << " form_time in fm = " << pIn[0].form_time()/fmToGeVinv ;
  INFO << " color = " << pIn[0].color() << " anti-color = " << pIn[0].anti_color();

  
  for (int i=0;i<pIn.size();i++)
    {
      TakeResponsibilityFor ( pIn[i] ); // Generate error if another module already has responsibility.
      INFO << " Parton Q2= " << pIn[i].t();
      INFO << " Parton Id= " << pIn[i].pid() << " and mass= " << pIn[i].restmass();
      
      double velocity[4];
      velocity[0] = 1.0;
      for(int j=1;j<=3;j++){
	velocity[j] = pIn[i].p(j)/pIn[i].e();
      }
      
      if (pIn[i].form_time()<0.0) { /// A parton without a virtuality or formation time, must set...
	pIn[i].set_jet_v(velocity);
	pIn[i].set_t(QS*2.); 
	pIn[i].set_mean_form_time();
	INFO << " UPDATED Parton Q2= " << pIn[i].t();
      }	
	
      if (pIn[i].t()>=QS ) {
	INFO << " ************ ";
	INFO << " DOING ELOSS  ";
	INFO << " ************ ";

	
	// Split once
	// ----------
	
	// daughter virtualities
	double tQd1 = QS-0.00001;
	double tQd2 = QS-0.00001;

	// Always just radiate a gluon
	int d1_pid = pIn[i].pid();
	int d2_pid = gid;

	double newp[4];
	double newx[4];

	double z1 = 0.6;
	double z2 = 1.0-z1;

	newx[0] = Time + deltaTime ;
	for (int j=1;j<=3;j++) {
	  newx[j] = pIn[i].x_in().comp(j) + (Time + deltaTime - pIn[i].x_in().comp(0) );
	}
	
	// First daughter
	newp[0] = z1 * pIn[i].e();
	newp[1] = z1 * pIn[i].px() * ( 1.1 );
	newp[2] = z1 * pIn[i].py() * ( 0.9 );
	newp[3] = z1 * pIn[i].pz() * ( 1.1 );
	pOut.push_back(Parton(0,d1_pid,0,newp,newx ));
	pOut.back().set_jet_v(pIn[i].jet_v());
	pOut.back().set_t(tQd1);
	pOut.back().set_mean_form_time();
	pOut.back().set_form_time(10);
	      
	      
	// Second daughter
	newp[0] = z2 * pIn[i].e();
	newp[1] = z2 * pIn[i].px() * ( 0.9 );
	newp[2] = z2 * pIn[i].py() * ( 1.1 );
	newp[3] = z2 * pIn[i].pz() * ( 0.8 );
	pOut.push_back(Parton(0,d2_pid,0,newp,newx ));
	pOut.back().set_jet_v(pIn[i].jet_v());
	pOut.back().set_t(tQd2);
	pOut.back().set_mean_form_time();
	pOut.back().set_form_time(10);
      }      
    }
       
}
