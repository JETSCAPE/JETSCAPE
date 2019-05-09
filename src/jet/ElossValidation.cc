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
  JSINFO<<"Intialize ElossValidate ...";

  std::string s = GetXMLElementText({"Eloss", "ElossValidate", "name"});
  JSINFO << s << " to be initializied ...";
  
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
    
  JSINFO << " the time in fm is " << time << " The time in GeV-1 is " << Time ;
  JSINFO << " pid = " << pIn[0].pid() << " E = " << pIn[0].e() << " px = " << pIn[0].p(1) << " py = " << pIn[0].p(2) << "  pz = " << pIn[0].p(3) << " virtuality = " << pIn[0].t() << " form_time in fm = " << pIn[0].form_time()/fmToGeVinv ;
  JSINFO << " color = " << pIn[0].color() << " anti-color = " << pIn[0].anti_color();

  
  for (int i=0;i<pIn.size();i++)
    {
      TakeResponsibilityFor ( pIn[i] ); // Generate error if another module already has responsibility.
      JSINFO << " Parton Q2= " << pIn[i].t();
      JSINFO << " Parton Id= " << pIn[i].pid() << " and mass= " << pIn[i].restmass();
      
      double velocity[4];
      velocity[0] = 1.0;
      for(int j=1;j<=3;j++){
	velocity[j] = pIn[i].p(j)/pIn[i].e();
      }
      
      if (pIn[i].form_time()<0.0) { /// A parton without a virtuality or formation time, must set...
	pIn[i].set_jet_v(velocity);
	pIn[i].set_t(QS*2.); 
	pIn[i].set_mean_form_time();
	JSINFO << " UPDATED Parton Q2= " << pIn[i].t();
      }	
	
      if (pIn[i].t()>=QS ) {
	JSINFO << " ************ ";
	JSINFO << " DOING ELOSS  ";
	JSINFO << " ************ ";

	
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
