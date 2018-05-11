/*************************************************************************************
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
*
* AdSCFT Module: Daniel Pablos (May 2017)
* Implementation of energy loss rate derived in JHEP1605(2016)098 [arXiv:1511.07567], 
* see also Phys.Rev.D90(2014)no.2,025033 [arXiv:1402.6756] 
*************************************************************************************/

#include "AdSCFT.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include <string>

#include "tinyxml2.h"
#include <iostream>
#include <fstream>

#include "FluidDynamics.h"

#define MAGENTA "\033[35m"

AdSCFT::AdSCFT() 
{
  SetId("AdSCFT");
  kappa=-99;
  T0=-99;
  Q0=-99; 
 
  VERBOSE(8);
}

AdSCFT::~AdSCFT()
{
  VERBOSE(8);
}

void AdSCFT::Init()
{
  INFO<<"Intialize AdSCFT ...";

  // Redundant (get this from Base) quick fix here for now
  tinyxml2::XMLElement *eloss= JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("Eloss" );  
  tinyxml2::XMLElement *adscft=eloss->FirstChildElement("AdSCFT");
  
  if (adscft)
    {   
      string s = adscft->FirstChildElement( "name" )->GetText();
      JSDEBUG << s << " to be initilizied ...";

      //Kappa
      adscft->FirstChildElement("kappa")->QueryDoubleText(&kappa);
      INFO<<"AdSCFT kappa = "<<kappa;
     
      //T0 [GeV]
      adscft->FirstChildElement("T0")->QueryDoubleText(&T0); 
      INFO << "AdSCFT T0 = " << T0;

      //Q0 [GeV]
      adscft->FirstChildElement("Q0")->QueryDoubleText(&Q0);
      INFO << "AdSCFT Q0 = " << Q0;

      //Vac or Med
      in_vac = false;
      int flagInt=-100;
      adscft->FirstChildElement("in_vac")->QueryIntText(&flagInt);
      in_vac = flagInt;
    }
  else
    {
      WARN << " : AdSCFT not properly initialized in XML file ...";
      exit(-1);
    }
}


void AdSCFT::WriteTask(weak_ptr<JetScapeWriter> w)
{
   VERBOSE(8);
   auto f = w.lock();
   if ( !f ) return;
   f->WriteComment("ElossModule Parton List: "+GetId());
   f->WriteComment("Energy loss to be implemented accordingly ...");
}

void AdSCFT::DoEnergyLoss(double deltaT, double time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut)
{
  
    VERBOSESHOWER(8)<< MAGENTA << "SentInPartons Signal received : "<<deltaT<<" "<<Q2<<" "<<&pIn;

    for (int i=0;i<pIn.size();i++)
    {
      JSDEBUG << " in AdS/CFT";
      JSDEBUG << " Parton Q2= " << pIn[i].t();
      JSDEBUG << " Parton Id= " << pIn[i].pid() << " and mass= " << pIn[i].restmass();

      //Parton 4-momentum
      double p[4];
      p[0]=pIn[i].px();
      p[1]=pIn[i].py();
      p[2]=pIn[i].pz();
      p[3]=pIn[i].e();

      //Parton velocity
      vector<double> w;
      for (unsigned int j =0; j<4; j++) w.push_back(p[j]/p[3]);
     
      //Parton 4-position
      double initR0=pIn[i].x_in().t();	//Time when the parton was last modified
      double x[4];
      x[0]=pIn[i].x_in().x()+(time-initR0)*w[0];
      x[1]=pIn[i].x_in().y()+(time-initR0)*w[1];
      x[2]=pIn[i].x_in().z()+(time-initR0)*w[2];
      x[3]=pIn[i].x_in().t()+(time-initR0)*w[3];
 
      //Extract fluid properties
      std::unique_ptr<FluidCellInfo> check_fluid_info_ptr;
      double tau=std::sqrt(x[3]*x[3]-x[2]*x[2]);
      double temp, vx, vy, vz;
      //Only get a temp!=0 if in_vac=0
      if (x[3]>=tStart && !in_vac)
      //if (tau>=tStart && !in_vac)  //Should use tau, not t, in absence of preequilibrium eloss
      {
        GetHydroCellSignal(x[3], x[0], x[1], x[2], check_fluid_info_ptr);
        VERBOSE(8)<< MAGENTA<<"Temperature from Brick (Signal) = "<<check_fluid_info_ptr->temperature;

        temp = check_fluid_info_ptr->temperature;
        vx = check_fluid_info_ptr->vx;
        vy = check_fluid_info_ptr->vy;
        vz = check_fluid_info_ptr->vz;
      }
      else temp=0., vx=0., vy=0., vz=0.;
      JSDEBUG << " system time= " << time << " parton time= " << x[3] << " tau= " << tau << " temp= " << temp << " vz= " << w[2];
      JSDEBUG << " energy= " << pIn[i].e(); 

      //Do eloss in AdS/CFT if:
      // *Virtuality below Qcut ( QS[GeV^2] , NOT taken from XML yet )
      // *Fluid temperature above Tcut ( T0 from XML )
      // *Parton is not completely quenched ( Ecut = 0.00001 )
      double QS=Q0*Q0;
      if (pIn[i].t()<=QS+rounding_error && temp>=T0 && pIn[i].e()>0.00001)
      {
        //cout << " ADS Q= " << pIn[i].t() << " Q0= " << Q0 << " temp= " << temp << " T0= " << T0 << endl;
	//cout << " ADS tau= " << tau << " x= " << x[0] << " y= " << x[1] << " z= " << x[2] << " t= " << x[3] << endl;
	TakeResponsibilityFor ( pIn[i] ); // Generate error if another module already has responsibility.

	VERBOSE(8) << " ************ \n \n";
        VERBOSE(8) << " DOING ADSCFT \n \n";
        VERBOSE(8) << " ************ \n \n";

	//Fluid quantities
        vector<double> v;
    	v.push_back(vx), v.push_back(vy), v.push_back(vz), v.push_back(1.);
    	double v2=std::pow(v[0],2.)+std::pow(v[1],2.)+std::pow(v[2],2.);
    	double lore = 1./sqrt(1.-v2);     //Gamma Factor

        //Color Factor (wrt quark)
        double CF;
        if (pIn[i].pid()==21) CF=std::pow(9./4.,1./3.);	//Gluon
        else CF=1.; 					//Quark

	//Energy of parton as it entered this module for the first time
	double ei=p[3];
	double l_dist=0., f_dist=0.;
	if (pIn[i].has_user_info<AdSCFTUserInfo>()) {
	  ei=pIn[i].user_info<AdSCFTUserInfo>().part_ei();
	  l_dist=pIn[i].user_info<AdSCFTUserInfo>().l_dist();
	  f_dist=pIn[i].user_info<AdSCFTUserInfo>().f_dist();
	}
	else
	{
	  pIn[i].set_user_info(new AdSCFTUserInfo(ei,f_dist,l_dist));
	}

	//JSDEBUG << " ei= " << ei;
	//JSDEBUG << " px= " << p[0] << " py= " << p[1] << " pz= " << p[2] << " en= " << p[3];
   	//JSDEBUG << " x= " << x[0] << " y= " << x[1] << " z= " << x[2] << " t= " << x[3];
	
        double w2=std::pow(w[0],2.)+std::pow(w[1],2.)+std::pow(w[2],2.);        
        //JSDEBUG << " w2= " << w2;
	double virt=std::sqrt(p[3]*p[3]-w2*p[3]*p[3]-std::pow(pIn[i].restmass(),2.));
	//JSDEBUG << " virt= " << virt;
  	
        //Needed for boosts (v.w)
        double vscalw=v[0]*w[0]+v[1]*w[1]+v[2]*w[2];

        //Distance travelled in LAB frame
        l_dist+=deltaT;

        //Distance travelled in FRF - accumulating steps from previous, different fluid cells
        f_dist+=deltaT*std::sqrt(w2+lore*lore*(v2-2.*vscalw+vscalw*vscalw));
	 
        //JSDEBUG << " l_dist= " << l_dist << " f_dist= " << f_dist; 
        //Initial energy of the parton in the FRF
        double Efs=ei*lore*(1.-vscalw);	  

	double newEn=p[3];
        if (temp>=0.) newEn=p[3]-AdSCFT::Drag(f_dist, deltaT, Efs, temp, CF);
        if (newEn<0.) newEn=p[3]/10000000.;
        double lambda=newEn/p[3];
	//JSDEBUG << " lambda= " << lambda;
        //JSDEBUG << " Elost= " << p[3]-newEn;

	//Update 4-momentum
	for (unsigned a=0; a<4; a++) p[a]*=lambda;
	pIn[i].reset_momentum(p[0],p[1],p[2],p[3]);

	//DON'T Update 4-position here, already done at beginning
	//for (unsigned a=0; a<4; a++) x[a]+=w[a]*deltaT;
	
        double fx[4];
	fx[0]=x[3], fx[1]=x[0], fx[2]=x[1], fx[3]=x[2];
	pIn[i].set_x(fx);
	//Feed into parton list
	pOut.push_back( pIn[i] );
	pOut.back().set_user_info(new AdSCFTUserInfo(ei,f_dist,l_dist));

	//Copy variables needed in case parton returns to MATTER in future steps
	double velocity_jet[4];
        velocity_jet[0]=1.0;
        velocity_jet[1]=pIn[i].jet_v().x();
        velocity_jet[2]=pIn[i].jet_v().y();
        velocity_jet[3]=pIn[i].jet_v().z();
        pOut[pOut.size()-1].set_jet_v(velocity_jet);
        pOut[pOut.size()-1].set_mean_form_time();
        double ft = pOut[pOut.size()-1].mean_form_time();
        pOut[pOut.size()-1].set_form_time(ft);

      } //End if do-eloss

    } //End pIn loop
}

double AdSCFT::Drag(double f_dist, double deltaT, double Efs, double temp, double CF)
{
  if (kappa!=0.)
  {
    double tstop=0.2*std::pow(Efs,1./3.)/(2.*std::pow(temp,4./3.)*kappa)/CF;	//Stopping distance in fm
    double beta=tstop/f_dist;						//Ratio of stopping distance to distance travelled
    if (beta>1.)		//If did not get to stopping distance
    {
      double intpiece=Efs*deltaT*4./(3.141592)*(1./(beta*tstop*std::sqrt(beta*beta-1.)));
      return intpiece;
    }
    else			//If reached stopping distance, subtract all energy
    {
      return 100000.;
    }
  }
  else return 0.;
}

void AdSCFT::Clear()
{

}
