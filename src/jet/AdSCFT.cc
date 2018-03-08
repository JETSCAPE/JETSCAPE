// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

// quick and dirty test class implementation for Eloss modules ...
// can be used as a user template ...

#include "AdSCFT.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include <string>

#include "tinyxml2.h"
#include <iostream>
#include <fstream>

#include "fluid_dynamics.h"

#define MAGENTA "\033[35m"

const double QS = 1.0 ;

// quick and dirty fix ...
//default_random_engine generator;
//uniform_real_distribution<double> distribution(0.0,1.0);

AdSCFT::AdSCFT() 
{
  SetId("AdSCFT");
  kappa=-99;
  
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

      //double m_qhat=-99.99;
      adscft->FirstChildElement("kappa")->QueryDoubleText(&kappa);
      //SetQhat(m_qhat);

      INFO<<"AdSCFT kappa = "<<kappa;
     
      //JSDEBUG  << s << " with qhat = "<<GetQhat();      	
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
   w.lock()->WriteComment("ElossModule Parton List: "+GetId());
   w.lock()->WriteComment("Energy loss to be implemented accordingly ...");
}

void AdSCFT::DoEnergyLoss(double deltaT, double time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut)
{
  
    VERBOSESHOWER(8)<< MAGENTA << "SentInPartons Signal received : "<<deltaT<<" "<<Q2<<" "<<&pIn;

    for (int i=0;i<pIn.size();i++)
    {
      VERBOSE(8) << " Parton Q2= " << pIn[i].t();
      VERBOSE(8) << " Parton Id= " << pIn[i].pid() << " and mass= " << pIn[i].restmass();
      if (pIn[i].t()<=QS && pIn[i].e()>0.00001)
      {
	VERBOSE(7) << " ************ \n \n";
        VERBOSE(7) << " DOING ADSCFT \n \n";
        VERBOSE(7) << " ************ \n \n";

	//Parton 4-position
	double x[4];
	x[0]=pIn[i].x_in().x();
	x[1]=pIn[i].x_in().y();
	x[2]=pIn[i].x_in().z();
	x[3]=pIn[i].x_in().t();
        
	//Extract fluid properties
	std::unique_ptr<FluidCellInfo> check_fluid_info_ptr;

	double tau=std::sqrt(x[3]*x[3]-x[2]*x[2]);
	double temp, vx, vy, vz;
	if (tau>=0.4) {
	  //Extract fluid properties
    	  GetHydroCellSignal(x[3], x[0], x[1], x[2], check_fluid_info_ptr);
	  VERBOSE(8)<< MAGENTA<<"Temperature from Brick (Signal) = "<<check_fluid_info_ptr->temperature;

	  temp = check_fluid_info_ptr->temperature;
    	  vx = check_fluid_info_ptr->vx;
    	  vy = check_fluid_info_ptr->vy;
    	  vz = check_fluid_info_ptr->vz;

	  // cout << " temp= " << temp << " vx= " << vx << " vy= " << vy << " vz= " << vz << endl;
	}
	else temp=0., vx=0., vy=0., vz=0.;
    	
	vector<double> v;
    	v.push_back(vx), v.push_back(vy), v.push_back(vz), v.push_back(1.);
    	double v2=std::pow(v[0],2.)+std::pow(v[1],2.)+std::pow(v[2],2.);
    	double lore = 1./sqrt(1.-v2);     //Gamma Factor

        //Color Factor (wrt quark)
        double CF;
        if (pIn[i].pid()==21) CF=std::pow(9./4.,1./3.);	//Gluon
        else CF=1.; 					//Quark

        //Parton 4-momentum
        double p[4];
        p[0]=pIn[i].px();
	p[1]=pIn[i].py(); 
	p[2]=pIn[i].pz();
	p[3]=pIn[i].e();
        
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
	  // pIn[i].set_user_info( new AdSCFTUserInfo (ei,f_dist,l_dist));
	  fjcore::SharedPtr<fjcore::PseudoJet::UserInfoBase> info_shared(new AdSCFTUserInfo (ei,f_dist,l_dist));
	  pIn[i].user_info_shared_ptr() = info_shared;
	}

	// cout << " ei= " << ei << endl;
	// cout << " px= " << p[0] << " py= " << p[1] << " pz= " << p[2] << " en= " << p[3] << endl;
   	// cout << " x= " << x[0] << " y= " << x[1] << " z= " << x[2] << " t= " << x[3] << endl;
	
	//Parton 4-velocity
        vector<double> w;
        for (unsigned int j =0; j<4; j++) w.push_back(p[j]/p[3]);
        double w2=std::pow(w[0],2.)+std::pow(w[1],2.)+std::pow(w[2],2.);        
        // cout << " w2= " << w2 << endl;
	double virt=std::sqrt(p[3]*p[3]-w2*p[3]*p[3]-std::pow(pIn[i].restmass(),2.));
	//cout << " virt= " << virt << endl;
  	
        //Needed for boosts (v.w)
        double vscalw=v[0]*w[0]+v[1]*w[1]+v[2]*w[2];

        //Distance travelled in LAB frame
        l_dist+=deltaT;

        //Distance travelled in FRF - accumulating steps from previous, different fluid cells
        f_dist+=deltaT*std::sqrt(w2+lore*lore*(v2-2.*vscalw+vscalw*vscalw));
	 
        // cout << " l_dist= " << l_dist << " f_dist= " << f_dist << endl; 
        //Initial energy of the parton in the FRF
        double Efs=ei*lore*(1.-vscalw);	  

	double newEn=p[3];
        if (temp>=0.) newEn=p[3]-AdSCFT::Drag(f_dist, deltaT, Efs, temp, CF);
        if (newEn<0.) newEn=p[3]/10000000.;
        double lambda=newEn/p[3];
	// cout << " lambda = " << lambda << endl;
        // cout << " Elost= " << p[3]-newEn << endl;

	//Update 4-momentum
	for (unsigned a=0; a<4; a++) p[a]*=lambda;
	pIn[i].reset_momentum(p[0],p[1],p[2],p[3]);

	//Update 4-position
	for (unsigned a=0; a<4; a++) x[a]+=w[a]*deltaT;
	double fx[4];
	fx[0]=x[3], fx[1]=x[0], fx[2]=x[1], fx[3]=x[2];
	pIn[i].set_x(fx);
	//Feed into parton list
	// pIn.back().set_user_info(new AdSCFTUserInfo(ei,f_dist,l_dist));
	fjcore::SharedPtr<fjcore::PseudoJet::UserInfoBase> info_shared(new AdSCFTUserInfo (ei,f_dist,l_dist));
	pIn.back().user_info_shared_ptr() = info_shared;
	pOut.push_back(pIn[i]);

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
