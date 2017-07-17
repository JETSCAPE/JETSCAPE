// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

// quick and dirty test class implementation for Eloss modules ...
// can be used as a user template ...

#include "ElossModulesTest.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include <string>

#include "tinyxml2.h"
#include<iostream>

#include "fluid_dynamics.h"

#define MAGENTA "\033[35m"

// quick and dirty fix ...
default_random_engine generator;
uniform_real_distribution<double> distribution(0.0,1.0);

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
 
  if (matter)
    {   
      string s = matter->FirstChildElement( "name" )->GetText();
    
      DEBUG << s << " to be initilizied ...";

      double m_qhat=-99.99;
      matter->FirstChildElement("qhat")->QueryDoubleText(&m_qhat);
      SetQhat(m_qhat);
      
      DEBUG  << s << " with qhat = "<<GetQhat();      	
    }
  else
    {
      WARN << " : Matter not properly initialized in XML file ...";
      exit(-1);
    }
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
void Matter::DoEnergyLoss(double deltaT, double t, double Q2, vector<Parton>& pIn, vector<Parton>& pOut)
{
  double z=0.5;
  
  if (Q2>5)
    {
      VERBOSESHOWER(8)<< MAGENTA << "SentInPartons Signal received : "<<deltaT<<" "<<Q2<<" "<<&pIn;
      
      FluidCellInfo* check_fluid_info_ptr = new FluidCellInfo;
      GetHydroCellSignal(1, 1.0, 1.0, 0.0, check_fluid_info_ptr);      
      VERBOSE(8)<< MAGENTA<<"Temperature from Brick (Signal) = "<<check_fluid_info_ptr->temperature;            
      delete check_fluid_info_ptr;

      double rNum;

      //DEBUG:
      //cout<<" ---> "<<pIn.size()<<endl;
      for (int i=0;i<pIn.size();i++)
	{	  
	  rNum=distribution(generator);
	  //DEBUG:
	  //cout<<i<<" "<<rNum<<endl;

	  // simulate a "random" split 50/50 in pT
	  if (rNum>0.7)
	    {
	      //cout<<pIn[i];
	      double newPt=pIn[i].pt()*z;
	      double newPt2=pIn[i].pt()*(1-z);
	    
	      //pOut.push_back(Parton(0,21,0,newPt,pIn[i].eta(),pIn[i].phi()+0.1,newPt));
	      //pOut.push_back(Parton(0,21,0,newPt2,pIn[i].eta(),pIn[i].phi()-0.1,newPt));

	      //tes case, perfectly collinear ...
	      pOut.push_back(Parton(0,21,0,newPt,pIn[i].eta(),pIn[i].phi(),newPt));
	      pOut.push_back(Parton(0,21,0,newPt2,pIn[i].eta(),pIn[i].phi(),newPt));
	      // DEBUG: dirty ...	      
	      //cout<<pOut[pOut.size()-1];
	      //cout<<pOut[pOut.size()-2];
	    }
	
	cout << "**********HERE is Delta TTT******** " << deltaT << "\n";
	cout << "**********HERE is TTT******** " << t << "\n";
	  
	}
      
      // Add a new root node ... (dummy ...)
      // Ahh stupid declared as const orginally (removed for test ...)
      // Maybe better a seperate vector !? (TBD)
      if (rNum>0.9)
	{
	  // quick and dirty ...
	  pIn.push_back(Parton(0,21,0,1.5,0,pIn[0].phi(),1.5));
	  //DEBUG:
	  //cout<<pIn.size()<<endl;
	}
      
    }
}
// obsolete in the future ...
/*
void Matter::Exec()
{
   INFO<<"Run Matter ...";
   DEBUG<<"Qhat = "<<GetQhat();
   
   DEBUG<<"Emit Signal: jetSignal(10,20.3)";
   jetSignal(10,20.3);
   double edensity=-1;
   edensitySignal(1,edensity);
   DEBUG<< MAGENTA<<"Received edensity = "<<edensity<<" for t="<<1;  
   
   if (GetShowerInitiatingParton())
     {
       //cout<<shared_from_this().get()<<endl;
       //cout<< *GetShowerInitiatingParton()<<endl;
       VERBOSEPARTON(6,*GetShowerInitiatingParton());
       //PrintShowerInitiatingParton();
     }
   
   FluidCellInfo* check_fluid_info_ptr = new FluidCellInfo;
   GetHydroCellSignal(1, 1.0, 1.0, 0.0, check_fluid_info_ptr);
   
   DEBUG<< MAGENTA<<"Temperature from Brick (Signal) = "<<check_fluid_info_ptr->temperature;
   
   //check_fluid_info_ptr->Print();
   
   delete check_fluid_info_ptr;
   
}
*/

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
}

//void Martini::DoEnergyLoss(double deltaT, double Q2, const vector<Parton>& pIn, vector<Parton>& pOut)
void Martini::DoEnergyLoss(double deltaT, double Q2, vector<Parton>& pIn, vector<Parton>& pOut)
{
  if (Q2<5)
    VERBOSESHOWER(8)<< MAGENTA << "SentInPartons Signal received : "<<deltaT<<" "<<Q2<<" "<<&pIn;
}

/*
void Martini::Exec()
{
  INFO<<"Run Martini ...";
  DEBUG<<"Qhat = "<<GetQhat();
  //DEBUG<<"Emit Signal: jetSignal(100,200.3)";
  //cout<<jetSignal.is_connected()<<endl;
  //jetSignal(100,200.3);

  //cout<<shared_from_this().get()<<endl;
  // Logger not working with overloaded << from Parton class ...
  // check and resolve ...
   if (GetShowerInitiatingParton())
     {
       cout<< *GetShowerInitiatingParton()<<endl;
       //PrintShowerInitiatingParton();
     }
   
   FluidCellInfo* check_fluid_info_ptr = new FluidCellInfo;
   GetHydroCellSignal(1, 1.0, 1.0, 0.0, check_fluid_info_ptr);  
   DEBUG<< MAGENTA<<"Temperature from Brick (Signal) = "<<check_fluid_info_ptr->temperature;
   
   //check_fluid_info_ptr->Print();
   
   delete check_fluid_info_ptr;
}
*/
