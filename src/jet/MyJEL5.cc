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
/*This is the altered MyJEL5 code for checking the reproduction of qhat plot for multiple scattering*/
#include "MyJEL5.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"
#include <string>
#include <fstream>
#include "tinyxml2.h"
#include <iostream>

//#include "scat.h"

#include "FluidDynamics.h"

#define MAGENTA "\033[35m"

using namespace Jetscape;
using namespace std;

const double QS = 1.0;

void trans(double v[4], double p[4]) {
  double vv = sqrt(v[1] * v[1] + v[2] * v[2] + v[3] * v[3]);
  double ga = 1.0 / sqrt(1.0 - vv * vv);
  double ppar = p[1] * v[1] + p[2] * v[2] + p[3] * v[3];
  double gavv = (ppar * ga / (1.0 + ga) - p[0]) * ga;
  p[0] = ga * (p[0] - ppar);
  p[1] = p[1] + v[1] * gavv;
  p[2] = p[2] + v[2] * gavv;
  p[3] = p[3] + v[3] * gavv;
}

void transback(double v[4], double p[4]) {
  double vv = sqrt(v[1] * v[1] + v[2] * v[2] + v[3] * v[3]);
  double ga = 1.0 / sqrt(1.0 - vv * vv);
  double ppar = p[1] * v[1] + p[2] * v[2] + p[3] * v[3];
  double gavv = (-ppar * ga / (1.0 + ga) - p[0]) * ga;
  p[0] = ga * (p[0] + ppar);
  p[1] = p[1] - v[1] * gavv;
  p[2] = p[2] - v[2] * gavv;
  p[3] = p[3] - v[3] * gavv;
}



std::ofstream file1("qhat_verification_elastic_part_gluon_multiple1.dat");
RegisterJetScapeModule<MyJEL5> MyJEL5::reg("CustomModuleMyJEL5");
MyJEL5::MyJEL5() {
  SetId("MyJEL5");
  VERBOSE(8);
}

MyJEL5::~MyJEL5() { VERBOSE(8); file1.close();}

void MyJEL5::Init() {

  //file1.open("qhat_verification_elastic_part0.2_quark.dat");
  //JSINFO << "Initialize MyJEL5 ...";
  std::string s = GetXMLElementText({"Eloss", "CustomModuleMyJEL5", "name"});
  //upperLimitE = GetXMLElementDouble({"Eloss","CustomModuleMyJEL5","upperLimitE"});
  //upperLimitT = GetXMLElementDouble({"Eloss","CustomModuleMyJEL5","upperLimitT"});
  //elasticCollision.setter(0.2,0.2,0.1,102,0.5,0.5);
  elasticCollision.setter(0.2,0.2,0.1,110.0,0.5,0.5);
  //JSINFO << s << " to be initializied ...";
}

void MyJEL5::WriteTask(weak_ptr<JetScapeWriter> w) {
  VERBOSE(8);
  auto f = w.lock();
  if (!f)
    return;
  f->WriteComment("ElossModule Parton List: " + GetId());
}


    


void MyJEL5::DoEnergyLoss(double deltaT, double time, double Q2,
                                 vector<Parton> &pIn, vector<Parton> &pOut) {

  VERBOSESHOWER(8) << MAGENTA << "SentInPartons Signal received : " << deltaT
                   << " " << Q2 << " " << &pIn;
  ///Users/ritobandatta/OneDrive - Wayne State University/Elastic_Scattering/differentProcessRates_fixedAlphas.dat;
  /*
  std::ofstream file0("differentProcessRates_fixedAlphas1.dat");
  JSINFO<<"Hurray";
  file0<<"The file has the scattering rate of differnt process printed out at a fixed alpha_s of 0.3 and temperature 0.2GeV "<<"\n";
  file0<<"Energy(GeV) "<<"Rate (GeV)\n";
  for (int j=0;j<9;j++){
  for (int i=1;i<=100;i++){
    file0<<i*0.01<<" "<<elasticCollision.scattering_obj.get_rate(0,i-1,j)<<"\n";
  }
  file0<<"\n";
  }
  file0.close();
  }
  */

  // Check hydro communication
  std::unique_ptr<FluidCellInfo> check_fluid_info_ptr;
  GetHydroCellSignal(1, 1.0, 1.0, 0.0, check_fluid_info_ptr);

  double delT = deltaT;
  double Time = time * fmToGeVinv;
  double deltaTime = delT  * fmToGeVinv;
  //elasticCollision.sampler2();
  //elasticCollision.sampler3();

  //JSINFO << " the time in fm is " << time << " The time in GeV-1 is " << Time;
  //JSINFO << " color = " << pIn[0].color()
  //       << " anti-color = " << pIn[0].anti_color();
  //JSINFO<<time<<" "<<pIn[0].pstat();
	if (int(time)==8 and pIn[0].pstat()==1){
      //JSINFO<<"yooo "<<time;
      file1<<GetCurrentEvent()<<" "<<time<<" "<<0.2<<" "<<pow(pIn[0].p(1),2)+pow(pIn[0].p(2),2)<<"\n";
  }
  //	JSINFO<<"size is "<<pIn.size();
  for (int i = 0; i < pIn.size(); i++) {
    //JSINFO<<" "<<pIn.size();

    TakeResponsibilityFor(
        pIn[i]); // Generate error if another module already has responsibility.
    //JSINFO << " Parton Q2= " << pIn[i].t();
    //JSINFO << " Parton Id= " << pIn[i].pid()
    //       << " and mass= " << pIn[i].restmass();
    if (pIn[i].pstat()==-1){return;}

    if (pIn[i].form_time() <
        0.0) { /// A parton without a virtuality or formation time, must set...
      //JSINFO<<"yoo";
      pIn[i].set_t(QS * 2.);
      pIn[i].set_mean_form_time();
      pIn[i].set_form_time(pIn[i].mean_form_time());
      //JSINFO << " UPDATED Parton Q2= " << pIn[i].t();
      //pIn[i].reset_p(0,0,pIn[i].e());
      //elasticCollision.sampler2();
      //if (GetCurrentEvent()%100000==0){
      int multiple_coeff = (GetCurrentEvent()/1000)+1;
      pIn[i].reset_momentum(0,0,multiple_coeff*4, multiple_coeff*4);
      //pIn[i].reset_p(0,0,multiple_coeff*4);
      //JSINFO<<multiple_coeff*4;
    }
    
    int pid0 = pIn[i].pid();  //ivan
    int pid2 = -999;  //ivan
    int pid3 = -999; //ivan

    double vc0[4]={0.0};
    double pc0[4]={pIn[i].e(),pIn[i].px(),pIn[i].py(),pIn[i].pz()};
    double pc2[4] = {0.0}; // final recoil thermal parton //can be removed
    double pc3[4] = {0.0}; // initial thermal parton  //can be removed
    double newx[4]={0.0};
    //trans(vc0,pc0);
    double tempLoc = 0.2;
    /*
    if (int(time)%2==0){
      JSINFO<<"time is here "<<time<<" pid "<<pIn[i].pstat();
      pOut.push_back(Parton(0, pid0, 2, pc0, newx));
      pOut.push_back(pIn[i]);
    }
    else{
      JSINFO<<"time is there "<<time<<" pid "<<pIn[i].pstat();
       return;
    }
    */
    
    /*
    if (elasticCollision.elastic_kinematics(tempLoc,pid0,pid2,pid3,pc0,pc2,pc3)){ 
            //JSINFO<<"YEAH!";
            //transback(vc0,pc0);   //leading ivan
            //transback(vc0,pc2);   //recoil ivan
            //transback(vc0,pc3);   //hole ivan
            if (pc0[0] < pc2[0]) { //disable switch for heavy quark, only allow switch for identical particles
                double p0temp[4] = {0.0};
                for (int k = 0; k <= 3; k++) {
                    p0temp[k] = pc2[k];
                    pc2[k] = pc0[k];
                    pc0[k] = p0temp[k];
                    pid0=pid2;
                }
            //JSINFO<<"YEAH";
            }   
            //JSINFO<<" "<<pc0[0]<<" "<<pc0[1]<<" "<<pc0[2]<<" "<<pc0[3];
            //file1<<GetCurrentEvent()<<" "<<time<<" "<<0.2<<" "<<pow(pc0[1],2)+pow(pc0[2],2)<<"\n";
            pOut.push_back(Parton(0, pid0, 1, pc0, newx));
            int iout = pOut.size() - 1;
            pOut[iout].set_form_time(10000.0);
    }
    */
  }

 }

