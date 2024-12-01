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

#include <iostream>

#include "PartonPrinter.h"
#include "GTL/graph.h"
#include <GTL/edge_map.h>
#include <GTL/node_map.h>
#include "JetClass.h"
#include "JetScapeLogger.h"
#include <string>

using namespace std;

namespace Jetscape {

// Register the module with the base class
RegisterJetScapeModule<PartonPrinter> PartonPrinter::reg("PartonPrinter");

PartonPrinter::PartonPrinter() {}

PartonPrinter::~PartonPrinter() {}

void PartonPrinter::Init() {
  this->SetId("PartonPrinter");
    JSINFO << "Initialize PartonPrinter ...";
    string filename = GetXMLElementText({"PartonPrinter","FileName"});
    double time_from_xml = GetXMLElementDouble({"Eloss", "maxT"});
    set_time (time_from_xml);
    double delta_t = GetXMLElementDouble({"Eloss", "deltaT"});
    set_delta_time (delta_t);
    
  dist_output.open(filename);
}

void PartonPrinter::Exec() {
  VERBOSE(2) << "Run PartonPrinter: print shower from event # "
             << GetCurrentEvent() << " ...";
}

void PartonPrinter::GetFinalPartons(
    shared_ptr<PartonShower>
        pShower /*, vector<shared_ptr<Parton>>& fPartons*/) {
  if (pShower) {
    //vector<shared_ptr<Parton>> vPin;
    for (unsigned int ipart = 0;
         ipart < pShower.get()->GetFinalPartons().size(); ++ipart) {
      //fPartons.push_back( pShower.get()->GetFinalPartons().at(ipart));
      if (std::abs(pShower.get()->GetFinalPartons().at(ipart)->pid()) <50)
      {
          double px = pShower.get()->GetFinalPartons().at(ipart)->px();
          double py = pShower.get()->GetFinalPartons().at(ipart)->py();
          double pz = pShower.get()->GetFinalPartons().at(ipart)->pz();
          double m = pShower.get()->GetFinalPartons().at(ipart)->restmass();
          
          double e = std::sqrt( px*px + py*py + pz*pz + m*m);
          double vx = px/e;
          double vy = py/e;
          double vz = pz/e;
          
          double x = pShower.get()->GetFinalPartons().at(ipart)->x_in().comp(1) ;
          double y = pShower.get()->GetFinalPartons().at(ipart)->x_in().comp(2) ;
          double z = pShower.get()->GetFinalPartons().at(ipart)->x_in().comp(3) ;
          double t = pShower.get()->GetFinalPartons().at(ipart)->x_in().comp(0) ;
          
          int choice = 1;
          double write_time = this->time() + this->delta_time();
          
          if (t > write_time)
          {
              JSWARN << " write out time is before split time, t = " << t << " this-> time = " << write_time ;
              JSWARN << " propagate backward? yes 1, no 0 "  ;
              cin >> choice;
          }
        
          if (choice==0) continue;
          
          x = x + vx*(write_time - t);
          y = y + vy*(write_time - t);
          z = z + vz*(write_time - t);
          
          double mu2 = pShower.get()->GetFinalPartons().at(ipart)->t() ;
          double tau = pShower.get()->GetFinalPartons().at(ipart)->form_time();
          if (std::isinf(tau) )
          {
              JSINFO << " tau is infinite = " << tau ;
              double blurb;
              cin >> blurb;
          }
          
          if (tau<=0.0)
          {
              JSINFO << " tau is negative or zero = " << " for parton id = " << pShower.get()->GetFinalPartons().at(ipart)->pid() << " = " <<  tau ;
              JSINFO << BOLDRED << " px =  " << px << " py = " << py << " pz = " << pz ;
              
              tau = hbarC*2*std::sqrt( e*e + mu2)/mu2;
              
              JSINFO << " fixed tau = " << tau ;
          }
          
          double t_size = 0.632/sqrt(mu2)*(write_time - t)/tau ;// 0.632 = \sqrt{10}*0.2 [fm GeV]
          
        dist_output << ipart << " "
                    << pShower.get()->GetFinalPartons().at(ipart)->pid() << " "
                    << pShower.get()->GetFinalPartons().at(ipart)->e() << " "
                    << pShower.get()->GetFinalPartons().at(ipart)->px() << " "
                    << pShower.get()->GetFinalPartons().at(ipart)->py() << " "
                    << pShower.get()->GetFinalPartons().at(ipart)->pz() << " "
                    << y  << " "
                    << z  << " "
                    << t_size << endl;
      }

      //vPin.push_back( pShower.get()->GetFinalPartons().at(ipart));
    }
    //this->pFinals.push_back(vPin);
  }
}

void PartonPrinter::GetFinalPartons2(shared_ptr<PartonShower> pShower) {
  if (pShower) {
    /*
    for(unsigned int ipart=0; ipart <  pShower.get()->GetFinalPartons().size(); ++ipart)
    {
      this->pFinals.push_back( pShower.get()->GetFinalPartons());
      //cout << "############### FINAL PARTON IN THE VECTOR NUMBER : " << ipart << " is " << this->pFinals.at(ipart) << "###################\n";
    }
*/
    //this->pFinals.clear();
    this->pFinals.push_back(pShower.get()->GetFinalPartons());
  }
}

void PartonPrinter::Clear() {
  //dist_output << " *********************************** " << endl;
  //dist_output.close();
  this->pFinals.clear();
}

void PartonPrinter::GetPartonsAtTime(shared_ptr<PartonShower> pShower,
                                     vector<shared_ptr<Parton>> &fPartons,
                                     double time) {}

} // end namespace Jetscape
