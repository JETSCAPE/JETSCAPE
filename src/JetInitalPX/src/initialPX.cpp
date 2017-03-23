#include <fstream>
#include <math.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

#include "initialPX.hpp"

void initialPX::next() {

      double pT, rapidity, phi;
      double max_pT = 100.0;
      double eta_cut = 1.0;
      double tempRand;
      const double maxN = 1.0*RAND_MAX;
      const double PI = 3.1415926;

      double parID,ppx,ppy,ppz,pp0,prx,pry,prz,pr0,mass;


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
      pT = max_pT*(rand()/maxN);

      phi = 2.0*PI*(rand()/maxN);
      rapidity=2.0*eta_cut*(rand()/maxN)-eta_cut;

      ppx = pT*cos(phi);
      ppy = pT*sin(phi);
      ppz = sqrt(pT*pT+mass*mass)*sinh(rapidity);
      pp0 = sqrt(pT*pT+mass*mass)*cosh(rapidity);

      prx = 0.0;
      pry = 0.0;
      prz = 0.0;
      pr0 = 0.0;

      set_px(ppx);     
      set_py(ppy);     
      set_pz(ppz);     
      set_e(pp0);     
      set_m(mass);     
      set_id(parID);     
      set_x(prx);     
      set_y(pry);     
      set_z(prz);     
      set_t(pr0);     

}

