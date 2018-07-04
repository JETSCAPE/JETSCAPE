// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

// quick and dirty test class for Eloss modules ...
// can be used as a user template ...


#ifndef LBTD_H
#define LBTD_H

#include <fstream>
#include <math.h>
#include "JetEnergyLossModule.h"
#include <boost/variant/variant.hpp>
#include "Rate.h"
#include <vector>

typedef Rate<2, 2, double(*)(const double, void*)> Rate22;
typedef Rate<3, 3, double(*)(const double*, void*)> Rate23;
typedef Rate<3, 4, double(*)(const double*, void*)> Rate32;
typedef boost::variant<Rate22, Rate23, Rate32> Process;

using namespace Jetscape;


class LBTD : public JetEnergyLossModule<LBTD>
{
 private:

 std::map<int, std::vector<Process> > AllProcesses;  

 public:

  LBTD();
  virtual ~LBTD();

  //main//
  void Init();
  void DoEnergyLoss(double deltaT, double Time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut);
  void WriteTask(weak_ptr<JetScapeWriter> w) {};
  int update_particle_momentum(double dt, double temp, std::vector<double> v3cell,
			int pid, double D_formation_t23, double D_formation_t32, fourvec incoming_p, std::vector<fourvec> & FS);
  void init_process(Process& r, std::string mode);
  int random_choose_Id();
};



#endif
