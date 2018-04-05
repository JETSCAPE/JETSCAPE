// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
//
// AdSCFT Module: Daniel Pablos (May 2017)
//
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...
//
// First test to include AdSCFT drag in JetScape ...
//

#ifndef ADSCFT_H
#define ADSCFT_H

#include "JetEnergyLossModule.h"
using namespace Jetscape;


class AdSCFTUserInfo: public fjcore::PseudoJet::UserInfoBase
{
 public :
 AdSCFTUserInfo(double ei, double f_dist, double l_dist) : _part_ei(ei), _f_dist(f_dist), _l_dist(l_dist){};
  double part_ei() const { return _part_ei; }  
  double f_dist() const { return _f_dist; }  
  double l_dist() const { return _l_dist; }
  double _part_ei;
  double _f_dist;
  double _l_dist;
  ~AdSCFTUserInfo(){};

};

class AdSCFT : public JetEnergyLossModule<AdSCFT> 
{  
 public:
  
  AdSCFT();
  virtual ~AdSCFT();

  void Init();
  void Clear();
  
  void DoEnergyLoss(double deltaT, double time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut);
  double Drag(double f_dist, double deltaT, double Efs, double temp, double CF);
  void WriteTask(weak_ptr<JetScapeWriter> w);
 
 private:

  double kappa;

};

#endif // ADSCFT_H

