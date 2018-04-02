/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * For the full list of contributors see AUTHORS.
 * Report issues at https://github.com/amajumder/JETSCAPE-COMP/issues
 * or via email to bugs.jetscape.org@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

// First test to include AdSCFT drag in JetScape

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

