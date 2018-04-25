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

  double tStart=0.6;		//Hydro starting time
  double T0;            	//End of quenching temperature
  bool in_vac;
 
 private:

  double kappa;

};

#endif // ADSCFT_H

