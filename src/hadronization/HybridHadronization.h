/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2019
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

#ifndef HYBRIDHADRONIZATION_H
#define HYBRIDHADRONIZATION_H

#include "HadronizationModule.h"
#include "JetScapeLogger.h"
#include "Pythia8/Pythia.h"

#include <cmath>
#include <random>
#include <vector>

using namespace Jetscape;

class HybridHadronization : public HadronizationModule<HybridHadronization>
{  
 public:

  HybridHadronization();
  virtual ~HybridHadronization();
  
  void Init();
  void DoHadronization(vector<vector<shared_ptr<Parton>>>& shower, vector<shared_ptr<Hadron>>& hOut, vector<shared_ptr<Parton>>& pOut);
  void WriteTask(weak_ptr<JetScapeWriter> w);

 private:
  // Allows the registration of the module so that it is available to be used by the Jetscape framework.
  static RegisterJetScapeModule<HybridHadronization> reg;
  
  double SigM2_calc(double R2chg, double qm1, double qm2, double qq1, double qq2);
  double SigBR2_calc(double R2chg, double qm1, double qm2, double qm3, double qq1, double qq2, double qq3);
  double SigBL2_calc(double SigBR2, double qm1, double qm2, double qm3);
  
  //double sigma_pi, sigma_k, sigma_nuc, maxE_level, gmax, xmq, xms, hbarc, dist2cut, sh_recofactor, th_recofactor, SigRB, SigLB;
  double maxM_level, maxB_level, gmax, xmq, xms, xmc, xmb, hbarc, dist2cut, sh_recofactor, th_recofactor, p_fake, had_prop, part_prop;
  int number_p_fake;
  double SigNucR2,SigNucL2,SigOmgR2,SigOmgL2,SigXiR2,SigXiL2,SigSigR2,SigSigL2,
	SigOcccR2,SigOcccL2,SigOccR2,SigOccL2,SigXiccR2,SigXiccL2,SigOcR2,SigOcL2,SigXicR2,SigXicL2,SigSigcR2,SigSigcL2,
	SigObbbR2,SigObbbL2,SigObbcR2,SigObbcL2,SigObbR2,SigObbL2,SigXibbR2,SigXibbL2,
	SigObccR2,SigObccL2,SigObcR2,SigObcL2,SigXibcR2,SigXibcL2,SigObR2,SigObL2,SigXibR2,SigXibL2,SigSigbR2,SigSigbL2;
  double SigPi2, SigPhi2, SigK2, SigJpi2, SigDs2, SigD2, SigUps2, SigBc2, SigB2;
  const double pi = 3.1415926535897932384626433832795;
  int attempts_max;
  unsigned int rand_seed;
  int reco_hadrons_pythia;
  bool goldstonereco;
  std::vector<int> IDs;
  
  //variables for recombination color structure
  vector<vector<vector<int>>> Tempjunctions; // vector of all tempjunctions
  vector<vector<int>> JunctionInfo; // vector of one junction's color tag and particle info
  vector<int> IdColInfo1; vector<int> IdColInfo2; vector<int> IdColInfo3; vector<int> IdColInfo4;
  
  std::mt19937_64 eng; //RNG - Mersenne Twist - 64 bit
  double ran();
  
  //4-vector boost
  static FourVector HHboost(FourVector B, FourVector vec_in){
	double xlam[4][4], beta2;
		beta2 = B.x()*B.x() + B.y()*B.y() + B.z()*B.z();
		B.Set(B.x(),B.y(),B.z(),1./(sqrt(1. - beta2)));
	
	double beta2inv;
	if(beta2 > 0.){beta2inv = 1./beta2;}
	else{beta2inv = 0.;}

	xlam[0][0] = B.t();
	xlam[0][1] = -B.t()*B.x();
	xlam[0][2] = -B.t()* B.y();
	xlam[0][3] = -B.t()*B.z();
	xlam[1][0] = xlam[0][1];
	xlam[1][1] = 1.+(B.t() - 1.)*(B.x()*B.x())*beta2inv;
	xlam[1][2] = (B.t()-1.)*B.x()*B.y()*beta2inv;
	xlam[1][3] = (B.t()-1.)*B.x()*B.z()*beta2inv;
	xlam[2][0] = xlam[0][2];
	xlam[2][1] = xlam[1][2];
	xlam[2][2] = 1.+(B.t()-1.)*(B.y()*B.y())*beta2inv;
	xlam[2][3] = (B.t()-1.)*B.y()*B.z()*beta2inv;
	xlam[3][0] = xlam[0][3];
	xlam[3][1] = xlam[1][3];
	xlam[3][2] = xlam[2][3];
	xlam[3][3] = 1.+(B.t()-1.)*(B.z()*B.z())*beta2inv;

	double t_out = vec_in.t()*xlam[0][0] + vec_in.x()*xlam[0][1] + vec_in.y()*xlam[0][2] + vec_in.z()*xlam[0][3];
	double x_out = vec_in.t()*xlam[1][0] + vec_in.x()*xlam[1][1] + vec_in.y()*xlam[1][2] + vec_in.z()*xlam[1][3];
	double y_out = vec_in.t()*xlam[2][0] + vec_in.x()*xlam[2][1] + vec_in.y()*xlam[2][2] + vec_in.z()*xlam[2][3];
	double z_out = vec_in.t()*xlam[3][0] + vec_in.x()*xlam[3][1] + vec_in.y()*xlam[3][2] + vec_in.z()*xlam[3][3];
	FourVector vec_out(x_out,y_out,z_out,t_out);
	
	return vec_out;
  }
  
  //3-vec (4-vec w/3 components) diff^2 function
  static double dif2(FourVector vec1, FourVector vec2){return (vec2.x()-vec1.x())*(vec2.x()-vec1.x()) + (vec2.y()-vec1.y())*(vec2.y()-vec1.y()) + (vec2.z()-vec1.z())*(vec2.z()-vec1.z());}

  //parton class
  //class parton{
  class HHparton : public Parton{
  protected:
	//is shower/thermal originating parton; has been used; is a decayed gluon; is a remnant parton; used for reco; used for stringfrag; is endpoint of string
	bool is_shower_, is_thermal_, is_used_, is_decayedglu_, is_remnant_, used_reco_, used_str_, is_strendpt_, used_junction_;
	//id: particle id, ?orig: particle origin(shower, thermal...)?, par: parent parton (perm. set for gluon decays), status: status of particle (is being considered in loop 'i')
	//string_id: string identifier, pos_str: position of parton in string, endpt_id: denotes which endpoint this parton is (of the string), if it is one
	int alt_id_, orig_, par_, string_id_, pos_str_, endpt_id_, sibling_, PY_par1_, PY_par2_, PY_dau1_, PY_dau2_, PY_stat_, PY_origid_, PY_tag1_, PY_tag2_, PY_tag3_;
	//parton mass
	double alt_mass_;
	
  public:
	//default constructor
	HHparton() : Parton::Parton(1,1,1,0.,0.,0.,0.) {
		is_shower_ = false; is_thermal_ = false; is_used_ = false; is_decayedglu_ = false; is_remnant_ = false; used_reco_ = false; used_str_ = false; is_strendpt_ = false; used_junction_ = false;
		alt_id_ = 0; orig_ = 0; par_ = -1; string_id_ = 0; pos_str_ = 0; endpt_id_ = 0; sibling_ = 0; alt_mass_ = 0.;
		PY_origid_ = 0; PY_par1_ = -1; PY_par2_ = -1; PY_dau1_ = -1; PY_dau2_ = -1; PY_stat_ = 23; /*PY_stat_ = 11;*/ PY_tag1_ = 0; PY_tag2_ = 0; PY_tag3_ = 0;
		set_color(0); set_anti_color(0); set_stat(0);
	}
	
	//getter functions
	double  x() {return x_in().x();} double  y() {return x_in().y();} double  z() {return x_in().z();} double x_t() {return x_in().t();}
	double px() {return p_in().x();} double py() {return p_in().y();} double pz() {return p_in().z();} double   e() {return p_in().t();}
	bool is_shower() {return is_shower_;} bool is_thermal() {return is_thermal_;} bool is_used() {return is_used_;} bool is_decayedglu() {return is_decayedglu_;}
	bool is_remnant() {return is_remnant_;} bool used_reco() {return used_reco_;} bool used_str() {return used_str_;} bool is_strendpt() {return is_strendpt_;}
	bool used_junction() {return used_junction_;}
	int id() {return alt_id_;} int orig() {return orig_;} int par() {return par_;} int string_id() const {return string_id_;} int pos_str() const {return pos_str_;}
	int endpt_id() {return endpt_id_;} int sibling() {return sibling_;} int PY_par1() {return PY_par1_;} int PY_par2() {return PY_par2_;}
	int PY_dau1() {return PY_dau1_;} int PY_dau2() {return PY_dau2_;} int PY_stat() {return PY_stat_;} int PY_origid() {return PY_origid_;}
	int PY_tag1() {return PY_tag1_;} int PY_tag2() {return PY_tag2_;} int PY_tag3() {return PY_tag3_;}
	double mass() {return alt_mass_;}
	FourVector pos() {return x_in();}
	FourVector P()   {return p_in();}
	int status() {return pstat();} int col() {return color();} int acol() {return anti_color();}
	
	//setter functions
	void  x(double val) {x_in_.Set(val,x_in().y(),x_in().z(),x_in().t());}
	void  y(double val) {x_in_.Set(x_in().x(),val,x_in().z(),x_in().t());}
	void  z(double val) {x_in_.Set(x_in().x(),x_in().y(),val,x_in().t());}
	void x_t(double val){x_in_.Set(x_in().x(),x_in().y(),x_in().z(),val);}
	void pos(FourVector val) {x_in_.Set(val.x(),val.y(),val.z(),val.t());}
	void px(double val) {reset_momentum(val,py(),pz(),e());}
	void py(double val) {reset_momentum(px(),val,pz(),e());}
	void pz(double val) {reset_momentum(px(),py(),val,e());}
	void  e(double val) {reset_momentum(px(),py(),pz(),val);}
	void  P(FourVector val) {reset_momentum(val);}
	
	void is_shower(bool val) {is_shower_ = val;} void is_thermal(bool val) {is_thermal_ = val;} void is_used(bool val) {is_used_ = val;} void is_decayedglu(bool val) {is_decayedglu_ = val;}
	void is_remnant(bool val) {is_remnant_ = val;} void used_reco(bool val) {used_reco_ = val;} void used_str(bool val) {used_str_ = val;} void is_strendpt(bool val) {is_strendpt_ = val;}
	void used_junction(bool val) {used_junction_ = val;}
	void id(int val) {alt_id_ = val;} void orig(int val) {orig_ = val;} void par(int val) {par_ = val;}
	void string_id(int val) {string_id_ = val;} void pos_str(int val) {pos_str_ = val;} void endpt_id(int val) {endpt_id_ = val;} void sibling(int val) {sibling_ = val;}
	void PY_par1(int val) {PY_par1_ = val;} void PY_par2(int val) {PY_par2_ = val;} void PY_dau1(int val) {PY_dau1_ = val;} void PY_dau2(int val) {PY_dau2_ = val;}
	void PY_stat(int val) {PY_stat_ = val;} void PY_origid(int val) {PY_origid_ = val;}
	void PY_tag1(int val) {PY_tag1_ = val;} void PY_tag2(int val) {PY_tag2_ = val;} void PY_tag3(int val) {PY_tag3_ = val;}
	void mass(double val) {alt_mass_ = val;}
	void status(int val) {set_stat(val);} void col(int val) {set_color(val);} void acol(int val) {set_anti_color(val);}
	
	//boost functions
	FourVector boost_P(FourVector B){return HHboost(B,p_in());}
	FourVector boost_P(double vx, double vy, double vz){FourVector B(vx,vy,vz,0.); return HHboost(B,p_in());}
	FourVector boost_pos(FourVector B){return HHboost(B,x_in());}
	FourVector boost_pos(double vx, double vy, double vz){FourVector B(vx,vy,vz,0.); return HHboost(B,x_in());}
	
	double pDif2(HHparton comp){return dif2(p_in(),comp.p_in());}
	double posDif2(HHparton comp){return dif2(x_in(),comp.x_in());}
	
  };

	//hadron class
  class HHhadron : public Hadron{
  protected:
	bool is_excited_, is_shsh_, is_shth_, is_thth_, is_recohad_, is_strhad_, is_final_;
	//id: particle id, ?orig: particle origin(sh-sh, sh-th, th-th, recombined, stringfragmented...),
	//par_str: parent string number, parh: parent hadron (decayed), status: status of particle (final, used, decayed, etc...)
	int orig_, parstr_, parh_;
	//hadron mass
	double alt_mass_;
	//vector of partonic parents
	//std::vector<int> parents;
	
  public:
	//default constructor
	HHhadron() : Hadron::Hadron(1,1,1,0.,0.,0.,0.) {
		is_excited_ = false; is_shsh_ = false; is_shth_ = false; is_thth_ = false; is_recohad_ = false; is_strhad_ = false; is_final_ = false;
		orig_ = 0; parstr_ = 0; parh_ = -1; alt_mass_ = 0.; set_stat(0);
	}
	
	//vector of partonic parents
	std::vector<int> parents;
	
	//getter/setter for parents
	int par(int i){if((i>=0) && (i<parents.size())){return parents[i];}else{return 999999;}}
	void add_par(int i){parents.push_back(i);}
	
	//vector of colors (from partons that formed it)
	std::vector<int> cols;
	
	//getter/setter for colors
	int col(int i){if((i>=0) && (i<cols.size())){return cols[i];}else{return -1;}}
	void add_col(int i){cols.push_back(i);}
	
	//getter functions
	double  x() {return x_in().x();} double  y() {return x_in().y();} double  z() {return x_in().z();} double x_t() {return x_in().t();}
	double px() {return p_in().x();} double py() {return p_in().y();} double pz() {return p_in().z();} double   e() {return p_in().t();}
	bool is_excited() {return is_excited_;} bool is_shsh() {return is_shsh_;} bool is_shth() {return is_shth_;} bool is_thth() {return is_thth_;}
	bool is_recohad() {return is_recohad_;} bool is_strhad() {return is_strhad_;} bool is_final() {return is_final_;}
	int orig() {return orig_;} int parstr() {return parstr_;} int parh() {return parh_;}
	double mass() {return alt_mass_;}
	FourVector pos() {return x_in();}
	FourVector P()   {return p_in();}
	int id() {return pid();} int status() {return pstat();}
	
	//setter functions
	void  x(double val) {x_in_.Set(val,x_in().y(),x_in().z(),x_in().t());}
	void  y(double val) {x_in_.Set(x_in().x(),val,x_in().z(),x_in().t());}
	void  z(double val) {x_in_.Set(x_in().x(),x_in().y(),val,x_in().t());}
	void x_t(double val){x_in_.Set(x_in().x(),x_in().y(),x_in().z(),val);}
	void pos(FourVector val) {x_in_.Set(val.x(),val.y(),val.z(),val.t());}
	void px(double val) {reset_momentum(val,py(),pz(),e());}
	void py(double val) {reset_momentum(px(),val,pz(),e());}
	void pz(double val) {reset_momentum(px(),py(),val,e());}
	void  e(double val) {reset_momentum(px(),py(),pz(),val);}
	void  P(FourVector val) {reset_momentum(val);}
	
	void is_excited(bool val) {is_excited_ = val;} void is_shsh(bool val) {is_shsh_ = val;} void is_shth(bool val) {is_shth_ = val;} void is_thth(bool val) {is_thth_ = val;}
	void is_recohad(bool val) {is_recohad_ = val;} void is_strhad(bool val) {is_strhad_ = val;} void is_final(bool val) {is_final_ = val;}
	void orig(int val) {orig_ = val;} void parstr(int val) {parstr_ = val;} void parh(int val) {parh_ = val;}
	void mass(double val) {alt_mass_ = val;}
	void id(int val) {set_id(val);} void status(int val) {set_stat(val);}
	
	//Lorentz boost functions
	FourVector boost_P(FourVector B){return HHboost(B,p_in());}
	FourVector boost_P(double vx, double vy, double vz){FourVector B(vx,vy,vz,0.); return HHboost(B,p_in());}
	FourVector boost_pos(FourVector B){return HHboost(B,x_in());}
	FourVector boost_pos(double vx, double vy, double vz){FourVector B(vx,vy,vz,0.); return HHboost(B,x_in());}
  };

	//class holding a collection of partons
  class parton_collection{
  public:
	//the actual collection of partons
	std::vector<HHparton> partons;
	//default constructor
	parton_collection() {partons.clear();}
	//constructor given a single initial parton
	parton_collection(HHparton par) {partons.push_back(par);}
	//add a parton to the collection
	void add(HHparton par) {partons.push_back(par);}
	//add a collection of partons to the collection
	void add(parton_collection pars) {for (int i=0; i<pars.num(); ++i){partons.push_back(pars[i]);}}
	//get the number of partons in the collection
	int num() {return partons.size();}
	//empty the collection
	void clear() {partons.clear();}
	//insert a parton at position i
	void insert(int i, HHparton newparton) {partons.insert(partons.begin()+i, newparton);}
	//remove the parton at position i
	void remove(int i) {partons.erase(partons.begin()+i);}
	//swap partons at positions i, j
	void swap(int i, int j) {if((i >= 0) && (i < partons.size()) && (j >= 0) && (j < partons.size())) {std::swap(partons[i],partons[j]);}}
	//random access iterator pointing to first element in partons
	std::vector<HHparton>::iterator begin() {return partons.begin();}
	//random access iterator pointing to last element in partons
	std::vector<HHparton>::iterator end() {return partons.end();}
	
	//overloading a few operators; ++, --, []
	//++ adds an empty particle, -- removes last particle, [i] allows access to i'th particle directly
/*	parton_collection& operator++(){HHparton par; partons.push_back(par);}
	parton_collection operator++(int){
		parton_collection temp(*this);
		++(*this);
		return temp;
	}
	parton_collection& operator--(){partons.pop_back();}
	parton_collection operator--(int){
		parton_collection temp(*this);
		--(*this);
		return temp;
	}*/
	HHparton& operator[](int i) {return partons[i];}
	const HHparton& operator[](int i) const {return partons[i];}
  };

	//class holding a collection of hadrons
  class hadron_collection{
  public:
	//the actual collection of hadrons
	std::vector<HHhadron> hadrons;
	//default constructor
	hadron_collection() {hadrons.clear();}
	//constructor given an initial hadron
	hadron_collection(HHhadron had) {hadrons.push_back(had);}
	//add a hadron to the collection
	void add(HHhadron had) {hadrons.push_back(had);}
	//add a collection of partons to the collection
	void add(hadron_collection hads) {for (int i=0; i<hads.num(); ++i){hadrons.push_back(hads[i]);}}
	//get the number of hadrons in the collection
	int num() {return hadrons.size();}
	//empty the collection
	void clear() {hadrons.clear();}
	
	//overloading a few operators; ++, --, []
	//++ adds an empty particle, -- removes last particle, [i] allows access to i'th particle directly
/*	hadron_collection& operator++(){HHhadron had; hadrons.push_back(had);}
	hadron_collection operator++(int){
		hadron_collection temp(*this);
		++(*this);
		return temp;
	}
	hadron_collection& operator--(){hadrons.pop_back();}
	hadron_collection operator--(int){
		hadron_collection temp(*this);
		--(*this);
		return temp;
	}*/
	HHhadron& operator[](int i) {return hadrons[i];}
//	const HHhadron& operator[](int i) {return hadrons[i];}
	const HHhadron& operator[](int i) const {return hadrons[i];}
  };
  
  //used classes
  parton_collection HH_shower, HH_thermal;
  parton_collection HH_showerptns, HH_remnants, HH_pyremn;
  hadron_collection HH_hadrons;

  //function to form strings out of original shower
  void stringform();

  //recombination module
  void recomb();

  //functions to set hadron id based on quark content, mass, and if it's in an excited state
  void set_baryon_id(parton_collection& qrks,HHhadron& had);
  void set_meson_id(parton_collection& qrks,HHhadron& had, int l, int k);

  //gluon to q-qbar splitting function - for recombination use
  void gluon_decay(HHparton& glu, parton_collection& qrks);
  
  //finding a sibling thermal parton for the "ithm'th" thermal parton
  int findthermalsibling(int ithm, parton_collection& therm);
  int findcloserepl(HHparton ptn, int iptn, bool lbt, bool thm, parton_collection& sh_lbt, parton_collection& therm);
  void findcloserepl_glu(HHparton ptn, int iptn, bool lbt, bool thm, parton_collection& sh_lbt, parton_collection& therm, int sel_out[]);
  
  //function to prepare strings for input into Pythia8
  void stringprep(parton_collection& SP_remnants, parton_collection& SP_prepremn, bool cutstr);

  //function to hand partons/strings and hadron resonances (and other color neutral objects) to Pythia8
  bool invoke_py();

  // function to set the spacetime information for the hadrons coming from pythia
  void set_spacetime_for_pythia_hadrons(Pythia8::Event &event, int &size_input, std::vector<int> &eve_to_had, int pythia_attempt);

  // function to 'shake' the event a bit to bring the hadrons to their mass shell
  void bring_hadrons_to_mass_shell(hadron_collection& HH_hadrons);

  void set_initial_parton_masses(parton_collection& HH_showerptns);

  protected:
	static Pythia8::Pythia pythia;
	
};


#endif // HYBRIDHADRONIZATION_H
