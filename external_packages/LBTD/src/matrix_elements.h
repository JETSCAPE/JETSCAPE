#ifndef MATRIX_ELEMENTS_H
#define MATRIX_ELEMENTS_H
#include <cstdlib>


//======running coupling=======================================================
double alpha_s(double Q2, double T);
double f_LPM(double x);
//======Debye mass class=======================================================
class Debye_mass{
private:
	const double TL, TH;
	const size_t NT;
	const double dT;
	const unsigned int type;
	double * mD2;
public:
	Debye_mass(const unsigned int _type);
	~Debye_mass(){delete[] mD2;};
	double get_mD2(double T);
};

//=====For external initialization of debye mass==============================
void initialize_mD_and_scale(const unsigned int type, double scale);

//========resonance?? Q+q-->D-->Q+q============
double dX_res_dt(const double t, void * params);

//=============Baisc function for Q+q --> Q+q==================================
double M2_Qq2Qq(const double t, void * params);
double M2_Qq2Qq_rad(const double t, void * params);
double dX_Qq2Qq_dt(const double t, void * params);

//=============Baisc function for Q+g --> Q+g==================================
double M2_Qg2Qg(const double t, void * params);
double M2_Qg2Qg_rad(const double t, void * params);
double dX_Qg2Qg_dt(const double t, void * params);

//=============Baisc function for Q+q --> Q+q+g==================================
double M2_Qq2Qqg(const double * x_, void * params_);
//=============Baisc function for Q+g --> Q+g+g==================================
double M2_Qg2Qgg(const double * x_, void * params_);

//=============Baisc function for Q+q+g --> Q+q==================================
double Ker_Qqg2Qq(const double * x_,  void * params_);
//=============Baisc function for Q+g+g --> Q+g==================================
double Ker_Qgg2Qg(const double * x_, void * params_);


#endif
