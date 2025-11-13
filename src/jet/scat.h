#ifndef SCAT_H
#define SCAT_H


#include <cmath>
 
#include "Math/DistSampler.h"
#include "Math/Factory.h"
#include "Math/IntegratorMultiDim.h"
#include "TF1.h"



class Scat {

	public:
    	Scat();
		void setter(double higher_temp0,double lower_temp0,double div_temp0,double max_energy0,double low_energy0,double energy_grid0);
		double get_energy(int energy_index);
		double get_temp(int temp_index);
		double get_rate(int temp_index,int energy_index,int process_index);
		void get_sample(int temp_index,int energy_index,int process_index,double (&V)[4]);
		void initialize_samplers(double b_val);		
		double Integrator_LQ(double e, double t,double proc,double b_val);
		double Integrator_HQ(double e, double t,double proc,double b_val,double mc_sq);
		double Integrator_(double e, double t,double proc);
		void Sampler_();
    ~Scat();

    private:
		double obj_low_temp,obj_hig_temp,obj_grid_temp,obj_low_energy,obj_hig_energy,obj_grid_energy;
		int energy_index_range,temp_index_range;
		std::vector<std::vector<std::vector<double>>> rates;
		std::vector<std::vector<std::vector<ROOT::Math::DistSampler*>>> samplers;
		std::vector<double>energy_marker;
		std::vector<double>temp_marker;		
};

double functionToIntegrateHQ(double x[], double* params);
double functionToIntegrate(double *x, double *params);
long double ThermalDistribution(double x, double Temp, double stat,double b_val);
double q1q1b_to_q2q2b(double s, double t, double u);
double q1bq1_to_q2bq2(double s, double t, double u);
double q1q2_to_q1q2(double s, double t, double u);
double q1bq2b_to_q1bq2b(double s, double t, double u);
double q1q1b_to_q1q1b(double s, double t, double u);
double q1bq1_to_q1bq1(double s, double t, double u);
double q1q1_to_q1q1(double s, double t, double u);
double q1bq1b_to_q1bq1b(double s, double t, double u);
double q1q1b_to_gg(double s, double t, double u);
double q1bq1_to_gg(double s, double t, double u);
double q1g_to_q1g(double s, double t, double u);
double q1bg_to_q1bg(double s, double t, double u);
double gq1_to_gq1(double s, double t, double u);
double gg_to_q1q1b(double s, double t, double u);
double gg_to_gg(double s, double t, double u);

#endif // SCATTERINGRATE_H
