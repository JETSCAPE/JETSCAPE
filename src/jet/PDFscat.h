#ifndef PDFSCAT_H
#define PDFSCAT_H

#include <cmath>
 
#include "Math/DistSampler.h"
#include "Math/Factory.h"
#include "Math/IntegratorMultiDim.h"
#include "TF1.h"

class PDFScat {
	public:
    	PDFScat();
    	~PDFScat();
		void setter(double max_energy0, double low_energy0, double energy_grid0);

		double get_energy(int energy_index);
		double get_rate(int energy_index, int process_index);
		void get_sample(int energy_index, int process_index, double (&V)[4]);

		void initialize_samplers(double b_val);
		double Integrator_LQ(double E, double proc, double b_val);
		double Integrator_HQ(double E, double proc, double b_val, double msq);
		// double Integrator_(double E, double proc);
		// void Sampler_();

    private:
		double obj_low_energy, obj_hig_energy, obj_grid_energy;
		int energy_index_range;
		std::vector<std::vector<double>> rates;
		std::vector<std::vector<ROOT::Math::DistSampler*>> samplers;
		Pythia8::PDFPtr pythiaPDF;
		// std::vector<double>energy_marker;
};

double functionToIntegrate_HQ(double *x, double *params);
double functionToIntegrate_LQ(double *x, double *params);

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

#endif // PDFSCAT_H
