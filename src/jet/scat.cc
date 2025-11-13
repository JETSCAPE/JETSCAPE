#include "scat.h"


Scat::Scat(){}

void Scat::setter(double higher_temp0,double lower_temp0,double div_temp0,double max_energy0,double low_energy0,double energy_grid0){
    obj_low_temp=lower_temp0;
    obj_hig_temp=higher_temp0;
    obj_grid_temp=div_temp0;
    obj_hig_energy=max_energy0;
    obj_low_energy=low_energy0;
    obj_grid_energy=energy_grid0;
}

double Scat::get_energy(int energy_index){
        if (energy_index>energy_index_range){
                energy_index=energy_index_range; 
        }
        if (energy_index<=0){
                energy_index=0;
        }
        //std::cout<<"energy index is "<<energy_index<<std::endl;
	return energy_marker[energy_index];
}

double Scat::get_temp(int temp_index){
        if (temp_index>temp_index_range){
                temp_index=temp_index_range; 
        }
        if (temp_index<0){
                temp_index=0;
        }        
	return temp_marker[temp_index];
}

double Scat::get_rate(int temp_index,int energy_index,int process_index){
        if (energy_index>energy_index_range){
                energy_index=energy_index_range; 
        }
        if (energy_index<0){
                energy_index=0;
        }
        if (temp_index>temp_index_range){
                temp_index=temp_index_range; 
        }
        if (temp_index<0){
                temp_index=0;
        }     
        //std::cout<<"energy index is "<<energy_index<<" temp index is "<<temp_index<<" process index is "<<process_index<<std::endl;
	return rates[temp_index][process_index][energy_index];
}
Scat::~Scat() {}

void Scat::get_sample(int temp_index,int energy_index,int process_index,double (&V)[4]){
        samplers[temp_index][process_index][energy_index]->Sample(V);
} 

void Scat::initialize_samplers(double b_val){
        double chm_massSq=1.6129;
        double btm_massSq=17.4724;
        double local_temp, local_energy,temp_rate;
    	double proc_[13]={1,2,3,4,5,6,7,8,9,1,2,1,2};
	energy_index_range = int((obj_hig_energy - obj_low_energy) / obj_grid_energy);
	temp_index_range = int((obj_hig_temp - obj_low_temp) / obj_grid_temp);
    	for (int i2 = 0; i2 <= temp_index_range; i2++) {
		std::vector<std::vector<double>> rates_layer;
		std::vector<std::vector<ROOT::Math::DistSampler*>> samplers_layer;
		local_temp=obj_low_temp+i2*obj_grid_temp;
		temp_marker.push_back(local_temp);
		for (int i1= 0; i1<13; i1++){
			std::vector<double> rates_row;
			std::vector<ROOT::Math::DistSampler*> samplers_row;						
			for (int i0 = 0; i0<=energy_index_range; i0++) {	
				local_energy=obj_low_energy+i0*obj_grid_energy;
				if(i1<9){       //Massless Partons
				        temp_rate=Integrator_LQ(local_energy, local_temp, proc_[i1],b_val);
				        TF1 *f1 = new TF1("myfunc",functionToIntegrate, 0, 1, 4);
				        ROOT::Math::DistSampler *sampler = ROOT::Math::Factory::CreateDistSampler("Foam");
                                        double xmin[4]= {0.0,0.0,0.0,0.0};
                                        double xmax[] = {M_PI, M_PI, 2 * M_PI, local_energy + 1};
        		                double par0[] = {local_energy, local_temp, proc_[i1],b_val};
                		        f1->SetParameters(par0);
				        sampler->SetFunction(*f1, 4);
        		                sampler->SetRange(xmin, xmax);
				        bool ret=sampler->Init();
				        samplers_row.push_back(sampler);
                                        rates_row.push_back(temp_rate);
                                }
                                else if(i1>=9 and i1<11){       //Charm Quark
                                        temp_rate=Integrator_HQ(local_energy, local_temp, proc_[i1],b_val,chm_massSq);
                                        TF1 *f1 = new TF1("myfunc",functionToIntegrateHQ, 0, 1, 5);
                                        ROOT::Math::DistSampler *sampler = ROOT::Math::Factory::CreateDistSampler("Foam");
                                        double xmin[] = {0.0, 0.0, 0.0, 0.0};
                                        double xmax[] = {local_temp*15, M_PI, M_PI, 2.0*M_PI};
                                        double par0[] = {local_energy, local_temp, proc_[i1],b_val,chm_massSq};
        		                f1->SetParameters(par0);
				        sampler->SetFunction(*f1, 4);
        		                sampler->SetRange(xmin, xmax);
				        bool ret=sampler->Init();
				        samplers_row.push_back(sampler);
                                        rates_row.push_back(temp_rate);                                        
                                }
                                else{   //Bottom Quark
                                        temp_rate=Integrator_HQ(local_energy, local_temp, proc_[i1],b_val,btm_massSq);
                                        TF1 *f1 = new TF1("myfunc",functionToIntegrateHQ, 0, 1, 5);
                                        ROOT::Math::DistSampler *sampler = ROOT::Math::Factory::CreateDistSampler("Foam");
                                        double xmin[] = {0.0, 0.0, 0.0, 0.0};
                                        double xmax[] = {local_temp*15, M_PI, M_PI, 2.0*M_PI};
                                        double par0[] = {local_energy, local_temp, proc_[i1],b_val,btm_massSq};
        		                f1->SetParameters(par0);
				        sampler->SetFunction(*f1, 4);
        		                sampler->SetRange(xmin, xmax);
				        bool ret=sampler->Init();
				        samplers_row.push_back(sampler);
                                        rates_row.push_back(temp_rate);                                
                                }


				if(i2==int((obj_hig_temp - obj_low_temp) / obj_grid_temp) && i1==8){
					energy_marker.push_back(local_energy);
				}
			}
			rates_layer.push_back(rates_row);
			samplers_layer.push_back(samplers_row);
		}
		rates.push_back(rates_layer);
		samplers.push_back(samplers_layer);
	}
}

double Scat::Integrator_LQ(double e, double t,double proc,double b_val){
	TF1 *f1=new TF1("myfunc",functionToIntegrate,0,1,4);
	//double par0[]={e,t,proc};
	double xmin[]={0.0, 0.0, 0.0, 0.0};
	double xmax[]={M_PI, M_PI, 2*M_PI, e*1.2};
	double par0[]={e,t,proc,b_val};	

	f1->SetParameters(par0);
	ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kVEGAS,1.E-10,1.E-10,50000);
	ig.SetFunction(*f1,4);
	//double xmin[]={0.0,0.0,0.0,0.0};
	double int_result=ig.Integral(xmin,xmax);
	delete f1;
	//~ROOT::Math::IntegratorMultiDim();
	//~TF1();

	return int_result;
}

double Scat::Integrator_HQ(double e, double t,double proc,double b_val,double mc_sq){
  TF1 *f1=new TF1("myfunc",functionToIntegrateHQ,0,1,5);
  double xmin[]={0.0, 0.0, 0.0, 0.0};
  double xmax[]={15*t, M_PI, M_PI, 2.0*M_PI};
  double par0[]={e,t,proc,b_val,mc_sq}; 

  f1->SetParameters(par0);
  ROOT::Math::IntegratorMultiDim ig(ROOT::Math::IntegrationMultiDim::kVEGAS,1.E-10,1.E-10,50000);
  ig.SetFunction(*f1,4);
  double int_result=ig.Integral(xmin,xmax);
  delete f1;

  return int_result;
}

// Member functions for scattering rates
inline double q1q1b_to_q2q2b(double s, double t, double u) {
    return (4.0/9)*(pow(t,2)+pow(u,2))/pow(s,2);
}

inline double q1bq1_to_q2bq2(double s, double t, double u) {
    return (4.0/9)*(pow(t,2)+pow(u,2))/pow(s,2);
}

inline double q1q2_to_q1q2(double s, double t, double u) {
    return (4.0/9)*(pow(s,2)+pow(u,2))/pow(t,2);
}

inline double q1bq2b_to_q1bq2b(double s, double t, double u) {
    return (4.0/9)*(pow(s,2)+pow(u,2))/pow(t,2);
}

inline double q1q1b_to_q1q1b(double s, double t, double u) {
    return (4.0/9)*((pow(s,2)+pow(u,2))/pow(t,2)+(pow(t,2)+pow(u,2))/pow(s,2)-2*pow(u,2)/(3*s*t));
}

inline double q1bq1_to_q1bq1(double s, double t, double u) {
    return (4.0/9)*((pow(s,2)+pow(u,2))/pow(t,2)+(pow(t,2)+pow(u,2))/pow(s,2)-2*pow(u,2)/(3*s*t));
}

inline double q1q1_to_q1q1(double s, double t, double u) {
    return (4.0/9)*((pow(u,2)+pow(s,2))/pow(t,2)+(pow(t,2)+pow(s,2))/pow(u,2)-2*pow(s,2)/(3*u*t));
}

inline double q1bq1b_to_q1bq1b(double s, double t, double u) {
    return (4.0/9)*((pow(u,2)+pow(s,2))/pow(t,2)+(pow(t,2)+pow(s,2))/pow(u,2)-2*pow(s,2)/(3*u*t));
}

inline double q1q1b_to_gg(double s, double t, double u) {
    return (32.0/27)*(u/t + t/u -(9.0/4)*(pow(t,2)+pow(u,2))/pow(s,2));
}

inline double q1bq1_to_gg(double s, double t, double u) {
    return (32.0/27)*(u/t + t/u -(9.0/4)*(pow(t,2)+pow(u,2))/pow(s,2));
}

inline double q1g_to_q1g(double s, double t, double u) {
    return (4.0/9)*(((-1*u)/s)+((-1*s)/u)+(9.0/4)*(pow(s,2)+pow(u,2))/pow(t,2));
}

inline double q1bg_to_q1bg(double s, double t, double u) {
    return (4.0/9)*(((-1*u)/s)+((-1*s)/u)+(9.0/4)*(pow(s,2)+pow(u,2))/pow(t,2));
}

inline double gq1_to_gq1(double s, double t, double u) {
    return (4.0/9)*(((-1*u)/s)+((-1*s)/u)+(9.0/4)*(pow(s,2)+pow(u,2))/pow(t,2));
}

inline double gg_to_q1q1b(double s, double t, double u) {
    return (1/6.0)*(u/t+t/u-(9.0/4)*(pow(t,2)+pow(u,2))/pow(s,2));
}

inline double gg_to_gg(double s, double t, double u) {
    return (9.0/2)*(3-t*u/pow(s,2)-s*u/pow(t,2)-s*t/pow(u,2));
}

inline double cq_to_cq(double s,double t,double u,double mc_sq){
        return ((4.0/9.0)*(pow((mc_sq-u),2)+pow(s-mc_sq,2)+2*mc_sq*t)/pow(t,2));
}

inline double cg_to_cg(double s,double t,double u,double mc_sq){     //problem w matrix ele
        return ((2.0*(s-mc_sq)*(mc_sq-u))/pow(t,2)+
                (4.0/9.0)*(((s-mc_sq)*(mc_sq-u)+2.0*mc_sq*(s+mc_sq))/(pow(s-mc_sq,2))+((s-mc_sq)*(mc_sq-u)+2.0*mc_sq*(u+mc_sq))/(pow(mc_sq-u,2))+
        (mc_sq*(4.0*mc_sq-t))/(4.0*(s-mc_sq)*(mc_sq-u)))
                +(1.0)*(((s-mc_sq)*(mc_sq-u)+mc_sq*(s-u))/(t*(s-mc_sq))-
        ((s-mc_sq)*(mc_sq-u)-mc_sq*(s-u))/(t*(mc_sq-u))));
}

long double ThermalDistribution(double x, double Temp, double stat,double b_val)
{
        double a,b;
        a=1;
        b=b_val;//40(works);
        long double EXP;
        if (b_val>1000){
            EXP = exp(x/(a*Temp));
        }
        else{   
             EXP = exp(x/(a*Temp)+x/(b*pow(Temp,2)));    //exp(x/Temp); exp(-x/(a*Temp)) //exp(x/(a*Temp)-x/(b*pow(Temp,2)))
        }
        EXP = exp(x/(a*Temp));
        //std::cout<<"b_val is "<<b_val<<"\n";
        if (stat == 0)
        {
                return 1.0/(1+EXP); //Fermi Dirac
        }
        else
        {
                return 1.0/(EXP-1); //Bose
        }
}

double functionToIntegrate(double *x, double *params) {
        double E1 = params[0];
        double Temp = params[1];
        double Process = params[2];
        double b_val = params[3];

        int g_b;//generacy factor
        double c23 = cos(x[0])*cos(x[1])+sin(x[0])*sin(x[1])*cos(x[2]);         //cos(theta_23)

        double E2 = E1*x[3]*(1-cos(x[1]))/(E1*(1-cos(x[0]))-x[3]*(1-c23)); // Thermal parton energy

        double s = 2*E1*E2*(1-cos(x[0]));
        double t = -2*E1*x[3]*(1-cos(x[1]));
        double u = -s-t;

        double EScale2 = std::max(std::pow(std::abs(t),2.0),std::pow(M_PI*Temp,2.0));
        double Nf = 3.0;
        double QCDScale2 = std::pow(0.2,2.0); //MS bas scheme 0.33; 0.2 LBT value
        double alphaS = 0.3;//4.0*M_PI/(11.0-2.0*Nf/3.0)*1/std::log(EScale2/QCDScale2);

        double debye = alphaS*pow(Temp,2.0)*(4.0*M_PI)*1.5;
        if(s <= 2*debye || t <= -s+debye || t >= -debye){
                return 0;
        }

        double matrix_element;
        double Stat;
        switch ((int)Process)
        {
                case 1:
                        matrix_element = 0.5*q1q1b_to_q2q2b(s,t,u); // 0.5*q1q1b_to_q2q2b(s,t,u);     //problem
                        Stat = 0;
                        g_b = 6;
                        break;
        
                case 2:
                        matrix_element = 0.5*q1q1b_to_q1q1b(s,t,u); //0.5*q1q1b_to_q1q1b(s,t,u);     //problem
                        Stat = 0;
                        g_b = 6;
                        break;

                case 3:
                        matrix_element = 0.5*q1q1_to_q1q1(s,t,u);       //correct
                        Stat = 0;
                        g_b = 6;
                        break;
                case 4:
                        matrix_element = 0.5*q1q1b_to_gg(s,t,u); //correct
                        Stat = 0;
                        g_b = 6;
                        break;
                case 5:
                        matrix_element = q1g_to_q1g(s,t,u);     //correct
                        Stat = 1;
                        g_b = 16;
                        break;
                case 6:
                        matrix_element = 4*q1q2_to_q1q2(s,t,u); //4 q and qb contribution of both considered    //correct
                        Stat = 0;
                        g_b = 6;
                        break;
                case 7:
                        matrix_element = gg_to_q1q1b(s,t,u);            //correct
                        Stat = 1;
                        g_b = 16;
                        break;
                case 8:
                        matrix_element = 0.5*gg_to_gg(s,t,u);           //correct
                        Stat = 1;
                        g_b = 16;
                        break;
                case 9:
                        matrix_element = 6*gq1_to_gq1(s,t,u);        //both q and qb contribution //correct
                        Stat=0;
                        g_b=6;
                        break;
       }
        double thermal_dist = 0;
        thermal_dist = ThermalDistribution(E2,Temp,Stat,b_val);

        return std::pow(alphaS,2.0)*pow(4*M_PI,2)*g_b/(pow(2*M_PI,4)*16.0*E1)*thermal_dist*matrix_element*pow(E2,2)*2*x[3]*sin(x[0])*sin(x[1])/abs(t);

}

double functionToIntegrateHQ(double x[], double* params){
        double E1 = params[0];
        double Temp = params[1];
        double Process = params[2];    
        double b_val = params[3];
        double mc_sq = params[4];
        int g_b;

        //x[0]=E2; x[1]=theta2; x[2]=theta4; x[3]=phi4
        //on-shell
        double c24=sin(x[1])*sin(x[2])*cos(x[3])+cos(x[1])*cos(x[2]);
        double p1=sqrt(E1*E1-mc_sq);
        double E4=(E1*x[0]- p1*x[0]*cos(x[1]))/(E1-p1*cos(x[2])+x[0]-x[0]*c24);
        
        double s = mc_sq+ 2*(E1*x[0]-p1*x[0]*cos(x[1]));
        double u = mc_sq-2*(E1*E4-p1*E4*cos(x[2]));
        double t = 2*mc_sq-s-u;

        double debye = 0.3*pow(Temp,2)*(4*M_PI)*(4.5/3.0);
        //std::cout<<"Values fom scattering.cpp s "<<s<<" t "<<t<<" u "<<u<<std::endl;
        if(s <= 2*debye || t <= -s+debye || t >= -debye){
                //std::cout<<"RETURN 0!";
                return 0;
        }
        double matrix_element;
        double Stat;
        switch ((int)Process)
        {
                case 1:
                        matrix_element = 6*cq_to_cq(s,t,u,mc_sq);
                        Stat = 0;
                        g_b = 6;
                        break;  
                case 2:
                        matrix_element = cg_to_cg(s,t,u,mc_sq);
                        Stat = 1;
                        g_b = 16;
                        break;  
        }
        long double thermal_dist = ThermalDistribution(x[0],Temp,Stat,b_val);
        long double extra_stat = ThermalDistribution(E4,Temp,Stat,b_val);
        if (Stat==0){extra_stat=extra_stat*(-1.0);}
        //std::cout<<"HURRAY";
        double answer=pow(0.3,2)*pow(4*M_PI,2)*g_b/(pow(2*M_PI,4)*16.0*E1)*thermal_dist*matrix_element*(1+extra_stat)*(x[0]*E4*sin(x[1])*sin(x[2]))/(E1-p1*cos(x[2])+x[0]-x[0]*c24);
        if (answer<0){std::cout<<"ALERT!";}
        return answer;         
}
