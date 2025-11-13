#ifndef ELASTICCOLLISION_H
#define ELASTICCOLLISION_H

#include "scat.h"
#include "TVector3.h"
#include <TRotation.h>
#include <random>
#include "FourVector.h"

class ElasticCollision{
    public:
        ElasticCollision();
		void setter(double higher_temp0,double lower_temp0,double div_temp0,double max_energy0,double low_energy0,double energy_grid0);
        int elastic_kinematics(double temperature,int& pid0,int& pid2, int& pid3 ,double (&pc0)[4],double (&pc2)[4],double (&pc3)[4]);
        double V[4];
    	TVector3 iZ_Vector;
        Scat scattering_obj; 
    private:
        double obj_low_temp,obj_hig_temp,obj_grid_temp,obj_low_energy,obj_hig_energy,obj_grid_energy;
        int pid_list[6]={1,2,3,-1,-2,-3};
	    std::default_random_engine generator;
	    //std::uniform_real_distribution<> dis; //between [0,1)              
};

#endif
