#include "PDFElasticCollision.h"

void Rectify1(TVector3 &v6){
	double threshold = 1e-10;
	if (std::abs(v6.X()) < threshold){ v6.SetX(0); }
	if (std::abs(v6.Y()) < threshold){ v6.SetY(0); } 
	if (std::abs(v6.Z()) < threshold){ v6.SetZ(0); }	
}

void Rectify_Momentum(double (&pc0)[4], const double E1, const double msq){
	double currentMagnitude = std::sqrt(pc0[1] * pc0[1] + pc0[2] * pc0[2] + pc0[3] * pc0[3] + msq);
	pc0[1] = E1*pc0[1]/currentMagnitude;
	pc0[2] = E1*pc0[2]/currentMagnitude;
	pc0[3] = E1*pc0[3]/currentMagnitude;
}

std::uniform_real_distribution<> uniform_rand(0.0,1.0);
PDFElasticCollision::PDFElasticCollision() {};

void PDFElasticCollision::setter(double max_energy0, double low_energy0, double energy_grid0) {
    obj_hig_energy = max_energy0;
    obj_low_energy = low_energy0;
    obj_grid_energy = energy_grid0;
    iZ_Vector.SetXYZ(0,0,1);
    scattering_obj.setter(max_energy0, low_energy0, energy_grid0);
	//scattering_obj.setter(0.2,0.2,0.1,102,0.5,0.5);
	//scattering_obj.setter(0.2,0.2,0.1,102,0.5,1);
    scattering_obj.initialize_samplers(9999); //the param here is b_val
}
/*The MATTER module had pre-defination  of calculating the probability of scattring.*/

int PDFElasticCollision::elastic_kinematics(int &pid0, int &pid2, int &pid3, double (&pc0)[4], double (&pc2)[4], double (&pc3)[4]) {
	double r0,r1,r2,r3,r4,r5;
    double E1,E2,E3,E4,c23,c24,p1,p3,THETA2,THETA3,THETA4,PHI3,PHI4,PHI2,PHI23;
    double mc_sq;

    int parent_pid,hole_pid,daughter1_pid,daughter2_pid,pid_index;
    int energy_index;//, temp_index;
    int parton_type = -1;
    int proc; //heavy(1) or light(0)

    TVector3 P1;
    TVector3 P2;
    TVector3 P3;
    TVector3 P4;
    TRotation r;

    E1 = pc0[0];

    //find which E1 bin to use and switch over to that discretized value
    energy_index = round((E1- obj_low_energy)/obj_grid_energy);
    E1 = scattering_obj.get_energy(energy_index);

    parent_pid = pid0;
    pid_index = floor(uniform_rand(generator)*6); //ONLY UDS FOR NOW
    if (parent_pid == 21) {
		r0 = scattering_obj.get_rate(energy_index,6); //has color abiguity for matter
		r1 = scattering_obj.get_rate(energy_index,7);
		r2 = scattering_obj.get_rate(energy_index,8);
        if (exp(-(r0+r1+r2))<=((double)rand())/RAND_MAX) {
            parton_type = 0;
    		std::discrete_distribution<int> distribution0{r0,r1,r2};
            switch (distribution0(generator)) {
                case 0: //g g -> q qbar
                    scattering_obj.get_sample(energy_index,6,V);
                    hole_pid      =  21;
                    daughter1_pid =  pid_list[pid_index];
                    daughter2_pid = -pid_list[pid_index];
                    break;

                case 1: //g g -> g g
                    scattering_obj.get_sample(energy_index,7,V);
                    hole_pid      = 21;
                    daughter1_pid = 21;
                    daughter2_pid = 21;
                    break;

                case 2: //g q -> g q
                    scattering_obj.get_sample(energy_index,8,V);
                    hole_pid      = pid_list[pid_index];
                    daughter1_pid = 21;
                    daughter2_pid = pid_list[pid_index];

                default:
                    //never gets here
                    //raise error?
                    break;
            }
        }   
    }
    else if (abs(parent_pid)>=1 && abs(parent_pid)<=3) {
		r0 = scattering_obj.get_rate(energy_index,0);
		r1 = scattering_obj.get_rate(energy_index,1);
		r2 = scattering_obj.get_rate(energy_index,2);
		r3 = scattering_obj.get_rate(energy_index,3);
		r4 = scattering_obj.get_rate(energy_index,4);
		r5 = scattering_obj.get_rate(energy_index,5);
        if (exp(-(r0+r1+r2+r3+r4+r5))<=((double)rand())/RAND_MAX) {
            parton_type = 0;
    		std::discrete_distribution<int> distribution0{r0,r1,r2,r3,r4,r5};
            switch (distribution0(generator)) {
                case 0: //q1 q1bar -> q2 q2bar
                    scattering_obj.get_sample(energy_index,0,V);
                    hole_pid      = -parent_pid;
                    do { //keep sampling until q2 != q1
                        pid_index = floor(uniform_rand(generator)*6);
                        daughter1_pid = pid_list[pid_index];
                    } while (daughter1_pid == parent_pid);                
                    daughter2_pid = -daughter1_pid;
                    break;

                case 1: //q1 q1bar -> q1 q1bar
                    scattering_obj.get_sample(energy_index,1,V);
                    hole_pid      = -parent_pid;
                    daughter1_pid =  parent_pid;
                    daughter2_pid = -parent_pid;
                    break;

                case 2: //q1 q1 -> q1 q1
                    scattering_obj.get_sample(energy_index,2,V);
                    hole_pid      = parent_pid;
                    daughter1_pid = parent_pid;
                    daughter2_pid = parent_pid;
                    break;

                case 3: //q1 q1bar -> g g
                    scattering_obj.get_sample(energy_index,3,V);
                    hole_pid      = -parent_pid;
                    daughter1_pid = 21;
                    daughter2_pid = 21;
                    break;

                case 4: //q1 g -> q1 g
                    scattering_obj.get_sample(energy_index,4,V);
                    hole_pid      = 21;
                    daughter1_pid = parent_pid;
                    daughter2_pid = 21;
                    break;

                case 5: //q1 q2 -> q1 q2
                    scattering_obj.get_sample(energy_index,5,V);
                    daughter1_pid = parent_pid;
                    do { //keep sampling until q2 != q1 or q1bar
                        pid_index = floor(uniform_rand(generator)*6);
                        daughter2_pid = pid_list[pid_index];
                    } while (daughter2_pid == parent_pid || daughter2_pid == -parent_pid);
                    hole_pid=daughter2_pid;     
                    break;

                default:
                    //never gets here
                    break;
            }
        } 
    }
     /*Heavy scattering part will be added here*/
    else if (abs(parent_pid)==4) {
        mc_sq = 1.6129;

		r0 = scattering_obj.get_rate(temp_index,energy_index,9);
		r1 = scattering_obj.get_rate(temp_index,energy_index,10);
        if (exp(-(r0+r1))<=((double)rand())/RAND_MAX) {
            parton_type = 1;
    		std::discrete_distribution<int> distribution0{r0,r1};
            switch (distribution0(generator)) {
                case 0: //TODO: should be q1 q1bar -> q2 q2bar?? 
                    scattering_obj.get_sample(temp_index,energy_index,9,V);
                    hole_pid      = pid_list[pid_index];
                    daughter1_pid = parent_pid;
                    daughter2_pid = pid_list[pid_index];
                    break;

                case 1: //TODO: should be q1 q1bar -> q1 q1bar?
                    scattering_obj.get_sample(temp_index,energy_index,10,V);
                    hole_pid      = 21;
                    daughter1_pid = parent_pid;
                    daughter2_pid = 21;
                    break;

                default:
                    //never gets here
                    break;
            }
        }
    }
    else {
        mc_sq = 17.4724;

		r0 = scattering_obj.get_rate(temp_index,energy_index,11);
		r1 = scattering_obj.get_rate(temp_index,energy_index,12);
        if (exp(-(r0+r1))<=((double)rand())/RAND_MAX) {
            parton_type = 1;
    		std::discrete_distribution<int> distribution0{r0,r1};
            switch (distribution0(generator)) {
                case 0: //TODO: should be q1 q1bar -> q2 q2bar?? 
                    scattering_obj.get_sample(temp_index,energy_index,11,V);
                    hole_pid      = pid_list[pid_index];
                    daughter1_pid = parent_pid;
                    daughter2_pid = pid_list[pid_index];
                    break;

                case 1: //TODO: should be q1 q1bar -> q1 q1bar?
                    scattering_obj.get_sample(temp_index,energy_index,12,V);
                    hole_pid      = 21;
                    daughter1_pid = parent_pid;
                    daughter2_pid = 21;
                    break;

                default:
                    //never gets here
                    break;
            }
        }
    }

    if (parton_type == 0) { //light parton
        THETA2 = V[0];
        THETA3 = V[1];
        PHI23  = V[2];
        E3     = V[3];

		c23 = cos(THETA2)*cos(THETA3)+sin(THETA2)*sin(THETA3)*cos(PHI23);
		
        E2 = (E1*E3*(1-cos(THETA3)))/(E1*(1-cos(THETA2))-E3*(1-c23));
		E4 = E1+E2-E3;
		
        PHI2 = 2*M_PI*uniform_rand(generator);
		PHI3 = PHI2-PHI23;

		THETA4= acos((E1+E2*cos(THETA2)-E3*cos(THETA3))/E4);

		double value = (E2*sin(THETA2)*cos(PHI2)-E3*sin(THETA3)*cos(PHI3))/(E4*sin(THETA4));
		value = std::max(-1.0, std::min(1.0, value)); // Clamping the value
		PHI4 = acos(value);

        //std::cout<<"LQ E1 "<<E1<<" E2 "<<E2<<" E3 "<<E3<<" E4 "<<E4<<std::endl;
        P1.SetXYZ(pc0[1],pc0[2],pc0[3]);
        P2.SetXYZ(E2*sin(THETA2)*cos(PHI2),E2*sin(THETA2)*sin(PHI2),E2*cos(THETA2));
		P3.SetXYZ(E3*sin(THETA3)*cos(PHI3),E3*sin(THETA3)*sin(PHI3),E3*cos(THETA3));
		P4.SetXYZ(E4*sin(THETA4)*cos(PHI4),E4*sin(THETA4)*sin(PHI4),E4*cos(THETA4));
		/*
		double s0=P1.Angle(iZ_Vector);
		TVector3 s1;
		s1=P1.Cross(iZ_Vector);
        
		r.Rotate(s0,s1);
        P2=r*P2;
        P3=r*P3;
        P4=r*P4;
		*/
        //std::cout<<"LQ  after one roatation E1 "<<E1<<" E2 "<<E2<<" E3 "<<E3<<" E4 "<<E4<<std::endl;
    }
    else if (parton_type==1) { //heavy
		E2     = V[0];
        THETA2 = V[1];
        THETA4 = V[2];
        PHI4   = V[3];

		//std::cout<<"E2 "<<E2<<" theta2 "<<THETA2<<" theta4 "<<THETA4<<" phi4 before "<<PHI4<<std::endl;
        c24 = sin(THETA2)*sin(THETA4)*cos(PHI4)+cos(THETA2)*cos(THETA4);

		if (E1*E1<mc_sq){return -1;}
		Rectify_Momentum(pc0,E1,mc_sq);

		//std::cout<<"E1 "<<E1*E1<<" mcsq "<<mc_sq<<std::endl;
        p1 = sqrt(E1*E1-mc_sq);
        E4 = (E1*E2-p1*E2*cos(THETA2))/(E1-p1*cos(THETA4)+E2-E2*c24);
        E3 = E1+E2-E4;
		if (E3*E3<mc_sq){return -1;}
        p3 = sqrt(E3*E3-mc_sq);

		//std::cout<<"p1 "<<p1<<" E4 "<<E4<<" E3 "<<E3<<" p3 after "<<p3<<std::endl;

		double value = (p1+E2*cos(THETA2)-E4*cos(THETA4))/p3;
		if (value>1 || value<-1) {std::cout<<"value cos "<<value<<std::endl;}
		value = std::max(-1.0, std::min(1.0, value));
        THETA3 = acos(value);

		value = -(E4*sin(PHI4)*sin(THETA4))/(p3*sin(THETA3));
		if (value>1 || value<-1){std::cout<<"value sin "<<value<<std::endl;}
		value = std::max(-1.0, std::min(1.0, value));
        PHI3 = asin(value);

        PHI2 = 2*M_PI*uniform_rand(generator);     
        r.RotateZ(PHI2);
        //std::cout<<"HQ E1 "<<E1<<" E2 "<<E2<<" E3 "<<E3<<" E4 "<<E4<<std::endl;
        P1.SetXYZ(pc0[1],pc0[2],pc0[3]);
        P2.SetXYZ(E2*sin(THETA2),0,E2*cos(THETA2));
        P3.SetXYZ(p3*sin(THETA3)*cos(PHI3),p3*sin(THETA3)*sin(PHI3),p3*cos(THETA3));
        P4.SetXYZ(E4*sin(THETA4)*cos(PHI4),E4*sin(THETA4)*sin(PHI4),E4*cos(THETA4));

        P2=r*P2;
        P3=r*P3;
        P4=r*P4;
 		Rectify1(P2);
 		Rectify1(P3);
 		Rectify1(P4);	
        //std::cout<<"HQ after one rotation E1 "<<E1<<" E2 "<<E2<<" E3 "<<E3<<" E4 "<<E4<<std::endl;
    }
    else { return -1; }

    //TODO: LEARN ROTATIONS ETC
    double s0 = P1.Angle(iZ_Vector);
    TVector3 s1;
    //s1=P1.Cross(iZ_Vector);
    s1=iZ_Vector.Cross(P1);
    TRotation w0;
    w0.Rotate(s0,s1);
    
    P2=w0*P2;
    P3=w0*P3;
    P4=w0*P4;
 	Rectify1(P2);
 	Rectify1(P3);
 	Rectify1(P4);

    pc3[0]=E2;
    pc3[1]=P2.x();
    pc3[2]=P2.y();
    pc3[3]=P2.z();
    pid3=hole_pid;

    pc0[0]=E3;
    pc0[1]=P3.x();
    pc0[2]=P3.y();
    pc0[3]=P3.z();
    pid0=daughter1_pid;

    pc2[0]=E4;
    pc2[1]=P4.x();
    pc2[2]=P4.y();
    pc2[3]=P4.z();   
    pid2=daughter2_pid; 
    return parton_type;
}