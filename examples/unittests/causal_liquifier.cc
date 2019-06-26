/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
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

#include "CausalLiquefier.h"
#include "gtest/gtest.h"

using namespace Jetscape;

// check coordinate transformation
TEST(CausalLiquifierTest, TEST_COORDINATES){
    
    CausalLiquefier lqf;
    
    // for the transformation from tau-eta to t-z (configuration space)
    EXPECT_DOUBLE_EQ(0.0, lqf.get_t(0.0,0.5));
    EXPECT_DOUBLE_EQ(5.0, lqf.get_t(5.0,0.0));
    EXPECT_DOUBLE_EQ(0.0, lqf.get_z(0.0,0.5));
    EXPECT_DOUBLE_EQ(0.0, lqf.get_z(5.0,0.0));
    
    for(double tau=0.1;tau<0.3;tau+=0.1){
        for(double eta=0.0;eta<0.2;eta+=0.1){
            double t=lqf.get_t(tau,eta);
            double z=lqf.get_z(tau,eta);
            EXPECT_DOUBLE_EQ(tau*tau, t*t-z*z);
            EXPECT_DOUBLE_EQ(exp(2.0*eta), (t+z)/(t-z));
        }
    }
    
    // for the transformation from t-z to tau-eta (momentum)
    EXPECT_DOUBLE_EQ(0.0, lqf.get_ptau(0.0,0.0,0.5));
    EXPECT_DOUBLE_EQ(5.0, lqf.get_ptau(5.0,0.0,0.0));
    EXPECT_DOUBLE_EQ(0.0, lqf.get_peta(0.0,0.0,0.5));
    EXPECT_DOUBLE_EQ(5.0, lqf.get_peta(0.0,5.0,0.0));

}

// check causality
TEST(CausalLiquifierTest, TEST_CAUSALITY){
    
    CausalLiquefier lqf;

    // check zero source outside of the cousal area
    double time = lqf.tau_delay;
    double r_bound = lqf.c_diff * time;
    for(int scale = 1.0; scale < 5.0; scale++ ){
        double r_out = r_bound * scale;
        EXPECT_DOUBLE_EQ(0.0, lqf.causal_diffusion_smooth(time,r_out));
        EXPECT_DOUBLE_EQ(0.0, lqf.causal_diffusion_delta(time,r_out));
        EXPECT_DOUBLE_EQ(0.0, lqf.causal_diffusion_kernel(time,r_out));
    }
    
    // check zero source before the deposition time
    time = 0.99*(lqf.tau_delay-0.5*lqf.dtau);
    std::array<Jetscape::real, 4> jmu = {0.0,0.0,0.0,0.0};
    std::array<Jetscape::real, 4> droplet_xmu = {0.0, 0.0, 0.0, 0.0};
    std::array<Jetscape::real, 4> droplet_pmu = {1.0, 1.0, 0.0, 0.0};
    Droplet drop_i(droplet_xmu, droplet_pmu);
    lqf.smearing_kernel( time, 0.0, 0.0, 0.0, drop_i, jmu);
    EXPECT_DOUBLE_EQ(0.0, jmu[0]);
    
    // check zero source after the deposition time
    time = 1.01*(lqf.tau_delay+0.5*lqf.dtau);
    lqf.smearing_kernel( time, 0.0, 0.0, 0.0, drop_i, jmu);
    EXPECT_DOUBLE_EQ(0.0, jmu[0]);
}



// check conservation law
TEST(CausalLiquifierTest, TEST_CONSEVATION){
    
    CausalLiquefier lqf;

    double dr = 0.005; //in [fm]
    double dt = 0.005; //in [fm]
    
    for(double t=0.15; t < 0.5; t+=dt){
        double integrated_value = 0.0;
        for(double r = 0.5*dr; r<1.5*t; r+=dr){
            integrated_value
            += 4.0*M_PI*r*r*dr*lqf.causal_diffusion_kernel(t,r);
        }
        // require conservation within 5%
        EXPECT_NEAR(1.0,integrated_value,0.05);
    }

}
