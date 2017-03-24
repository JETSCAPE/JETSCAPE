//
//  main01.cpp
//  
//
//  Created by Abhijit Majumder on 11/25/16.
//  A simple program to create a jet using a predetermined
//  algorithm for parton propagation and splitting. 

#include <iostream>
#include <complex>
#include <fstream>
#include <cmath>
#include <cassert>
#include "../JetClass.hpp"
#include "gtest/gtest.h"

TEST(TEST_JET_CLASS, TEST_TRUE){
    double virt;
    double pAssign[4], xLoc[4];
    int i;
    
    double deltaTimeStep;
    
    deltaTimeStep = 0.1;
    
    for (i=0;i<=3; i++) {
        xLoc[i] = 0.0;
    };
    
    pAssign[0] = 11.0;
    
    pAssign[3] = 10.0;
    
    pAssign[1] = 1.0;
    pAssign[2] = 2.0;
    
    FourVector p_in(pAssign);
    
    EXPECT_EQ(p_in.t(), pAssign[0]);
    EXPECT_EQ(p_in.x(), pAssign[1]);
    EXPECT_EQ(p_in.y(), pAssign[2]);
    EXPECT_EQ(p_in.z(), pAssign[3]);
    
    Jet jet(p_in);

    double pmod = std::sqrt(p_in.x() * p_in.x()
                           +p_in.y() * p_in.y()
                           +p_in.z() * p_in.z());

    // test |p| = sqrt(px * px + py * py + pz * pz);
    EXPECT_EQ(pmod, jet.get_jet_p());

    
    // test pseudo-rapidity calculation
    EXPECT_EQ(jet.get_jet_eta(), 0.5*std::log((pmod
              + p_in.z())/(pmod - p_in.z())));
}
