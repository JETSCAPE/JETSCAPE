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

#include <string>
#include <fstream>
#include "IPGlasmaWrapper.h"

// Register the module with the base class
RegisterJetScapeModule<IPGlasmaWrapper> IPGlasmaWrapper::reg("IPGlasma");


IPGlasmaWrapper::IPGlasmaWrapper() {
    SetId("IPGlasma");
    event_id_ = 0;
}


IPGlasmaWrapper::~IPGlasmaWrapper() {}


void IPGlasmaWrapper::InitTask() {
    IPGlasma_ptr_ = std::unique_ptr<IPGlasma>(
                                new IPGlasma(0, 1, 1, "ipglasma.input"));
}


void IPGlasmaWrapper::Exec() {
    Clear();
    Jetscape::JSINFO << "Run IPGlasma ...";
    try {
        IPGlasma_ptr_->generateAnEvent(event_id_);
        event_id_++;
    } catch (std::exception &err) {
        Jetscape::JSWARN << err.what();
        std::exit(-1);
    }
    ReadNbcList("NcollList0.dat");
}


void IPGlasmaWrapper::Clear() {
    Jetscape::JSINFO << "clear initial condition vectors";
}


void IPGlasmaWrapper::ReadNbcList(std::string filename) {
    Jetscape::JSINFO << "Read in binary collision list from "
                     << filename << "...";
    std::ifstream infile(filename.c_str());
    if (!infile.good()) {
        Jetscape::JSWARN << "Can not open " << filename;
        exit(1);
    }

    double x, y;
    infile >> x >> y;
    while (!infile.eof()) {
        binary_collision_x_.push_back(x);
        binary_collision_y_.push_back(y);
        infile >> x >> y;
    }
    infile.close();
    ncoll_ = binary_collision_x_.size();
    rand_int_ptr_ = (
        std::make_shared<std::uniform_int_distribution<int>>(0, ncoll_-1));
    Jetscape::JSINFO << "Ncoll = " << ncoll_;
}


void IPGlasmaWrapper::SampleABinaryCollisionPoint(double &x, double &y) {
    int rand_idx = (*rand_int_ptr_)(*GetMt19937Generator());
    x = binary_collision_x_[rand_idx];
    y = binary_collision_y_[rand_idx];
}


void IPGlasmaWrapper::Write(weak_ptr<JetScapeWriter> w) {}
