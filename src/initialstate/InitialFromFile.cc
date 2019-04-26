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

#include "InitialFromFile.h"

InitialFromFile::InitialFromFile() {
    SetId("InitialFromFile");
    event_id_ = -1;
}


InitialFromFile::~InitialFromFile() {
}

void InitialFromFile::InitTask() {
}

void InitialFromFile::Exec() {  
    Clear();
    Jetscape::JSINFO << "Read initial condition from file";
    try {
        auto * xml_path = GetIniStateXML()->FirstChildElement("initial_profile_path");
        if (!xml_path) {
            throw("Not a valid JetScape IS::initial_profile_path XML section in file!");
        } else {
            event_id_++;
            std::ostringstream path_with_filename;
            path_with_filename << xml_path->GetText() << "/event-" << event_id_
                               << "/initial.hdf5";
            Jetscape::JSINFO << "External initial profile path is"
                           << path_with_filename.str();

            herr_t status;
            std::ostringstream event_group;


            event_group << "/event_0";
            Jetscape::JSINFO << "event_group=" << event_group.str().c_str();
            H5file_ptr_ = H5Fopen(path_with_filename.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
            H5group_ptr_ = H5Gopen(H5file_ptr_, event_group.str().c_str(), H5P_DEFAULT);

            ReadConfigs();
            ReadNbcDist();
            ReadEntropyDist();

            status = H5Gclose(H5group_ptr_);
            status = H5Fclose(H5file_ptr_);
        }
    } catch(std::exception & err) {
        Jetscape::JSWARN << err.what();
        std::exit(-1);
    }
}

void InitialFromFile::ReadConfigs(){
    Jetscape::JSINFO << "Read initial state configurations from file";
    double grid_step = h5_helper_->readH5Attribute_double(H5group_ptr_, "dxy");
    dim_x_ = h5_helper_->readH5Attribute_int(H5group_ptr_, "Nx");
    dim_y_ = h5_helper_->readH5Attribute_int(H5group_ptr_, "Ny");
    double xmax = dim_x_ * grid_step / 2;
    SetRanges(xmax, xmax, 0.0);
    SetSteps(grid_step, grid_step, 0.0);
    Jetscape::JSINFO << "xmax = " << xmax;

    npart = h5_helper_->readH5Attribute_double(H5group_ptr_, "npart");
    ncoll = h5_helper_->readH5Attribute_double(H5group_ptr_, "ncoll");
    totalentropy = h5_helper_->readH5Attribute_double(H5group_ptr_, "mult");
}

void InitialFromFile::ReadNbcDist(){
    Jetscape::JSINFO << "Read number of binary collisions from file";
    auto dataset = H5Dopen(H5group_ptr_, "Ncoll_density", H5P_DEFAULT);
    int dimx = dim_x_;
    int dimy = dim_y_;
    double temp_data[dimx][dimy];
    auto status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                          H5P_DEFAULT, temp_data);
    for (int i = 0; i < dimx; i++) {
        for (int j = 0; j < dimy; j++) {
            num_of_binary_collisions_.push_back(temp_data[i][j]);
        }
    }
    status = H5Dclose(dataset);
}

void InitialFromFile::ReadEntropyDist(){
    Jetscape::JSINFO << "Read initial entropy density distribution from file";
    auto dataset = H5Dopen(H5group_ptr_, "matter_density", H5P_DEFAULT);
    int dimx = dim_x_;
    int dimy = dim_y_;
    double temp_data[dimx][dimy];
    auto status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
                          H5P_DEFAULT, temp_data);
    for (int i = 0; i < dimx; i++) {
        for (int j = 0; j < dimy; j++) {
            entropy_density_distribution_.push_back(temp_data[i][j]);
        }
    }
    status = H5Dclose(dataset);
}

void InitialFromFile::Clear() {
    Jetscape::JSINFO << "clear initial condition vectors";
    entropy_density_distribution_.clear();
    num_of_binary_collisions_.clear();
}


void InitialFromFile::Write(weak_ptr<JetScapeWriter> w){
}
