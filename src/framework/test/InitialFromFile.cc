// Copyright @ LongGang Pang
#include "InitialFromFile.h"

InitialFromFile::InitialFromFile() {
    SetId("InitialFromFile");
}


InitialFromFile::~InitialFromFile() {
}

void InitialFromFile::InitTask() {
}

void InitialFromFile::Exec() {  
    Jetscape::INFO << "Read initial condition from file";
    try {
        auto * xml_path = GetIniStateXML()->FirstChildElement("initial_profile_path");
        if (!xml_path) {
            throw("Not a valid JetScape IS::initial_profile_path XML section in file!");
        } else {
            std::string path = xml_path->GetText();
            Jetscape::INFO << "External initial profile path is" << path;
            int event_id = 0;

            herr_t status;
            std::ostringstream event_group;
            event_group << "/event_" << event_id;
            Jetscape::INFO << "event_group=" << event_group.str().c_str();
            H5file_ptr_ = H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
            H5group_ptr_ = H5Gopen(H5file_ptr_, event_group.str().c_str(), H5P_DEFAULT);

            read_nbc_dist_();
            read_entropy_dist_();
            read_configs_();

            status = H5Gclose(H5group_ptr_);
            status = H5Fclose(H5file_ptr_);
        }
    } catch(std::exception & err) {
        Jetscape::WARN << err.what();
        std::exit(-1);
    }
}

void InitialFromFile::read_configs_(){
    Jetscape::INFO << "Read initial state configurations from file";
    double grid_step = h5_helper_->readH5Attribute_double(H5group_ptr_, "dxy");
    int total_grids = entropy_density_distribution_.size();
    int nx = int(std::sqrt(total_grids));
    double xmax = nx * grid_step / 2;
    set_ranges(xmax, xmax, 0.0);
    set_steps(grid_step, grid_step, 0.0);
    Jetscape::INFO << "xmax = " << xmax;
}

void InitialFromFile::read_nbc_dist_(){
    Jetscape::INFO << "Read number of binary collisions from file";
    auto dataset = H5Dopen(H5group_ptr_, "Ncoll_density", H5P_DEFAULT);
    int dimx = 100;
    int dimy = 100;
    double temp_data[dimx][dimy];
    auto status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_data);
    for(int i=0; i<dimx; i++)
        for(int j=0; j<dimy; j++) {
            num_of_binary_collisions_.push_back(temp_data[i][j]);
        }
    status = H5Dclose(dataset);
}

void InitialFromFile::read_entropy_dist_(){
    Jetscape::INFO << "Read initial entropy density distribution from file";
    auto dataset = H5Dopen(H5group_ptr_, "matter_density", H5P_DEFAULT);
    int dimx = 100;
    int dimy = 100;
    double temp_data[dimx][dimy];
    auto status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp_data);
    for(int i=0; i<dimx; i++)
        for(int j=0; j<dimy; j++) {
            entropy_density_distribution_.push_back(temp_data[i][j]);
        }
    status = H5Dclose(dataset);
}

void InitialFromFile::Clear() {
    Jetscape::INFO << "clear initial condition vectors";
    entropy_density_distribution_.clear();
    num_of_binary_collisions_.clear();
}


void InitialFromFile::Write(weak_ptr<JetScapeWriter> w){
}
