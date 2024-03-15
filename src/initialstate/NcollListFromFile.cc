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
#include "NcollListFromFile.h"

// Register the module with the base class
RegisterJetScapeModule<NcollListFromFile> NcollListFromFile::reg(
        "NcollListFromFile");

NcollListFromFile::NcollListFromFile() {
  SetId("NcollListFromFile");
  event_id_ = -1;
}

NcollListFromFile::~NcollListFromFile() {}

void NcollListFromFile::Exec() {
  Clear();
  Jetscape::JSINFO << "Read binary collision list from file ...";
  try {
    std::string initialProfilePath =
        GetXMLElementText({"IS", "initial_Ncoll_list"});

    event_id_++;
    std::ostringstream path_with_filename;
    path_with_filename << initialProfilePath << "/event-" << event_id_
                       << "/NcollList.dat";
    JSINFO << "External initial profile path is" << path_with_filename.str();

    ReadNbcList(path_with_filename.str());
  } catch (std::exception &err) {
    Jetscape::JSWARN << err.what();
    std::exit(-1);
  }
}


void NcollListFromFile::Clear() {
  Jetscape::JSINFO << "clear initial condition vectors";
  binary_collision_x_.clear();
  binary_collision_y_.clear();
}


void NcollListFromFile::ReadNbcList(std::string filename) {
  Jetscape::JSINFO << "Read in binary collision list ...";
  std::ifstream infile(filename.c_str());
  if (!infile.good()) {
    Jetscape::JSWARN << "Can not open " << filename;
    exit(1);
  }

  std::string dummy;
  std::getline(infile, dummy);
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
  Jetscape::JSINFO << "done ...";
}


void NcollListFromFile::SampleABinaryCollisionPoint(double &x, double &y) {
  int rand_idx = (*rand_int_ptr_)(*GetMt19937Generator());
  x = binary_collision_x_[rand_idx];
  y = binary_collision_y_[rand_idx];
}
