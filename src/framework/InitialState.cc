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

#include "InitialState.h"
#include "JetScapeWriter.h"
#include <iostream>

namespace Jetscape {

InitialState::~InitialState() {}

void InitialState::Init() {
  JetScapeModuleBase::Init();

  JSINFO << "Initialize InitialState ... " << GetId() << " ...";

  grid_max_x_ = GetXMLElementDouble({"IS", "grid_max_x"});
  grid_max_y_ = GetXMLElementDouble({"IS", "grid_max_y"});
  grid_max_z_ = GetXMLElementDouble({"IS", "grid_max_z"});
  grid_step_x_ = GetXMLElementDouble({"IS", "grid_step_x"});
  grid_step_y_ = GetXMLElementDouble({"IS", "grid_step_y"});
  grid_step_z_ = GetXMLElementDouble({"IS", "grid_step_z"});
  JSINFO << "x range for bulk evolution = [" << -grid_max_x_ << ", "
         << grid_max_x_ << "]";

  InitTask();

  JetScapeTask::InitTasks();
}

void InitialState::Exec() {
  // Do whatever is needed to figure out the internal temp...
}

void InitialState::Clear() {}

void InitialState::Write(weak_ptr<JetScapeWriter> w) {
  //Write out the original vertex so the writer can keep track of it...
  // auto f = w.lock();
  // if ( f ) f->Write(make_shared<Vertex>(initialVtx));
}

void InitialState::CollectHeader(weak_ptr<JetScapeWriter> w) {
  auto f = w.lock();
  if (f) {
    auto &header = f->GetHeader();
    header.SetNpart(GetNpart());
    header.SetNcoll(GetNcoll());
    header.SetEventCentrality(GetEventCentrality());
    header.SetTotalEntropy(GetTotalEntropy());
  }
}

std::tuple<double, double, double> InitialState::CoordFromIdx(int idx) {
  int nx = GetXSize();
  int ny = GetYSize();
  int nz = GetZSize();

  int ix = idx / (ny * nz);
  int iy = (idx - (ny * nz * ix))/ nz;
  int ieta = idx - (ny * nz * ix) - (nz * iy); 

  return std::make_tuple(-grid_max_x_ + ix * grid_step_x_,
                         -grid_max_y_ + iy * grid_step_y_,
                         -grid_max_z_ + ieta * grid_step_z_);
}


void InitialState::SampleABinaryCollisionPoint(double &x, double &y) {
  if (num_of_binary_collisions_.size() == 0) {
    JSWARN << "num_of_binary_collisions is empty, setting the starting "
              "location to 0. Make sure to add e.g. trento before PythiaGun.";
  } else {
    std::discrete_distribution<> dist(
        begin(num_of_binary_collisions_),
        end(num_of_binary_collisions_)); // Create the distribution
    // Now generate values
    auto idx = dist(*GetMt19937Generator());
    auto coord = CoordFromIdx(idx);
    x = std::get<0>(coord);
    y = std::get<1>(coord);
  }
}

} // end namespace Jetscape
