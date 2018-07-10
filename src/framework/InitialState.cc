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

InitialState::InitialState(){
    //FourVector v(0,0,0,0);
    //initialVtx.set_location(v);
    //Init();
}

InitialState::~InitialState()
{
}

void InitialState::Init(){
  JetScapeModuleBase::Init();

  JSINFO<<"Intialize InitialState ... " << GetId() << " ...";

  // JetScapeXML::Instance() returns a singleton if jetscape_init.xml is read in.
  xml_ = JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("IS" );  

  if ( !xml_ ) {
    JSWARN << " : Not a valid JetScape Initial State XML section in file!";
    exit(-1);
  } else {
    xml_->FirstChildElement("grid_max_x")->QueryDoubleText(&grid_max_x_);
    xml_->FirstChildElement("grid_max_y")->QueryDoubleText(&grid_max_y_);
    xml_->FirstChildElement("grid_max_z")->QueryDoubleText(&grid_max_z_);
    xml_->FirstChildElement("grid_step_x")->QueryDoubleText(&grid_step_x_);
    xml_->FirstChildElement("grid_step_y")->QueryDoubleText(&grid_step_y_);
    xml_->FirstChildElement("grid_step_z")->QueryDoubleText(&grid_step_z_);
    JSINFO<<"x range for bulk evolution = ["<< -grid_max_x_ <<", "<<grid_max_x_ << "]";
  }

  InitTask();

  JetScapeTask::InitTasks();
}

void InitialState::Exec(){
    // Do whatever is needed to figure out the internal temp...
    
}

void InitialState::Clear(){
}

void InitialState::Write(weak_ptr<JetScapeWriter> w){
  //Write out the original vertex so the writer can keep track of it...
  // auto f = w.lock();
  // if ( f ) f->Write(make_shared<Vertex>(initialVtx));
}

void InitialState::CollectHeader( weak_ptr<JetScapeWriter> w ){
  auto f = w.lock();
  if ( f ){
    auto& header = f->GetHeader();
    header.SetNpart( GetNpart() );
    header.SetNcoll( GetNcoll() );
    header.SetTotalEntropy( GetTotalEntropy() );
  }
}

std::tuple<double, double, double> InitialState::CoordFromIdx(int idx) {
    int nx = GetXSize();
    int ny = GetYSize();
    int nz = GetZSize();

    int page = idx / (nx * ny);
    int row = (idx - page * nx * ny) / nx;
    int col = idx - page * nx * ny - row * nx;

    return std::make_tuple(-grid_max_x_ + col * grid_step_x_,
                           -grid_max_y_ + row * grid_step_y_,
                           -grid_max_z_ + page * grid_step_z_);
}

} // end namespace Jetscape
