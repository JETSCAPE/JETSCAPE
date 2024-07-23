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
// -----------------------------------------
// This is a wrapper for iSpectraSampler (iSS) with the JETSCAPE framework
// -----------------------------------------

#include "JetScapeLogger.h"
#include "EMProbeWrapper.h"

#include "FluidDynamics.h"
#include "FluidCellInfo.h"
#include "SurfaceCellInfo.h"
#include "FluidEvolutionHistory.h"
#include "JetScapeSignalManager.h"



#include <memory>
#include <string>
#include <fstream>

using namespace Jetscape;

// Register the module with the base class
RegisterJetScapeModule<EMProbeWrapper> EMProbeWrapper::reg("Dilepton");

EMProbeWrapper::EMProbeWrapper() {
    SetId("Dilepton");
    statusCode_ = 0;
}

EMProbeWrapper::~EMProbeWrapper() {}

void EMProbeWrapper::InitTask() {

   JSINFO << "Initialize EM probe";
   bulk_info_array.clear();
   std::string input_file =
       GetXMLElementText({"EMProbe", "Dilepton", "Dilepton_input_file"});
   double dilepton_alpha_s =
       GetXMLElementDouble({"EMProbe", "Dilepton", "alpha_s"});
   std::string dilepton_flag_hydro =
       GetXMLElementText({"EMProbe", "Dilepton", "flag_hydro"});
   std::string dilepton_flag_prehydro =
       GetXMLElementText({"EMProbe", "Dilepton", "flag_prehydro"});
   std::string include_baryondiff_deltaf =
       GetXMLElementText({"EMProbe", "Dilepton", "include_baryondiff_deltaf"});
   std::string include_shearvisc_deltaf =
       GetXMLElementText({"EMProbe", "Dilepton", "include_shearvisc_deltaf"});
   std::string turn_on_muB =
       GetXMLElementText({"EMProbe", "Dilepton", "turn_on_muB"});


  Dilepton_ptr_ = std::unique_ptr<Dilepton::JS_dilepton>(
           new Dilepton::JS_dilepton(input_file));
  Dilepton_ptr_->paraRdr->setVal("alpha_s",dilepton_alpha_s);
  Dilepton_ptr_->paraRdr->setVal("flag_hydro",1);
  Dilepton_ptr_->paraRdr->setVal("flag_prehydro",0);
  Dilepton_ptr_->paraRdr->setVal("include_baryondiff_deltaf",0);
  Dilepton_ptr_->paraRdr->setVal("include_shearvisc_deltaf",0);
  Dilepton_ptr_->paraRdr->setVal("turn_on_muB",0);
}

void EMProbeWrapper::Exec() {
   JSINFO << "running Dilepton ...";
   bulk_info_array.resize(0);
   getBulkInforfromJetScape();
   dilepton_spec = std::make_shared<std::vector<float>>(Dilepton_ptr_->run(bulk_info_array));

   JSINFO << "Dilepton calculation finished.";
}

void EMProbeWrapper::Clear() {

}


void EMProbeWrapper::WriteTask(weak_ptr<JetScapeWriter> w) {
   VERBOSE(4) << "In EMProbeWrapper::WriteTask";
   auto f = w.lock();
   if (!f)
     return;

   f->WriteComment("JetScape module: " + GetId());
   f->WriteWhiteSpace("[" + to_string(dilepton_spec->size()) + "] EM ");
   f->Write(dilepton_spec);
}


void EMProbeWrapper::getBulkInforfromJetScape() {
   std::shared_ptr<FluidDynamics> hydro_ptr = JetScapeSignalManager::Instance()->GetHydroPointer().lock();
   const EvolutionHistory& bulk_info = hydro_ptr->get_bulk_info();
   int turn_on_rhob = 0;
   int turn_on_shear = 0;
   int turn_on_bulk = 0;
   int turn_on_diff = 0;

   const int nVar_per_cell = (11 + turn_on_rhob *2 + turn_on_shear*5
                                 + turn_on_bulk*1 +  turn_on_diff*3);
   int ncell_above_T0 = 0;
   bulk_info_array.push_back(bulk_info.tau_min);
   bulk_info_array.push_back(bulk_info.dtau);
   bulk_info_array.push_back(bulk_info.nx);
   bulk_info_array.push_back(bulk_info.dx);
   bulk_info_array.push_back(bulk_info.x_min);
   bulk_info_array.push_back(bulk_info.ny);
   bulk_info_array.push_back(bulk_info.dy);
   bulk_info_array.push_back(bulk_info.y_min);
   bulk_info_array.push_back(bulk_info.neta);
   bulk_info_array.push_back(bulk_info.deta);
   bulk_info_array.push_back(bulk_info.eta_min);
   
   bulk_info_array.push_back(turn_on_rhob);
   bulk_info_array.push_back(turn_on_shear);
   bulk_info_array.push_back(turn_on_bulk);
   bulk_info_array.push_back(turn_on_diff);
   bulk_info_array.push_back(nVar_per_cell);
   float hbarc = 0.19732;
   for(int icell = 0; icell < bulk_info.data.size(); icell++ ){
        
    int id_eta = icell % bulk_info.neta;
    int id_y = (icell / bulk_info.neta) % bulk_info.ny;
    int id_x = (icell / (bulk_info.neta * bulk_info.ny)) % bulk_info.nx;
    int id_tau = icell / (bulk_info.neta * bulk_info.ny * bulk_info.nx);
    float Temp_tem = bulk_info.data[icell].temperature;
    if( Temp_tem > 0.1){
        
        bulk_info_array.push_back(id_tau);
        bulk_info_array.push_back(id_y);
        bulk_info_array.push_back(id_x);
        bulk_info_array.push_back(id_eta);
        bulk_info_array.push_back(bulk_info.data[icell].energy_density);
        bulk_info_array.push_back(bulk_info.data[icell].pressure);
        bulk_info_array.push_back(bulk_info.data[icell].temperature);
        bulk_info_array.push_back(0.0);
        float vx = bulk_info.data[icell].vx;
        float vy = bulk_info.data[icell].vy;
        float vz = bulk_info.data[icell].vz;

        float gamma = 1.0/sqrt(1.0 - (vx*vx+vy*vy+vz*vz));
        float u0 = gamma;
        float ux = gamma*vx;
        float uy = gamma*vy;
        float uz = gamma*vz;

        float icell_eta = bulk_info.EtaCoord(id_eta);
        float cosheta = cosh(icell_eta);
        float sinheta = sinh(icell_eta);

        
        float utau = cosheta*u0 - sinheta*uz;
        float ueta = -sinheta*u0 + cosheta*uz;

        bulk_info_array.push_back(ux);
        bulk_info_array.push_back(uy);
        bulk_info_array.push_back(ueta);
        ncell_above_T0++;
        

        
    
    }
      
   }
   


}
