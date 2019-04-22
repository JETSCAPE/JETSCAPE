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
// This is a wrapper for iS3D with the JETSCAPE framework
// -----------------------------------------

#include "JetScapeLogger.h"
#include "iS3DWrapper.h"

#include <string>

using namespace Jetscape;

iS3DWrapper::iS3DWrapper() 
{
    SetId("iS3D");
}

iS3DWrapper::~iS3DWrapper() 
{
}

void iS3DWrapper::InitTask() 
{
    JSINFO << "Initialize a particle sampler (iS3D)";
}

void iS3DWrapper::Exec() 
{
  // run_particlization(0) : read surface from memory
  // run_particlization(1) : read fo surface from file
  iS3D.run_particlization(1);
  PassHadronListToJetscape();
}

void iS3DWrapper::Clear() 
{
    VERBOSE(2) << "Finish iS3D particle sampling";
}

void iS3DWrapper::PassHadronListToJetscape()
{
  VERBOSE(2) << "Passing all sampled hadrons to the JETSCAPE framework";
  
  for (unsigned int iev = 0; iev < iS3D.final_particles_.size(); iev++)
    {
      std::vector<shared_ptr<Hadron>> hadrons;
      unsigned int nparticles = ( iS3D.final_particles_[iev].size() );
      //VERBOSE(4) << "number of particles = " << nparticles;
      for (unsigned int ipart = 0; ipart < nparticles; ipart++)
	{
	  Sampled_Particle current_hadron = iS3D.final_particles_[iev][ipart];
	  int hadron_label = 0;
	  int hadron_status = -1;
	  int hadron_id = current_hadron.mcID;
	  double hadron_mass = current_hadron.mass;
	  FourVector hadron_p(current_hadron.px, current_hadron.py, current_hadron.pz, current_hadron.E);
	  FourVector hadron_x(current_hadron.x, current_hadron.y, current_hadron.z, current_hadron.t);

	  // create a JETSCAPE Hadron
	  hadrons.push_back(make_shared<Hadron>(hadron_label, hadron_id, hadron_status, hadron_p, hadron_x, hadron_mass));
	  //Hadron* jetscape_hadron = new Hadron(hadron_label, hadron_id, hadron_status, hadron_p, hadron_x, hadron_mass);
	  //(*Hadron_list_)[iev]->push_back(*jetscape_hadron);
	}
      Hadron_list_.push_back(hadrons);
    }
  
  VERBOSE(4) << "JETSCAPE received " << Hadron_list_.size() << " events.";
  for (unsigned int iev = 0; iev < Hadron_list_.size(); iev++)
  {
    VERBOSE(4) << "In event " << iev << " JETSCAPE received " << Hadron_list_.at(iev).size() << " particles.";
  }
}

void iS3DWrapper::WriteTask(weak_ptr<JetScapeWriter> w)
{
  
  VERBOSE(4)<<"In iS3DWrapper::WriteTask";
  auto f = w.lock();
  if ( !f ) return;
  unsigned int nev = Hadron_list_.size();
  // write the number of samples as first line of header
  f->WriteComment( to_string( nev ) );
  f->WriteComment("JetScape module: " + GetId());
  if (nev > 0) {
    f->WriteComment("Final State Bulk Hadrons");
    for (unsigned int j = 0; j < nev; j++){
      vector<shared_ptr<Hadron>> hadVec = Hadron_list_.at(j);
      for (unsigned int i = 0; i < hadVec.size(); i++) {
	f->WriteWhiteSpace("["+to_string(i)+"] H");
	f->Write(hadVec.at(i));
      }
    }
  } else {
    f->WriteComment("There are no bulk Hadrons");
  }
 
}
