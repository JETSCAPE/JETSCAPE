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

#ifndef JETSCAPEWRITERHEPMC_H
#define JETSCAPEWRITERHEPMC_H

#include <fstream>
#include <string>

#include "JetScapeWriter.h"
#include "PartonShower.h"

#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/Print.h"

// using namespace HepMC;
using HepMC3::GenEvent;
using HepMC3::GenVertex;
using HepMC3::GenParticle;
using HepMC3::GenVertexPtr;
using HepMC3::GenParticlePtr;

namespace Jetscape {

class JetScapeWriterHepMC : public JetScapeWriter, public HepMC3::WriterAscii {

public:
  JetScapeWriterHepMC() : HepMC3::WriterAscii("") { SetId("HepMC writer"); };
  JetScapeWriterHepMC(string m_file_name_out)
      : JetScapeWriter(m_file_name_out), HepMC3::WriterAscii(m_file_name_out) {
    SetId("HepMC writer");
  };
  virtual ~JetScapeWriterHepMC();

  void Init();
  void Exec();

  bool GetStatus() { return failed(); }
  void Close() { close(); }

  // // NEVER use this!
  // // Can work with only one writer, but with a second one it gets called twice
  // void WriteTask(weak_ptr<JetScapeWriter> w);

  // overload write functions
  void WriteEvent();

  // At parton level, we should never accept anything other than a full shower
  // void Write(weak_ptr<Vertex> v);
  void Write(weak_ptr<PartonShower> ps);
  void Write(weak_ptr<Hadron> h);
  void WriteHeaderToFile();

private:
  HepMC3::GenEvent evt;
  vector<HepMC3::GenVertexPtr> vertices;
  HepMC3::GenVertexPtr hadronizationvertex;

  inline HepMC3::GenVertexPtr
  castVtxToHepMC(const shared_ptr<Vertex> vtx) const {
    double x = vtx->x_in().x();
    double y = vtx->x_in().y();
    double z = vtx->x_in().z();
    double t = vtx->x_in().t();
    HepMC3::FourVector vtxPosition(x, y, z, t);
    // if ( t< 1e-6 ) t = 1e-6; // could do this. Exact 0 is bit quirky but works for hepmc
    return make_shared<GenVertex>(vtxPosition);
  }

  inline HepMC3::GenParticlePtr
  castPartonToHepMC(const shared_ptr<Parton> pparticle) const {
    return castPartonToHepMC(*pparticle);
  }

  inline HepMC3::GenParticlePtr
  castPartonToHepMC(const Parton &particle) const {
    HepMC3::FourVector pmom(particle.px(), particle.py(), particle.pz(),
                            particle.e());
    return make_shared<GenParticle>(pmom, particle.pid(), particle.pstat());
  }

  inline HepMC3::GenParticlePtr
  castHadronToHepMC(const shared_ptr<Hadron> pparticle) const {
    return castHadronToHepMC(*pparticle);
  }

  inline HepMC3::GenParticlePtr
  castHadronToHepMC(const Hadron &particle) const {
    HepMC3::FourVector pmom(particle.px(), particle.py(), particle.pz(),
                            particle.e());
    return make_shared<GenParticle>(pmom, particle.pid(), particle.pstat());
  }

  //int m_precision; //!< Output precision
};

} // end namespace Jetscape

#endif
