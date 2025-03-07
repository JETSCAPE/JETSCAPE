/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion
 *collisions
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

#include "HepMC3/GenEvent.h"
#include "HepMC3/Print.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "JetScapeWriter.h"
#include "PartonShower.h"

// using namespace HepMC;
using HepMC3::GenEvent;
using HepMC3::GenParticle;
using HepMC3::GenParticlePtr;
using HepMC3::GenVertex;
using HepMC3::GenVertexPtr;

namespace Jetscape {

/**
 * @class JetScapeWriterHepMC
 * @brief A class for writing JetScape events to HepMC format.
 *
 * This class inherits from JetScapeWriter and HepMC3::WriterAscii to provide
 * functionality for writing JetScape events to a HepMC file.
 */
class JetScapeWriterHepMC : public JetScapeWriter, public HepMC3::WriterAscii {
 public:
  /**
   * @brief Default constructor.
   *
   * Initializes the HepMC writer with an empty file name and sets the writer
   * ID.
   */
  JetScapeWriterHepMC() : HepMC3::WriterAscii("") { SetId("HepMC writer"); };

  /**
   * @brief Constructor with file name.
   *
   * @param m_file_name_out The output file name for the HepMC writer.
   *
   * Initializes the JetScapeWriter and HepMC writer with the specified file
   * name and sets the writer ID.
   */
  JetScapeWriterHepMC(string m_file_name_out)
      : JetScapeWriter(m_file_name_out), HepMC3::WriterAscii(m_file_name_out) {
    SetId("HepMC writer");
  };

  /**
   * @brief Destructor.
   *
   * Cleans up resources used by the JetScapeWriterHepMC.
   */
  virtual ~JetScapeWriterHepMC();

  /**
   * @brief Initializes the writer.
   *
   * This function should be called to initialize any necessary resources before
   * writing events.
   */
  void Init();

  /**
   * @brief Executes the writer.
   *
   * @note This function currently does not perform any operations.
   */
  void Exec();

  /**
   * @brief Gets the status of the writer.
   *
   * @return True if the writer has failed, false otherwise.
   */
  bool GetStatus() { return failed(); }

  /**
   * @brief Closes the writer.
   *
   * This function should be called to close the writer and release any
   * resources.
   */
  void Close() { close(); }

  // NEVER use this!
  // Can work with only one writer, but with a second one it gets
  // called twice
  // void WriteTask(weak_ptr<JetScapeWriter> w);

  // overload write functions

  /**
   * @brief Writes an event to the HepMC file.
   *
   * This function should be called to write a single event to the HepMC file.
   */
  void WriteEvent();

  // At parton level, we should never accept anything other than a full shower
  // void Write(weak_ptr<Vertex> v);

  /**
   * @brief Writes a parton shower to the HepMC file.
   *
   * @param ps A weak pointer to the PartonShower to be written.
   */
  void Write(weak_ptr<PartonShower> ps);

  /**
   * @brief Writes a hadron to the HepMC file.
   *
   * @param h A weak pointer to the Hadron to be written.
   */
  void Write(weak_ptr<Hadron> h);

  /**
   * @brief Writes the header to the HepMC file.
   *
   * This function should be called to write the header information to
   * the HepMC file.
   */
  void WriteHeaderToFile();

 private:
  HepMC3::GenEvent evt;
  vector<HepMC3::GenVertexPtr> vertices;
  HepMC3::GenVertexPtr hadronizationvertex;

  /**
   * @note WriteEvent needs to know whether it should overwrite final
   * partons status to 1
   */
  bool hashadrons = false;

  /**
   * @brief Casts a JetScape Vertex to a HepMC GenVertex.
   *
   * @param vtx A shared pointer to the Vertex to be cast.
   * @return A shared pointer to the corresponding HepMC GenVertex.
   */
  inline HepMC3::GenVertexPtr castVtxToHepMC(
      const shared_ptr<Vertex> vtx) const {
    double x = vtx->x_in().x();
    double y = vtx->x_in().y();
    double z = vtx->x_in().z();
    double t = vtx->x_in().t();
    HepMC3::FourVector vtxPosition(x, y, z, t);
    // if ( t< 1e-6 ) t = 1e-6; // could do this. Exact 0 is bit quirky but
    // works for hepmc
    return make_shared<GenVertex>(vtxPosition);
  }

  /**
   * @brief Casts a JetScape Parton to a HepMC GenParticle.
   *
   * @param pparticle A shared pointer to the Parton to be cast.
   * @return A shared pointer to the corresponding HepMC GenParticle.
   */
  inline HepMC3::GenParticlePtr castPartonToHepMC(
      const shared_ptr<Parton> pparticle) const {
    return castPartonToHepMC(*pparticle);
  }

  /**
   * @brief Casts a JetScape Parton to a HepMC GenParticle.
   *
   * @param particle A reference to the Parton to be cast.
   * @return A shared pointer to the corresponding HepMC GenParticle.
   */
  inline HepMC3::GenParticlePtr castPartonToHepMC(
      const Parton &particle) const {
    HepMC3::FourVector pmom(particle.px(), particle.py(), particle.pz(),
                            particle.e());
    return make_shared<GenParticle>(pmom, particle.pid(), particle.pstat());
  }

  /**
   * @brief Casts a JetScape Hadron to a HepMC GenParticle.
   *
   * @param pparticle A shared pointer to the Hadron to be cast.
   * @return A shared pointer to the corresponding HepMC GenParticle.
   */
  inline HepMC3::GenParticlePtr castHadronToHepMC(
      const shared_ptr<Hadron> pparticle) const {
    return castHadronToHepMC(*pparticle);
  }

  /**
   * @brief Casts a JetScape Hadron to a HepMC GenParticle.
   *
   * @param particle A reference to the Hadron to be cast.
   * @return A shared pointer to the corresponding HepMC GenParticle.
   */
  inline HepMC3::GenParticlePtr castHadronToHepMC(
      const Hadron &particle) const {
    HepMC3::FourVector pmom(particle.px(), particle.py(), particle.pz(),
                            particle.e());
    return make_shared<GenParticle>(pmom, particle.pid(), particle.pstat());
  }

  // int m_precision; //!< Output precision
};

}  // end namespace Jetscape

#endif
