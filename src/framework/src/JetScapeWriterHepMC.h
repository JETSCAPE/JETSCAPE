// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#ifndef JETSCAPEWRITERHEPMC_H
#define JETSCAPEWRITERHEPMC_H

#include <fstream>
#include <string>

#include "JetScapeWriter.h"
#include "PartonShower.h"

#undef DEBUG
#include "HepMC/GenEvent.h"
#include "HepMC/ReaderAscii.h"
#include "HepMC/WriterAscii.h"
#include "HepMC/Print.h"
#undef DEBUG

using namespace HepMC;

namespace Jetscape {

class JetScapeWriterHepMC : public JetScapeWriter , public WriterAscii
{

 public:

 JetScapeWriterHepMC() : WriterAscii("") { SetId("HepMC writer"); };
 JetScapeWriterHepMC(string m_file_name_out) : JetScapeWriter(m_file_name_out), WriterAscii(m_file_name_out) { SetId("HepMC writer"); };
  virtual ~JetScapeWriterHepMC();

  void Init();
  void Exec();
  void WriteTask(weak_ptr<JetScapeWriter> w);

  bool GetStatus() {return failed();}
  void Close() {close();}

  // overload write functions ...
  void WriteEvent(); 
  void Write(weak_ptr<Vertex> v);
  void Write(weak_ptr<PartonShower> ps);

 private:

  bool vertexFlag;
  vector<HepMC::GenVertex*> vertices;

  inline HepMC::GenVertex* castVtxToHepMC(shared_ptr<Vertex> vtx){
      HepMC::FourVector vtxPosition(vtx->x_in().x(), vtx->x_in().y(), vtx->x_in().z(), vtx->x_in().t());
      HepMC::GenVertex *hepVtx = new HepMC::GenVertex(vtxPosition);
      return hepVtx;
  }

  inline HepMC::GenParticle* castParticleToHepMC(Parton &particle){
      HepMC::FourVector pmom(particle.px(), particle.py(), particle.pz(), particle.e());
      HepMC::GenParticle *p1 = new HepMC::GenParticle(pmom, particle.pid(), particle.pstat());
      return p1;
  }

  inline HepMC::GenParticle* castParticleToHepMC(shared_ptr<Parton> particle){
      HepMC::FourVector pmom(particle->px(), particle->py(), particle->pz(), particle->e());
      HepMC::GenParticle *p1 = new HepMC::GenParticle(pmom, particle->pid(), particle->pstat());
      return p1;
  }
  //int m_precision; //!< Output precision
  
};


} // end namespace Jetscape

#endif
