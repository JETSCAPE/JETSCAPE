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

// Jetscape final state {hadrons,kartons} writer ascii class
// Based on JetScapeWriterStream.
// author: Raymond Ehlers <raymond.ehlers@cern.ch>, ORNL

#ifndef JETSCAPEWRITERSTREAM_H
#define JETSCAPEWRITERSTREAM_H

#include <fstream>
#include <string>

#ifdef USE_GZIP
#include "gzstream.h"
#endif

#include "JetScapeWriter.h"

using std::ofstream;

namespace Jetscape {

template <class T>
class JetScapeWriterFinalStateStream : public JetScapeWriter {

public:
  JetScapeWriterFinalStateStream<T>(){};
  JetScapeWriterFinalStateStream<T>(string m_file_name_out);
  virtual ~JetScapeWriterFinalStateStream<T>();

  // This add "_final_state_*.dat" to the given filename.
  virtual void SetOutputFileName(string m_file_name_out);

  void Init();
  void Exec();

  virtual std::string GetName() { throw std::runtime_error("Don't use the base class"); }
  bool GetStatus() { return output_file.good(); }
  // Close is utilized to add the xsec and error.
  void Close();

  void Write(weak_ptr<PartonShower> ps);
  void Write(weak_ptr<Hadron> h);
  // We aren't interested in the individual partons or vertices, so skip them.

  void WriteHeaderToFile() { };
  void WriteEvent();

  void Write(string s) { output_file << s << endl; }
  // Intentionally make these no-ops since we want to fully control our output from this
  // class. Tasks will often directly call these functions, so we need to prevent them from doing so.
  void WriteComment(string s) { }
  void WriteWhiteSpace(string s) { }

protected:
  T output_file; //!< Output file
  std::vector<std::shared_ptr<JetScapeParticleBase>> particles;
};

template <class T>
class JetScapeWriterFinalStatePartonsStream : public JetScapeWriterFinalStateStream<T> {
  std::string GetName() { return "partons"; }
  // Don't collect the hadrons by making it a no-op
  void Write(weak_ptr<Hadron> h) { }
protected:
  // Allows the registration of the module so that it is available to be used by the Jetscape framework.
  static RegisterJetScapeModule<JetScapeWriterFinalStatePartonsStream<ofstream>> regParton;
  static RegisterJetScapeModule<JetScapeWriterFinalStatePartonsStream<ogzstream>> regPartonGZ;
};

template <class T>
class JetScapeWriterFinalStateHadronsStream : public JetScapeWriterFinalStateStream<T> {
  std::string GetName() { return "hadrons"; }
  // Don't collect the hadrons by making it a no-op
  void Write(weak_ptr<PartonShower> ps) { }
protected:
  // Allows the registration of the module so that it is available to be used by the Jetscape framework.
  static RegisterJetScapeModule<JetScapeWriterFinalStateHadronsStream<ofstream>> regHadron;
  static RegisterJetScapeModule<JetScapeWriterFinalStateHadronsStream<ogzstream>> regHadronGZ;
};

typedef JetScapeWriterFinalStatePartonsStream<ofstream> JetScapeWriterFinalStatePartonsAscii;
typedef JetScapeWriterFinalStateHadronsStream<ofstream> JetScapeWriterFinalStateHadronsAscii;
#ifdef USE_GZIP
typedef JetScapeWriterFinalStatePartonsStream<ogzstream> JetScapeWriterFinalStatePartonsAsciiGZ;
typedef JetScapeWriterFinalStateHadronsStream<ogzstream> JetScapeWriterFinalStateHadronsAsciiGZ;
#endif

} // end namespace Jetscape

#endif // JETSCAPEWRITERSTREAM_H
