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

// Jetscape Qn vector writer ascii class
// Based on JetScapeWriterStream, JetScapeWriterFinalStateStream
// author: Xiang-Yu Wu <xiangyu.wu2@mail.mcgill.ca>

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
class JetScapeWriterQnVectorStream : public JetScapeWriter {

public:
  JetScapeWriterQnVectorStream<T>(){};
  JetScapeWriterQnVectorStream<T>(string m_file_name_out);
  virtual ~JetScapeWriterQnVectorStream<T>();

  void Init();
  void Exec();

  virtual std::string GetName() { return "QnVector"; }
  bool GetStatus() { return output_file.good(); }
  // Close is utilized to add the xsec and error.
  void Close();

  void Write(weak_ptr<PartonShower> ps){ };
  void Write(weak_ptr<Hadron> h);
  // We aren't interested in the individual partons or vertices, so skip them.

  void WriteHeaderToFile() { };
  void WriteEvent();

  void Write(string s) { output_file << s << endl; }
  // Intentionally make these no-ops since we want to fully control our output from this
  // class. Tasks will often directly call these functions, so we need to prevent them from doing so.
  void WriteComment(string s) { }
  void WriteWhiteSpace(string s) { }

  void get_ptclist();
  int get_ch(const int pid);

protected:
  T output_file; //!< Output file
  std::vector<std::shared_ptr<Hadron>> particles;
  static RegisterJetScapeModule<JetScapeWriterQnVectorStream<ofstream>> regQnVector;
  static RegisterJetScapeModule<JetScapeWriterQnVectorStream<ogzstream>> regQnVectorGZ;
private:
  double pTmin_;
  double pTmax_;
  double rapmin_;
  double rapmax_;
  int npT_;
  int nrap_;
  int norder_;
  std::map< int, int > chpdg_;
  

};


typedef JetScapeWriterQnVectorStream<ofstream> JetScapeWriterQnVectorAscii;
#ifdef USE_GZIP
typedef JetScapeWriterQnVectorStream<ogzstream> JetScapeWriterQnVectorAsciiGZ;
#endif

} // end namespace Jetscape

#endif // JETSCAPEWRITERSTREAM_H
