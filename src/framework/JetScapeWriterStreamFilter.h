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

// jetscape writer ascii class, filter

#ifndef JETSCAPEWRITERSTREAMFILTER_H
#define JETSCAPEWRITERSTREAMFILTER_H

#include "JetScapeWriterStream.h"
#include "sigslot.h"

#define JETSCAPEWRITER_PARTONSHOWER 1
#define JETSCAPEWRITER_PARTON 2
#define JETSCAPEWRITER_VERTEX 4
#define JETSCAPEWRITER_HADRON 8

namespace Jetscape {

template <class T>
class JetScapeWriterStreamFilter : public JetScapeWriterStream<T> {
 public:
  JetScapeWriterStreamFilter<T>(){};
  JetScapeWriterStreamFilter<T>(string m_file_name_out, unsigned char filter)
      : JetScapeWriterStream<T>(m_file_name_out) {
    displayFilter = filter;
  }
  virtual ~JetScapeWriterStreamFilter<T>(){};

  sigslot::signal1<vector<shared_ptr<Hadron>>&> GetHadronList;

  // void Init();
  // void Exec();
  //
  // bool GetStatus() {return output_file.good();}
  // void Close() {output_file.close();}
  //
  // void WriteInitFileXML();

  void Write(weak_ptr<PartonShower> ps) {
    if (displayFilter & JETSCAPEWRITER_PARTONSHOWER) {
      JetScapeWriterStream<T>::Write(ps);
    }
  }

  void Write(weak_ptr<Parton> p) {
    if (displayFilter & JETSCAPEWRITER_PARTON) {
      JetScapeWriterStream<T>::Write(p);
    }
  }

  void Write(weak_ptr<Vertex> v) {
    if (displayFilter & JETSCAPEWRITER_VERTEX) {
      JetScapeWriterStream<T>::Write(v);
    }
  }

  void Write(weak_ptr<Hadron> h) {
    if (displayFilter & JETSCAPEWRITER_HADRON) {
      JetScapeWriterStream<T>::Write(h);
    }
  }

  // void WriteHeaderToFile();
  //
  // void Write(string s) {output_file<<s<<endl;}
  // void WriteComment(string s) {output_file<<"# "<<s<<endl;}
  // void WriteWhiteSpace(string s) {output_file<<s<<" ";}
  // void WriteEvent();

 protected:
  unsigned char displayFilter;
};

typedef JetScapeWriterStreamFilter<ofstream> JetScapeWriterAsciiFilter;
#ifdef USE_GZIP
typedef JetScapeWriterStreamFilter<ogzstream> JetScapeWriterAsciiGZFilter;
#endif

}  // end namespace Jetscape

#endif  // JETSCAPEWRITERSTREAM_H
