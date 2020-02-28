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

// generic jetscape writer base clase

#ifndef JETSCAPEWRITER_H
#define JETSCAPEWRITER_H

#include <string>
#include "JetScapeModuleBase.h"
#include "PartonShower.h"
#include "JetClass.h"
#include "JetScapeEventHeader.h"

using std::to_string;

namespace Jetscape {

class JetScapeWriter : public JetScapeModuleBase {

public:
  JetScapeWriter(){};
  JetScapeWriter(string m_file_name_out) { file_name_out = m_file_name_out; }
  virtual ~JetScapeWriter(){};

  void SetOutputFileName(string m_file_name_out) {
    file_name_out = m_file_name_out;
  }
  string GetOutputFileName() { return file_name_out; }

  virtual bool GetStatus() = 0;
  virtual void Close(){};
  virtual void Open(){};

  virtual void WriteInitFileXML(){};
  virtual void Write(weak_ptr<Parton> p){};
  virtual void Write(weak_ptr<Jet> j){};
  virtual void Write(weak_ptr<Vertex> v){};
  virtual void Write(weak_ptr<PartonShower> ps){};
  virtual void Write(string s){};
  virtual void WriteComment(string s){};
  virtual void WriteWhiteSpace(string s){};
  virtual void Write(ostream *o){};
  virtual void Write(weak_ptr<Hadron> h){};

  /// Gets called first, before all tasks write themselves
  virtual void WriteHeaderToFile(){};

  /// Gets called last, after all tasks have written themselves
  virtual void WriteEvent(){};

  virtual JetScapeEventHeader &GetHeader() { return header; };

protected:
  string file_name_out;
  JetScapeEventHeader header;
};

} // end namespace Jetscape

#endif
