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

#ifndef HADRONIZATION_H
#define HADRONIZATION_H

#include "JetScapeModuleBase.h"
#include "JetClass.h"
#include "JetScapeWriter.h"
#include "PartonShower.h"
//#include "MakeUniqueHelper.h"
#include <vector>
#include <random>

namespace Jetscape {

class Hadronization : public JetScapeModuleBase,
                      public std::enable_shared_from_this<Hadronization> {
public:
  Hadronization();
  virtual ~Hadronization();
  virtual shared_ptr<Hadronization> Clone() const { return nullptr; }
  virtual void Init();
  virtual void Exec();
  virtual void DoHadronization(vector<vector<shared_ptr<Parton>>> &pIn,
                               vector<shared_ptr<Hadron>> &hOut,
                               vector<shared_ptr<Parton>> &pOut){};
  virtual void WriteTask(weak_ptr<JetScapeWriter> w);
  virtual void Clear();

  void GetHadrons(vector<shared_ptr<Hadron>>& signal){signal = outHadrons;}
  sigslot::signal3<vector<vector<shared_ptr<Parton>>> &,
                   vector<shared_ptr<Hadron>> &, vector<shared_ptr<Parton>> &,
                   multi_threaded_local>
      TransformPartons;

  vector<shared_ptr<Hadron>> GetHadrons() { return outHadrons; }
  vector<shared_ptr<Parton>> GetOutPartons() { return outPartons; }

  void AddInPartons(vector<vector<shared_ptr<Parton>>> ip) { inPartons = ip; }

  void SetTransformPartonsConnected(bool m_TransformPartonsConnected) {
    TransformPartonsConnected = m_TransformPartonsConnected;
  }
  const bool GetTransformPartonsConnected() {
    return TransformPartonsConnected;
  }

  void AddInHadrons(vector<shared_ptr<Hadron>> ih) { outHadrons = ih; }

private:
  vector<vector<shared_ptr<Parton>>> inPartons;
  vector<shared_ptr<Hadron>> outHadrons;
  vector<shared_ptr<Parton>> outPartons;
  void DoHadronize();

  bool TransformPartonsConnected;
};

} // namespace Jetscape

#endif
