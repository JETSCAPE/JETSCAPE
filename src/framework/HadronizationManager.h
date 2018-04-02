/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * For the full list of contributors see AUTHORS.
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 * or via email to bugs.jetscape.org@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#ifndef HADRONIZATIONMANAGER_H
#define HADRONIZATIONMANAGER_H

#include "JetScapeTask.h"
#include "JetClass.h"
#include "sigslot.h"

#include <vector>

namespace Jetscape {

class HadronizationManager : public JetScapeTask, public std::enable_shared_from_this<HadronizationManager>
{
  
 public:
  
  HadronizationManager();
  virtual ~HadronizationManager();
  
  virtual void Init();
  virtual void Exec();
  virtual void Clear();
  virtual void WriteTask(weak_ptr<JetScapeWriter> w);
  
  int GetNumSignals();
  
  void CreateSignalSlots();

  //sigslot::signal1<vector<shared_ptr<Hadron>>& > GetHadronList;

  sigslot::signal1<vector<vector<shared_ptr<Parton>>>& > GetFinalPartonList;

  void SetGetFinalPartonListConnected(bool m_GetFinalPartonListConnected) {GetFinalPartonListConnected=m_GetFinalPartonListConnected;}
  const bool GetGetFinalPartonListConnected() {return GetFinalPartonListConnected;}

 private:

  bool GetFinalPartonListConnected;
  vector<vector<shared_ptr<Parton>>> hd;
  
};

} // end namespace Jetscape

#endif
