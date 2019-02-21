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

#ifndef OTHERPARTICLESSTORAGE_H
#define OTHERPARTICLESSTORAGE_H

#include "JetScapeTask.h"
#include "JetScapeParticles.h"
#include "sigslot.h"

namespace Jetscape {

class OtherParticlesStorage : public JetScapeTask {

  public:
    OtherParticlesStorage();
    ~OtherParticlesStorage();

    void Init();
    void Clear();
    void WriteTask(weak_ptr<JetScapeWriter> w);


    

  private:
    bool GetBosonListConnected;
    bool GetLeptonListConnected;
    bool GetPhotonListConnected;
    bool GetElectronListConnected;
    bool GetMuonListConnected;
    bool GetPositronListConnected;
    vector<shared_ptr<Boson>> bosons;
    vector<shared_ptr<Lepton>> leptons;
    vector<shared_ptr<Photon>> photons;
    vector<shared_ptr<Electron>> electrons;
    vector<shared_ptr<Muon>> muons;
    vector<shared_ptr<Positron>> positrons;

};

};

#endif //OTHERPARTICLESSTORAGE_H
