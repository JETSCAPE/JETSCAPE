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
// -----------------------------------------
// This is a wrapper for iS3D with the JETSCAPE framework
// -----------------------------------------

#ifndef IS3DWRAPPER_H
#define IS3DWRAPPER_H

#include "SoftParticlization.h"
#include "iS3D.h"

using namespace Jetscape;

class iS3DWrapper: public SoftParticlization {
 private:
    iS3D::IS3D iS3D_;

    int freezeoutSurfaceFromFile_;

  // Allows the registration of the module so that it is available to be used by the Jetscape framework.
  static RegisterJetScapeModule<iS3DWrapper> reg;

 public:
    iS3DWrapper();
    ~iS3DWrapper();

    void InitTask();
    void Exec();
    void Clear();
    void WriteTask(weak_ptr<JetScapeWriter> w);

    void PassHadronListToJetscape();
};

#endif  // IS3DWRAPPER_H
