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

#ifndef MILNECLOCK_H
#define MILNECLOCK_H

#include "ModuleClock.h"
#include "MainClock.h"
#include "RealType.h"
#include <memory>

using Jetscape::real;

namespace Jetscape {

class MilneClock : public ModuleClock {
 public:
    MilneClock();
    virtual ~MilneClock() {};

    void Info();
    void Transform(std::weak_ptr<MainClock> mainClock);

    void setEtaMax(real etaMax) {etaMax_ = etaMax;}
    real getTauMin() const {return(tauMin_);}
    real getTauMax() const {return(tauMax_);}

 private:
    real etaMax_;
    real tauMin_;
    real tauMax_;
};

}  // end namespace Jetscape

#endif  // MILNECLOCK_H
