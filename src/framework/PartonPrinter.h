/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * Modular, task-based framework
 * Intial Design: Joern Putschke, Kolja Kauder (Wayne State University)
 * For the full list of contributors see AUTHORS.

 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 * or via email to bugs.jetscape.org@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#ifndef PARTONPRINTER_H
#define PARTONPRINTER_H

#include "JetScapeModuleBase.h"
#include "PartonShower.h"
#include <vector>
#include <string>
#include<fstream>


namespace Jetscape {

class PartonPrinter : public JetScapeModuleBase
{

public:

PartonPrinter();
virtual ~PartonPrinter();

virtual void Init();
virtual void Exec() final;
virtual void Clear();
    std::ofstream dist_output; ///< the output stream where events are saved to file

void GetFinalPartons(shared_ptr<PartonShower> pShower, vector<shared_ptr<Parton>>& fPartons);

void GetFinalPartons2(shared_ptr<PartonShower> pShower);

void GetPartonsAtTime(shared_ptr<PartonShower> pShower, vector<shared_ptr<Parton>>& fPartons, double time);

void PrintFinalPartons(vector<vector<shared_ptr<Parton>>>& fPartons) {fPartons = pFinals;};

private:

vector<vector<shared_ptr<Parton>>> pFinals;

};

} // end namespace Jetscape

#endif
