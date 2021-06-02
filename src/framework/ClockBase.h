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

//Remark JP: Think about if this base class is truly necessary ...
//Anyways, keep for now in case changes are needed to current implementation idea ...
//Like putting more clock functions in to base class ...

#ifndef CLOCKBASE_H
#define CLOCKBASE_H

#include <string>
#include <memory>

using std::string;

namespace Jetscape {

class ClockBase
{

public:

	ClockBase();
	virtual ~ClockBase() {};

	virtual void Info();

	void SetId(string m_id) { id = m_id; }
	void SetTimeRefFrameId(string m_time_id) { time_id = m_time_id; }
	//void SetCurrentTime(double m_CurrentTime) {currentTime = m_CurrentTime;}

	const string GetId() const { return id; }
	const string GetTimeRefFrameId() const { return time_id; }

	virtual double GetCurrentTime() {return -99.;}  //not clear if needed ...

	//static bool ClockUsed() { return use_clock; }

private:

	string id;
  	string time_id;

  	//static bool use_clock; //better in time based module base ...
};

} // end namespace Jetscape

#endif