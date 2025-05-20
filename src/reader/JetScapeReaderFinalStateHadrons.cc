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

#include "JetScapeReaderFinalStateHadrons.h"
#include <sstream>

namespace Jetscape {

JetScapeReaderFinalStateHadrons::~JetScapeReaderFinalStateHadrons() { VERBOSE(8); }

void JetScapeReaderFinalStateHadrons::ClearTask() {
  //pShower->clear();//pShower=nullptr; //check ...
  hadrons.clear();

  sigmaGen = -1;
  sigmaErr = -1;
  eventWeight = -1;
  EventPlaneAngle = 0.0;
}

void JetScapeReaderFinalStateHadrons::AddHadron(string s) {
  string token;
  strT.set(s);

  vector<string> vS;
  double x[4];
  x[0] = x[1] = x[2] = x[3] = 0.0;
  while (!strT.done()) {
    token = strT.next();
    if (token.compare("H") != 0)
      vS.push_back(token);
  }

  double p[4] = {stod(vS[3]), stod(vS[4]), stod(vS[5]), stod(vS[6])};
  double mass = sqrt(p[3]*p[3] - p[4]*p[4] - p[5]*p[5] - p[6]*p[6]);

  hadrons.push_back(make_shared<Hadron>(0, stoi(vS[1]), stoi(vS[2]), p, x, mass));
}

void JetScapeReaderFinalStateHadrons::Next() {
  if (currentEvent > 1)
    ClearTask();

  //ReadEvent(currentPos);
  string line;
  string token;

  JSINFO << "Current Event = " << currentEvent;

  std::string EPAngleStr = "EventPlaneAngle";
  while (getline(inFile, line)) {
    strT.set(line);

    if (strT.isCommentEntry()) {

      // Cross section
      if (line.find("sigmaGen") != std::string::npos) {
        std::stringstream data(line);
        std::string dummy;
        data >> dummy >> dummy >> sigmaGen >> dummy >> sigmaErr >> dummy >> eventWeight;
        JSDEBUG << " sigma gen=" << sigmaGen;
        JSDEBUG << " sigma err=" << sigmaErr;
        JSDEBUG << " Event weight=" << eventWeight;
      }
      // New Event
      if (strT.isEventEntry()) {
        strT.next(); strT.next();
        int newEvent = stoi(strT.next());
        strT.next(); strT.next(); strT.next();
        EventPlaneAngle = stod(strT.next());

        if (currentEvent != newEvent && currentEvent > -1) {
          currentEvent++;
          break;
        }
        currentEvent = newEvent;
      }
      continue;
    }


    // rest is list entry == hadron entry
    // Some questionable nomenclature here - identifying all that begins with "[" as "GraphEntry"
    // oh well
    AddHadron(line);
  }

  // different from normal since event line reader starts at 1 instead of 0
  if (Finished())
    currentEvent += 0;
}


vector<fjcore::PseudoJet> JetScapeReaderFinalStateHadrons::GetHadronsForFastJet() {
  vector<fjcore::PseudoJet> forFJ;

  for (auto &h : hadrons) {
    forFJ.push_back(h->GetPseudoJet());
  }

  return forFJ;
}

void JetScapeReaderFinalStateHadrons::InitTask() {
  VERBOSE(8) << "Open Input File = " << file_name_in;
  JSINFO << "Open Input File = " << file_name_in;

  inFile.open(file_name_in.c_str());

  if (!inFile.good()) {
    JSWARN << "Corrupt input file!";
    exit(-1);
  } else
    JSINFO << "File opened";

  currentEvent = 1;
}

int JetScapeReaderFinalStateHadrons::TotalEventCount() {
    std::string filename = file_name_in.c_str();
    std::string search_string = "Event";

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return 1;
    }

    std::string line;
    std::string last_occurrence_line;
    size_t last_occurrence_pos = std::string::npos;
    int line_number = 0;
    int last_occurrence_line_number = -1;

    while (std::getline(file, line)) {
        line_number++;
        size_t pos = line.rfind(search_string);
        if (pos != std::string::npos) {
            last_occurrence_pos = pos;
            last_occurrence_line = line;
            last_occurrence_line_number = line_number;
        }
    }
    file.close();

    int totalEvents = 0;
    if (last_occurrence_pos != std::string::npos) {
        std::stringstream data(last_occurrence_line);
        std::string dummy;
        data >> dummy >> dummy >> totalEvents;
    } else {
        std::cout << "String \"" << search_string << "\" not found in file." << std::endl;
    }

    return totalEvents;
}


} // end namespace Jetscape
