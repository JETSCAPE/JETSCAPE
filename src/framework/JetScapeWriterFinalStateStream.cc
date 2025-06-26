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
// Jetscape final state {hadrons,partons} writer ascii class
// Based on JetScapeWriterStream.
// author: Raymond Ehlers <raymond.ehlers@cern.ch>, LBL/UCB

#include "JetScapeWriterFinalStateStream.h"

#include "JetScapeLogger.h"
#include "JetScapeXML.h"

namespace Jetscape {

namespace detail {

std::vector<int> stringToVector(const std::string& str) {
  std::vector<int> result;
  std::stringstream ss(str);
  std::string item;

  while (std::getline(ss, item, ',')) {
    // Trim whitespace before and after for safety
    item.erase(0, item.find_first_not_of(" \t"));
    item.erase(item.find_last_not_of(" \t") + 1);

    // Store if there are any characters left.
    // NOTE: There's no validation for e.g. if characters other than numbers are
    // provided.
    if (!item.empty()) {
      result.push_back(std::stoi(item));
    }
  }

  return result;
}

};  // end namespace detail

// Register the modules with the base class
template <>
RegisterJetScapeModule<JetScapeWriterFinalStatePartonsStream<ofstream>>
    JetScapeWriterFinalStatePartonsStream<ofstream>::regParton(
        "JetScapeWriterFinalStatePartonsAscii");
template <>
RegisterJetScapeModule<JetScapeWriterFinalStateHadronsStream<ofstream>>
    JetScapeWriterFinalStateHadronsStream<ofstream>::regHadron(
        "JetScapeWriterFinalStateHadronsAscii");
template <>
RegisterJetScapeModule<JetScapeWriterFinalStatePartonsStream<ogzstream>>
    JetScapeWriterFinalStatePartonsStream<ogzstream>::regPartonGZ(
        "JetScapeWriterFinalStatePartonsAsciiGZ");
template <>
RegisterJetScapeModule<JetScapeWriterFinalStateHadronsStream<ogzstream>>
    JetScapeWriterFinalStateHadronsStream<ogzstream>::regHadronGZ(
        "JetScapeWriterFinalStateHadronsAsciiGZ");

template <class T>
JetScapeWriterFinalStateStream<T>::JetScapeWriterFinalStateStream(
    string m_file_name_out)
    : particles{},
      writeCentrality{false},
      writePtHat{false},
      particleStatusToSkip{},
      output_file{},
      headerVersion{2},
      binaryOutput{false} {
  SetOutputFileName(m_file_name_out);
}

template <class T>
JetScapeWriterFinalStateStream<T>::~JetScapeWriterFinalStateStream() {
  VERBOSE(8);
  if (GetActive())
    Close();
}

template <class T>
void JetScapeWriterFinalStateStream<T>::WriteEvent() {
  // Write the entire event all at once.
  if (!binaryOutput) {
    // Optionally write event centrality to event header
    std::string centrality_text = "";
    if (writeCentrality) {
      centrality_text += "\tcentrality\t";
      centrality_text += std::to_string(GetHeader().GetEventCentrality());
    }

    // Optionally write pt-hat value to event header
    std::string pt_hat_text = "";
    if (writePtHat) {
      pt_hat_text += "\tpt_hat\t";
      pt_hat_text += std::to_string(GetHeader().GetPtHat());
    }

    // We cannot write directly to the file since we don't know how many
    // particles we have in the filtered case. To resolve this, we'll write to
    // an intermediate stringstream, and then write that to file.
    std::stringstream ss;

    // Next, write the particles. Will contain either hadrons or partons based
    // on the derived class.
    unsigned int ipart = 0;
    for (const auto& p : particles) {
      auto particle = p.get();
      // Skip particles with requested status codes (if enabled).
      if (particleStatusToSkip.size() > 0) {
        // Skip particles with status codes that are in the list to skip
        if (std::find(particleStatusToSkip.begin(), particleStatusToSkip.end(),
                      particle->pstat()) != particleStatusToSkip.end()) {
          continue;
        }
      }

      ss << ipart
         // << " " << particle->plabel()
         << " " << particle->pid() << " " << particle->pstat() << " "
         << particle->e() << " " << particle->px() << " " << particle->py()
         << " " << particle->pz() << "\n";
      ++ipart;
    }

    // Now that we've filtered the particles, we can finally write to the output
    // file. First, write header NOTE: Needs consistent "\t" between all entries
    // to simplify parsing later. NOTE: Could also add Npart, Ncoll, and
    // TotalEntropy. See the original stream writer. NOTE: ipart == total number
    // of (filtered) particles by the end of the loop, so we can just take
    // advantage of it.
    output_file << "#"
                << "\t"
                << "Event\t"
                << GetCurrentEvent() + 1  // +1 to index the event count from 1
                << "\t"
                << "weight\t" << std::setprecision(15)
                << GetHeader().GetEventWeight() << std::setprecision(6) << "\t"
                << "EPangle\t"
                << (GetHeader().GetEventPlaneAngle() > -999
                        ? GetHeader().GetEventPlaneAngle()
                        : 0)
                << "\t"
                << "N_" << GetName() << "\t" << ipart;
    if (headerVersion == 3) {
      output_file << "\t"
                  << "vertex_x\t" << GetHeader().GetVertexX() << "\t"
                  << "vertex_y\t" << GetHeader().GetVertexY() << "\t"
                  << "vertex_z\t" << GetHeader().GetVertexZ();
    }
    output_file << centrality_text << pt_hat_text << "\n";

    // Finally, write the particles
    output_file << ss.str();
  } else {
    const int event_number = GetCurrentEvent() + 1;
    const double event_weight =
        static_cast<double>(GetHeader().GetEventWeight());
    const float event_plane_angle =
        static_cast<float>(GetHeader().GetEventPlaneAngle() > -999
                               ? GetHeader().GetEventPlaneAngle()
                               : 0);
    const float vertex_x = static_cast<float>(GetHeader().GetVertexX());
    const float vertex_y = static_cast<float>(GetHeader().GetVertexY());
    const float vertex_z = static_cast<float>(GetHeader().GetVertexZ());
    const float centrality =
        static_cast<float>(GetHeader().GetEventCentrality());
    const float pt_hat = static_cast<float>(GetHeader().GetPtHat());

    // Write the event number, weight, and event plane angle, ...
    output_file.write(reinterpret_cast<const char*>(&event_number),
                      sizeof(int));
    output_file.write(reinterpret_cast<const char*>(&event_weight),
                      sizeof(double));
    output_file.write(reinterpret_cast<const char*>(&event_plane_angle),
                      sizeof(float));
    output_file.write(reinterpret_cast<const char*>(&vertex_x), sizeof(float));
    output_file.write(reinterpret_cast<const char*>(&vertex_y), sizeof(float));
    output_file.write(reinterpret_cast<const char*>(&vertex_z), sizeof(float));
    output_file.write(reinterpret_cast<const char*>(&centrality),
                      sizeof(float));
    output_file.write(reinterpret_cast<const char*>(&pt_hat), sizeof(float));

    // Create a vector of int to store the pid, pstat values
    std::vector<int> pid_values;
    std::vector<int> pstat_values;
    // Create a vector of float to store the e, px, py, pz values
    std::vector<float> e_values;
    std::vector<float> px_values;
    std::vector<float> py_values;
    std::vector<float> pz_values;
    for (const auto& p : particles) {
      auto particle = p.get();
      // Skip particles with requested status codes (if enabled).
      if (particleStatusToSkip.size() > 0) {
        // Skip particles with status codes that are in the list to skip
        if (std::find(particleStatusToSkip.begin(), particleStatusToSkip.end(),
                      particle->pstat()) != particleStatusToSkip.end()) {
          continue;
        }
      }

      // Store the values in the vectors
      pid_values.push_back(particle->pid());
      pstat_values.push_back(particle->pstat());
      e_values.push_back(static_cast<float>(particle->e()));
      px_values.push_back(static_cast<float>(particle->px()));
      py_values.push_back(static_cast<float>(particle->py()));
      pz_values.push_back(static_cast<float>(particle->pz()));
    }

    // Write the number of particles
    const int num_particles = static_cast<int>(pid_values.size());
    output_file.write(reinterpret_cast<const char*>(&num_particles),
                      sizeof(int));
    // Write the pid, pstat, e, px, py, pz values for each particle
    for (int i = 0; i < num_particles; ++i) {
      output_file.write(reinterpret_cast<const char*>(&pid_values[i]),
                        sizeof(int));
      output_file.write(reinterpret_cast<const char*>(&pstat_values[i]),
                        sizeof(int));
      output_file.write(reinterpret_cast<const char*>(&e_values[i]),
                        sizeof(float));
      output_file.write(reinterpret_cast<const char*>(&px_values[i]),
                        sizeof(float));
      output_file.write(reinterpret_cast<const char*>(&py_values[i]),
                        sizeof(float));
      output_file.write(reinterpret_cast<const char*>(&pz_values[i]),
                        sizeof(float));
    }
  }
  // Cleanup to be ready for the next event.
  particles.clear();
}

template <class T>
void JetScapeWriterFinalStateStream<T>::Init() {
  // Capitalize name for convenience
  std::string name = GetName();
  name[0] = toupper(name[0]);

  // Check if we want to write binary output
  binaryOutput = static_cast<bool>(JetScapeXML::Instance()->GetElementInt(
      {"Writer", "JetscapeWriterFinalStateBinaryOutput"}));
  // Whether to write the centrality and pt hat value for each event
  writeCentrality = static_cast<bool>(
      JetScapeXML::Instance()->GetElementInt({"write_centrality"}));
  // Whether to write the pt hat value for each event
  writePtHat = static_cast<bool>(
      JetScapeXML::Instance()->GetElementInt({"write_pthat"}));

  // Status codes to filter out from what is written (i.e. to be skipped)
  // NOTE: Need the capitalized name here!
  std::string s = JetScapeXML::Instance()->GetElementText(
      {"Writer", (std::string("FinalState") + name).c_str(), "statusToSkip"},
      false);
  if (s.size() > 0) {
    particleStatusToSkip = detail::stringToVector(s);
    // Remove the default value, if there
    particleStatusToSkip.erase(std::remove(particleStatusToSkip.begin(),
                                           particleStatusToSkip.end(), -9999),
                               particleStatusToSkip.end());
  }
  // Print the status codes that will be skipped for logging purposes to ensure
  // that it's clear that the values are propagated correctly.
  if (particleStatusToSkip.size() > 0) {
    std::stringstream ss;
    ss << "Filtering (i.e. removing) " << GetName() << " with status codes: ";
    unsigned int count = 0;
    for (const auto status : particleStatusToSkip) {
      if (count > 0) {
        ss << ", ";
      }
      ss << status;
      ++count;
    }
    JSINFO << ss.str();
  }

  if (GetActive()) {
    if (!binaryOutput) {
      // Writing ascii output
      // We need the header version to determine how to write.
      // NOTE: Don't require the version. If not specified, defaults to v2.
      int result = GetXMLElementInt(
          {"Writer", (std::string("FinalState") + name).c_str(),
           "headerVersion"},
          false);
      // If it fails to retrieve the value, it will return 0. The header version
      // must be >= 2, so only assign if the value is actually set.
      if (result) {
        headerVersion = static_cast<unsigned int>(result);
      }
      JSINFO << "JetScape Final State " << name << " Stream Writer v"
             << headerVersion
             << " initialized with output file = " << GetOutputFileName();
      // Normal ascii output
      output_file.open(GetOutputFileName().c_str());
      // NOTE: This header will only be printed once at the beginning on the
      // file.
      output_file << "#"
                  // The specifics the version number. For consistency in
                  // parsing, the string will always be "v<number>"
                  << "\t"
                  << "JETSCAPE_FINAL_STATE\t"
                  << "v" << headerVersion << "\t"
                  << "|"  // As a delimiter
                  << "\t"
                  << "N"
                  << "\t"
                  << "pid"
                  << "\t"
                  << "status"
                  << "\t"
                  << "E"
                  << "\t"
                  << "Px"
                  << "\t"
                  << "Py"
                  << "\t"
                  << "Pz"
                  << "\n";
    } else {
      // Binary output file
      output_file.open(GetOutputFileName().c_str(), std::ios::binary);
      const int binary_file_version = 3;
      output_file.write(reinterpret_cast<const char*>(&binary_file_version),
                        sizeof(int));
      // Status codes to skip in the writer, assume that there are at maximum 4
      // codes to skip. If there is no code to skip, write 666 instead. If there
      // are more than 4 codes to skip, write the first three and ignore the
      // rest (use JSWARN to message the user).
      if (particleStatusToSkip.size() > 0) {
        if (particleStatusToSkip.size() > 4) {
          JSWARN << "More than 4 status codes to skip. Only the first four "
                    "will be written to the file.";
        }
        int status_codes[4] = {666, 666, 666, 666};
        for (unsigned int i = 0; i < particleStatusToSkip.size() && i < 4;
             ++i) {
          status_codes[i] = particleStatusToSkip[i];
        }
        output_file.write(reinterpret_cast<const char*>(&status_codes),
                          sizeof(int) * 4);
      } else {
        int status_codes[4] = {666, 666, 666, 666};
        output_file.write(reinterpret_cast<const char*>(&status_codes),
                          sizeof(int) * 4);
      }
    }
  }
}

template <class T>
void JetScapeWriterFinalStateStream<T>::Exec() {
  // JSINFO<<"Run JetScapeWriterFinalStateStream<T>: Write event #
  // "<<GetCurrentEvent()<<" ...";

  // if (GetActive())
  //   WriteEvent();
}

template <class T>
void JetScapeWriterFinalStateStream<T>::Write(weak_ptr<PartonShower> ps) {
  auto pShower = ps.lock();
  if (!pShower)
    return;

  auto finalStatePartons = pShower->GetFinalPartons();

  // Store final state partons.
  for (const auto& parton : finalStatePartons) {
    particles.push_back(parton);
  }
}

template <class T>
void JetScapeWriterFinalStateStream<T>::Write(weak_ptr<Hadron> h) {
  auto hh = h.lock();
  if (hh) {
    particles.push_back(hh);
  }
}

template <class T>
void JetScapeWriterFinalStateStream<T>::Close() {
  if (!binaryOutput) {
    // Write xsec output at the end.
    // NOTE: Needs consistent "\t" between all entries to simplify parsing
    // later.
    output_file << "#"
                << "\t"
                << "sigmaGen\t" << GetHeader().GetSigmaGen() << "\t"
                << "sigmaErr\t" << GetHeader().GetSigmaErr() << "\n";
  } else {
    // Write xsec output at the end.
    const double sigmaGen = static_cast<double>(GetHeader().GetSigmaGen());
    output_file.write(reinterpret_cast<const char*>(&sigmaGen), sizeof(double));
    const double sigmaErr = static_cast<double>(GetHeader().GetSigmaErr());
    output_file.write(reinterpret_cast<const char*>(&sigmaErr), sizeof(double));
  }
  output_file.close();
}

template class JetScapeWriterFinalStatePartonsStream<ofstream>;
template class JetScapeWriterFinalStateHadronsStream<ofstream>;

#ifdef USE_GZIP
template class JetScapeWriterFinalStatePartonsStream<ogzstream>;
template class JetScapeWriterFinalStateHadronsStream<ogzstream>;
#endif

}  // end namespace Jetscape
