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
// Jetscape final state {hadrons,kartons} writer ascii class
// Based on JetScapeWriterStream.
// author: Raymond Ehlers <raymond.ehlers@cern.ch>, ORNL

#include "JetScapeWriterFinalStateStream.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"

namespace Jetscape {

// Register the modules with the base class
template <>
RegisterJetScapeModule<JetScapeWriterFinalStatePartonsStream<ofstream>>
    JetScapeWriterFinalStatePartonsStream<ofstream>::regParton("JetScapeWriterFinalStatePartonsAscii");
template <>
RegisterJetScapeModule<JetScapeWriterFinalStateHadronsStream<ofstream>>
    JetScapeWriterFinalStateHadronsStream<ofstream>::regHadron("JetScapeWriterFinalStateHadronsAscii");
template <>
RegisterJetScapeModule<JetScapeWriterFinalStatePartonsStream<ogzstream>>
    JetScapeWriterFinalStatePartonsStream<ogzstream>::regPartonGZ("JetScapeWriterFinalStatePartonsAsciiGZ");
template <>
RegisterJetScapeModule<JetScapeWriterFinalStateHadronsStream<ogzstream>>
    JetScapeWriterFinalStateHadronsStream<ogzstream>::regHadronGZ("JetScapeWriterFinalStateHadronsAsciiGZ");

template <class T>
JetScapeWriterFinalStateStream<T>::JetScapeWriterFinalStateStream(string m_file_name_out) {
  SetOutputFileName(m_file_name_out);
}

template <class T> JetScapeWriterFinalStateStream<T>::~JetScapeWriterFinalStateStream() {
  VERBOSE(8);
  if (GetActive())
    Close();
}

template <class T> void JetScapeWriterFinalStateStream<T>::WriteEvent() {
  // Write the entire event all at once.

  // Optionally write pt-hat value to event header
  std::string pt_hat_text = "";
  int write_pthat = JetScapeXML::Instance()->GetElementInt({"write_pthat"});
  if (write_pthat) {
    pt_hat_text += "\tpt_hat\t";
    pt_hat_text += std::to_string(GetHeader().GetPtHat());
  }

  // First, write header
  // NOTE: Needs consistent "\t" between all entries to simplify parsing later.
  // NOTE: Could also add Npart, Ncoll, and TotalEntropy. See the original stream writer.
  output_file << "#"
      << "\t" << "Event\t" << GetCurrentEvent() + 1  // +1 to index the event count from 1
      << "\t" << "weight\t" << GetHeader().GetEventWeight()
      << "\t" << "EPangle\t" << (GetHeader().GetEventPlaneAngle() > -999 ? GetHeader().GetEventPlaneAngle() : 0)
      << "\t" << "N_" << GetName() << "\t" << particles.size()
      << pt_hat_text
      <<  "\n";

  // Next, write the particles. Will contain either hadrons or partons based on the derived class.
  unsigned int ipart = 0;
  for (const auto & p : particles) {
    auto particle = p.get();
    output_file << ipart
        << " " << particle->pid()
        << " " << particle->pstat()
        << " " << particle->e()
        << " " << particle->px()
        << " " << particle->py()
        << " " << particle->pz()
        << "\n";
    ++ipart;
  }

  // Cleanup to be ready for the next event.
  particles.clear();
}

template <class T> void JetScapeWriterFinalStateStream<T>::Init() {
  if (GetActive()) {
    // Capitalize name
    std::string name = GetName();
    name[0] = toupper(name[0]);
    JSINFO << "JetScape Final State " << name << " Stream Writer initialized with output file = "
           << GetOutputFileName();
    output_file.open(GetOutputFileName().c_str());
    // NOTE: This header will only be printed once at the beginning on the file.
    output_file << "#"
        // The specifics the version number. For consistency in parsing, the string
        // will always be "v<number>"
        << "\t" << "JETSCAPE_FINAL_STATE\t" << "v2"
        << "\t" << "|"  // As a delimiter
        << "\t" << "N"
        << "\t" << "pid"
        << "\t" << "status"
        << "\t" << "E"
        << "\t" << "Px"
        << "\t" << "Py"
        << "\t" << "Pz"
        << "\n";
  }
}

template <class T> void JetScapeWriterFinalStateStream<T>::Exec() {
  // JSINFO<<"Run JetScapeWriterFinalStateStream<T>: Write event # "<<GetCurrentEvent()<<" ...";

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
  for (const auto parton : finalStatePartons) {
      particles.push_back(parton);
  }
}

template <class T> void JetScapeWriterFinalStateStream<T>::Write(weak_ptr<Hadron> h) {
  auto hh = h.lock();
  if (hh) {
    particles.push_back(hh);
  }
}

template <class T> void JetScapeWriterFinalStateStream<T>::Close() {
    // Write xsec output at the end.
    // NOTE: Needs consistent "\t" between all entries to simplify parsing later.
    output_file << "#" << "\t"
        << "sigmaGen\t" << GetHeader().GetSigmaGen() << "\t"
        << "sigmaErr\t" << GetHeader().GetSigmaErr() << "\n";
    output_file.close();
}

template class JetScapeWriterFinalStatePartonsStream<ofstream>;
template class JetScapeWriterFinalStateHadronsStream<ofstream>;

#ifdef USE_GZIP
template class JetScapeWriterFinalStatePartonsStream<ogzstream>;
template class JetScapeWriterFinalStateHadronsStream<ogzstream>;
#endif

} // end namespace Jetscape
