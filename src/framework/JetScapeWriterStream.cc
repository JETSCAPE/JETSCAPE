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
// jetscape writer ascii class

#include "JetScapeWriterStream.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"

namespace Jetscape {

// Register the modules with the base class
template <>
RegisterJetScapeModule<JetScapeWriterStream<ofstream>>
    JetScapeWriterStream<ofstream>::reg("JetScapeWriterAscii");
template <>
RegisterJetScapeModule<JetScapeWriterStream<ogzstream>>
    JetScapeWriterStream<ogzstream>::regGZ("JetScapeWriterAsciiGZ");

template <class T>
JetScapeWriterStream<T>::JetScapeWriterStream(string m_file_name_out) {
  SetOutputFileName(m_file_name_out);
}

template <class T> JetScapeWriterStream<T>::~JetScapeWriterStream() {
  VERBOSE(8);
  if (GetActive())
    Close();
}

template <class T> void JetScapeWriterStream<T>::WriteHeaderToFile() {
  VERBOSE(3) << "Run JetScapeWriterStream<T>: Write header of event # "
             << GetCurrentEvent() << " ...";
  Write(to_string(GetCurrentEvent()) + " Event");

  std::ostringstream oss;
  oss.str("");
  oss << GetId() << "sigmaGen " << GetHeader().GetSigmaGen();
  WriteComment(oss.str());
  oss.str("");
  oss << GetId() << "sigmaErr " << GetHeader().GetSigmaErr();
  WriteComment(oss.str());
  oss.str("");
  oss << GetId() << "weight " << GetHeader().GetEventWeight();
  WriteComment(oss.str());

  if (GetHeader().GetNpart() > -1) {
    oss.str("");
    oss << GetId() << "Npart " << GetHeader().GetNpart();
    WriteComment(oss.str());
  }
  if (GetHeader().GetNcoll() > -1) {
    oss.str("");
    oss << GetId() << "Ncoll " << GetHeader().GetNcoll();
    WriteComment(oss.str());
  }
  if (GetHeader().GetTotalEntropy() > -1) {
    oss.str("");
    oss << GetId() << "TotalEntropy " << GetHeader().GetTotalEntropy();
    WriteComment(oss.str());
  }

  if (GetHeader().GetEventPlaneAngle() > -999) {
    oss.str("");
    oss << GetId() << "EventPlaneAngle " << GetHeader().GetEventPlaneAngle();
    WriteComment(oss.str());
  }
}

template <class T> void JetScapeWriterStream<T>::WriteEvent() {
  // JSINFO<<"Run JetScapeWriterStream<T>: Write event # "<<GetCurrentEvent()<<" ...";
  // do nothing, the modules handle this
}

template <class T> void JetScapeWriterStream<T>::Write(weak_ptr<Parton> p) {
  auto pp = p.lock();
  if (pp) {
    output_file << *pp << endl;
  }
}

template <class T> void JetScapeWriterStream<T>::Write(weak_ptr<Vertex> v) {
  auto vv = v.lock();
  if (vv) {
    output_file << *vv << endl;
  }
}

template <class T> void JetScapeWriterStream<T>::Init() {
  if (GetActive()) {
    JSINFO << "JetScape Stream Writer initialized with output file = "
           << GetOutputFileName();
    output_file.open(GetOutputFileName().c_str());

    //Write Init Informations, like XML and ... to file ...
    //WriteInitFileXMLMain();
    //WriteInitFileXMLUser();
  }
}

template <class T> void JetScapeWriterStream<T>::Exec() {
  // JSINFO<<"Run JetScapeWriterStream<T>: Write event # "<<GetCurrentEvent()<<" ...";

  // if (GetActive())
  //   WriteEvent();
}

template <class T> void JetScapeWriterStream<T>::WriteInitFileXMLMain() {
  JSDEBUG << "Write XML Main to output file. XML file = "
          << JetScapeXML::Instance()->GetXMLMainFileName();
  tinyxml2::XMLPrinter printer;
  JetScapeXML::Instance()->GetXMLDocumentMain().Print(&printer);
  WriteComment("Init XML Main file used : " +
               JetScapeXML::Instance()->GetXMLMainFileName());
  output_file << printer.CStr();
}

template <class T> void JetScapeWriterStream<T>::WriteInitFileXMLUser() {
  JSDEBUG << "Write XML User to output file. XML file = "
          << JetScapeXML::Instance()->GetXMLUserFileName();
  tinyxml2::XMLPrinter printer;
  JetScapeXML::Instance()->GetXMLDocumentUser().Print(&printer);
  WriteComment("Init XML User file used : " +
               JetScapeXML::Instance()->GetXMLUserFileName());
  output_file << printer.CStr();
}

template <class T>
void JetScapeWriterStream<T>::Write(weak_ptr<PartonShower> ps) {
  auto pShower = ps.lock();
  if (!pShower)
    return;

  WriteComment(
      "Parton Shower in JetScape format to be used later by GTL graph:");

  // write vertices
  PartonShower::node_iterator nIt, nEnd;

  for (nIt = pShower->nodes_begin(), nEnd = pShower->nodes_end(); nIt != nEnd;
       ++nIt) {
    WriteWhiteSpace("[" + to_string(nIt->id()) + "] V");
    Write(pShower->GetVertex(*nIt));
  }

  PartonShower::edge_iterator eIt, eEnd;
  for (eIt = pShower->edges_begin(), eEnd = pShower->edges_end(); eIt != eEnd;
       ++eIt) {
    WriteWhiteSpace("[" + to_string(eIt->source().id()) + "]=>[" +
                    to_string(eIt->target().id()) + "] P");
    Write(pShower->GetParton(*eIt));
  }
}

template <class T> void JetScapeWriterStream<T>::Write(weak_ptr<Hadron> h) {
  auto hh = h.lock();
  if (hh) {
    output_file << *hh << endl;
  }
}


template <class T>
void JetScapeWriterStream<T>::Write(std::weak_ptr<Qvector> Qv) {
    auto Qv_tem = Qv.lock();
    if (Qv_tem) {
        for (int ipt = 0; ipt < Qv_tem->get_npt(); ipt++) {
            for (int iy = 0; iy < Qv_tem->get_ny(); iy++) {
                output_file << Qv_tem->get_header()<<" "
                            << Qv_tem->get_pdgcode()<<" "
			    << Qv_tem->get_pt(ipt)<<" "
			    << Qv_tem->get_y(iy)<<" "
		            << Qv_tem->get_value(ipt, iy, 1) <<" "
			    <<Qv_tem->get_value(ipt, iy, 3)<<" ";
                for (int icol = 4; icol < Qv_tem->get_ncols(); icol++) {
                    output_file << Qv_tem->get_value(ipt, iy, icol)<<" ";
                }
                output_file << std::endl;
            }
            //output_file << std::endl;
        }
    }
}


template class JetScapeWriterStream<ofstream>;

#ifdef USE_GZIP
template class JetScapeWriterStream<ogzstream>;
#endif

} // end namespace Jetscape
