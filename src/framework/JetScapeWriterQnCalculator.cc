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
// Jetscape Qn vector writer ascii class
// Based on JetScapeWriterStream, JetScapeWriterFinalStateStream
// author: Xiang-Yu Wu <xiangyu.wu2@mail.mcgill.ca>

#include "JetScapeWriterQnCalculator.h"
#include "JetScapeLogger.h"
#include "JetScapeXML.h"

namespace Jetscape {

// Register the modules with the base class
template <>
RegisterJetScapeModule<JetScapeWriterQnVectorStream<ofstream>>
    JetScapeWriterQnVectorStream<ofstream>::regQnVector("JetScapeWriterQnVectorAscii");
template <>
RegisterJetScapeModule<JetScapeWriterQnVectorStream<ogzstream>>
    JetScapeWriterQnVectorStream<ogzstream>::regQnVectorGZ("JetScapeWriterQnVectorAsciiGZ");

template <class T>
JetScapeWriterQnVectorStream<T>::JetScapeWriterQnVectorStream(string m_file_name_out):
  writeCentrality{false}
{
  SetOutputFileName(m_file_name_out);
}

template <class T> JetScapeWriterQnVectorStream<T>::~JetScapeWriterQnVectorStream() {
  VERBOSE(8);
  if (GetActive())
    Close();
}

template <class T> void JetScapeWriterQnVectorStream<T>::WriteEvent() {
  // Write the entire event all at once.

  // Optionally write event centrality to event header
  std::string centrality_text = "";
  if (writeCentrality) {
    centrality_text += "\tcentrality\t";
    centrality_text += std::to_string(GetHeader().GetEventCentrality());
  }

  // First, write header
  // NOTE: Needs consistent "\t" between all entries to simplify parsing later.
  // NOTE: Could also add Npart, Ncoll, and TotalEntropy. See the original stream writer.
  output_file << "#"
      << "\t" << "Event\t" << GetCurrentEvent() + 1;

  output_file
      << centrality_text
      <<  "\n";

  double oversamplenevent = JetScapeXML::Instance()->GetElementDouble({"SoftParticlization", "iSS", "number_of_repeated_sampling"});

  std::vector<int> pdglist={211,2212,321,-211,-2212,-321,9999};

  for(int id = 0; id < pdglist.size(); id++)
  {  
    int select_pid = pdglist[id];
    int rapidity_kind = 1;
    if(select_pid == 9999){
      rapidity_kind = 0;
    }

    std::shared_ptr<Qvector> Qvec_tem_ptr = std::make_shared<Qvector>(pTmin_,pTmax_,npT_,rapmin_,rapmax_,nrap_,norder_,select_pid,rapidity_kind);
    int num_pid = 0;
    for (const auto particle : particles) {

      if( particle->pstat()!= 11 && particle->pstat()!=27) continue;
      // only includes hadrons from iSS and SMASH
      PdgCode pdgcode_tem(std::to_string(particle->pid()));
              
      if (select_pid == 9999) {
        if (std::abs(pdgcode_tem.charge()) > 0) {
          Qvec_tem_ptr->fill_particle(particle);
          num_pid = num_pid + 1;
        }
      } 
      else if (pdgcode_tem.get_decimal() == select_pid) {
        Qvec_tem_ptr->fill_particle(particle);
        num_pid = num_pid + 1;

      }
    }

    double dpt = Qvec_tem_ptr->get_dpt();
    double dy = Qvec_tem_ptr->get_dy();

    for (int ipt = 0; ipt < Qvec_tem_ptr->get_npt(); ipt++) {
      for (int iy = 0; iy < Qvec_tem_ptr->get_ny(); iy++) {
          output_file << Qvec_tem_ptr->get_pdgcode()<<" ";
          double total_dN = Qvec_tem_ptr->get_value(ipt, iy, 5);
          double total_dN_mean = Qvec_tem_ptr->get_value(ipt, iy, 5)/oversamplenevent;
          double total_dN_mean_err = sqrt(total_dN_mean/oversamplenevent);
          
          double total_ET_mean ;
          total_ET_mean = Qvec_tem_ptr->get_value(ipt, iy, 4)/oversamplenevent;
          
          double mean_pT, mean_pT_err;
          double mean_y, mean_y_err;
          if(total_dN_mean>0.0){
              mean_pT = Qvec_tem_ptr->get_value(ipt,iy,0)/total_dN;
              mean_pT_err = sqrt(Qvec_tem_ptr->get_value(ipt,iy,1)/total_dN - mean_pT*mean_pT)/sqrt(total_dN);
              mean_y = Qvec_tem_ptr->get_value(ipt,iy,2)/total_dN;
              mean_y_err = sqrt(Qvec_tem_ptr->get_value(ipt,iy,3)/total_dN - mean_y*mean_y)/sqrt(total_dN);
              
          }
          else{
              mean_pT = Qvec_tem_ptr->get_pt(ipt);
              mean_pT_err = 0.0;
              mean_y = Qvec_tem_ptr->get_y(iy);
              mean_y_err = 0.0;
          }

          //mean_pT = Qvec_tem_ptr->get_pt(ipt);
          //mean_pT_err = 0.0;
          //mean_y = Qvec_tem_ptr->get_y(iy);
          //mean_y_err = 0.0;
          
          output_file<<mean_pT<<" "<<mean_pT_err<<" "
                     <<mean_y <<" "<<mean_y_err<<" "
                     <<total_ET_mean<<" "
                     <<total_dN_mean/(dy*dpt)<<" "<<total_dN_mean_err/(dy*dpt)<<" ";
          for(int iorder = 1; iorder < norder_;iorder++){
             double Qn_real_mean = 0.0;
             double Qn_imag_mean = 0.0;
             double Qn_real_err = 0.0;
             double Qn_imag_err = 0.0;
             if(total_dN_mean>0.0){
                 Qn_real_mean = Qvec_tem_ptr->get_value(ipt,iy,4*iorder+3)/total_dN; 
                 Qn_imag_mean = Qvec_tem_ptr->get_value(ipt,iy,4*iorder+4)/total_dN;
                 Qn_real_err = sqrt(Qvec_tem_ptr->get_value(ipt,iy,4*iorder+5)/total_dN - Qn_real_mean*Qn_real_mean)/sqrt(total_dN); 
                 if (std::isnan(Qn_real_err)) Qn_real_err = 0.0;
                 Qn_imag_err = sqrt(Qvec_tem_ptr->get_value(ipt,iy,4*iorder+6)/total_dN - Qn_imag_mean*Qn_imag_mean)/sqrt(total_dN);  
                 if (std::isnan(Qn_imag_err)) Qn_imag_err = 0.0;
             }
             output_file<<Qn_real_mean<<" "<<Qn_real_err<<" "
                        <<Qn_imag_mean<<" "<<Qn_imag_err<<" ";
             
          }
          output_file<<total_dN<<std::endl;    

      }
    }
  }

  // Cleanup to be ready for the next event.
  particles.clear();
}

template <class T> void JetScapeWriterQnVectorStream<T>::Init() {
  // Whether to write the centrality for each event
  writeCentrality = static_cast<bool>(JetScapeXML::Instance()->GetElementInt({"write_centrality"}));

  if (GetActive()) {
    // Capitalize name
    std::string name = GetName();
    name[0] = toupper(name[0]);
    JSINFO << "JetScape Final State " << name << " Stream Writer initialized with output file = "
           << GetOutputFileName();

    pTmin_ = JetScapeXML::Instance()->GetElementDouble({"QnVector_pTmin"});
    pTmax_ = JetScapeXML::Instance()->GetElementDouble({"QnVector_pTmax"});
    npT_ = JetScapeXML::Instance()->GetElementInt({"QnVector_NpT"});
    
    rapmin_ = JetScapeXML::Instance()->GetElementDouble({"QnVector_rapmin"});
    rapmax_ = JetScapeXML::Instance()->GetElementDouble({"QnVector_rapmax"});
    nrap_ = JetScapeXML::Instance()->GetElementInt({"QnVector_Nrap"});
    norder_ = JetScapeXML::Instance()->GetElementInt({"QnVector_Norder"});
    
    output_file.open(GetOutputFileName().c_str());
    // NOTE: This header will only be printed once at the beginning on the file.
    output_file << "#"
        // The specifics the version number. For consistency in parsing, the string
        // will always be "v<number>"
        << "\t" << "JETSCAPE_FINAL_STATE\t" << "QnVector"
        << "\t" << "\n"
        << "#"
        << "\t" << "pTmin\t" << pTmin_ << "\t" << "pTmax\t" << pTmax_ << "\t" << "NpT\t" << npT_ 
        << "\t" << "rapmin\t" << rapmin_ << "\t" << "rapmax\t" << rapmax_ << "\t" << "Nrap\t" << nrap_ 
        << "\t" << "Norder\t" << norder_ << "\n"
        << "#"
        << "\t" << "PID of Charged particle 9999 \t" << " y(Charged) = pseudo-rapidity \t y(PID) = rapidity"
        << "\t" << "\n"  
        << "#"
        << "\t" << "pid"
        << "\t" << "pT"
        << "\t" << "pT_err"
        << "\t" << "y"
        << "\t" << "y_err"
        << "\t" << "ET"
        << "\t" << "dNdpTdy"
        << "\t" << "dNdpTdy_err"
        << "\t" << "vncos"
        << "\t" << "vncos_err"
        << "\t" << "vnsin"
        << "\t" << "vnsin_err"
        << "\t" << "dN"
        << "\n";

  }
}

template <class T> void JetScapeWriterQnVectorStream<T>::Exec() {
  // JSINFO<<"Run JetScapeWriterFinalStateStream<T>: Write event # "<<GetCurrentEvent()<<" ...";

  // if (GetActive())
  //   WriteEvent();
}

template <class T> void JetScapeWriterQnVectorStream<T>::Write(weak_ptr<Hadron> h) {
  auto hh = h.lock();
  if (hh) {
    particles.push_back(hh);
  }
}

template <class T> void JetScapeWriterQnVectorStream<T>::Close() {
    // Write xsec output at the end.
    // NOTE: Needs consistent "\t" between all entries to simplify parsing later.
  output_file << "#"
      << "\t" << "Event\t" << GetCurrentEvent()   <<  " End \n";
    output_file.close();
}

template class JetScapeWriterQnVectorStream<ofstream>;
#ifdef USE_GZIP
template class JetScapeWriterQnVectorStream<ogzstream>;
#endif

} // end namespace Jetscape