// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

// jetscape writer ascii class

#ifndef ROOT_JetScapeWriterRoot
#define ROOT_JetScapeWriterRoot

#include <fstream>
#include <string>

#include "JetScapeWriter.h"
#include "JetScapeRootEvent.h"

#include "TFile.h"
#include "TTree.h"

using namespace Jetscape;

class JetScapeWriterRoot : public JetScapeWriter
{

 public:

  JetScapeWriterRoot() {};
  JetScapeWriterRoot(string m_file_name_out);
  virtual ~JetScapeWriterRoot();

  void Init();
  void Exec();

  bool GetStatus() {return output_file->GetErrno();}
  void Close() {output_file->Write();output_file->Close();}

  void WriteInitFileXML();
  void Clear();

  void Write(weak_ptr<PartonShower> pS);

  //void Write(string s) {}
  //void WriteComment(string s) {}
  //void WriteWhiteSpace(string s) {}

  void WriteEvent();

 private:

  std::unique_ptr<TFile> output_file; //!< Output file
  std::shared_ptr<TTree> evTree;

  // ***********************************************************
  //std::unique_ptr<JetScapeRootEvent> event;
  // REMARK: Not really clear why smart pointers dont work here,
  // but in the other classes !???
  // What am I missing? do want to keep new/delete clean !!!!
  // To be followed up!!!

  JetScapeRootEvent *event;

  // ***********************************************************

  ClassDef(JetScapeWriterRoot,0)

};

#endif
