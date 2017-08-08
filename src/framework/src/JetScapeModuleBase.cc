// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#include "JetScapeModuleBase.h"
#include "JetScapeXML.h"
#include "JetScapeTaskSupport.h"
#include "JetScapeLogger.h"

#include<iostream>

namespace Jetscape {

  int JetScapeModuleBase::current_event=0;

  // ---------------------------------------------------------------------------
  JetScapeModuleBase::JetScapeModuleBase()
  {
    xml_file_name = "";
  }

  // ---------------------------------------------------------------------------
  JetScapeModuleBase::JetScapeModuleBase(string m_name)
  {
    xml_file_name = m_name;
  }

  // ---------------------------------------------------------------------------
  JetScapeModuleBase::~JetScapeModuleBase()
  {
    disconnect_all();
  }

  // ---------------------------------------------------------------------------
  void JetScapeModuleBase::Init()
  {
    if (!JetScapeXML::Instance()->GetXMLRoot())
      {
	WARN << "Not a valide JetScape XML file or no XML file loaded!";
	exit(-1);
      }
  }

  // ---------------------------------------------------------------------------
  shared_ptr<std::mt19937> JetScapeModuleBase::get_mt19937_generator(){
    // Instantiate if it isn't there yet
    if ( !mt19937_generator_ ){
      mt19937_generator_ = JetScapeTaskSupport::Instance()->get_mt19937_generator( get_my_task_number() );
    }
    return mt19937_generator_;
  }
  

} // end namespace Jetscape
