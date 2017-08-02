// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#include "JetScapeTaskSupport.h"
#include "JetScapeLogger.h"
#include <stdlib.h>

using namespace std;

namespace Jetscape {


  JetScapeTaskSupport* JetScapeTaskSupport::m_pInstance = nullptr;

  JetScapeTaskSupport* JetScapeTaskSupport::Instance()
  {
    if (!m_pInstance) {
      m_pInstance = new JetScapeTaskSupport;
      INFO<<"Created JetScapeTaskSupport Instance";
    }
    
    return m_pInstance;
  }

			   
  int JetScapeTaskSupport::RegisterTask(){
    INFO << "JetScapeTaskSupport::RegisterTask called, answering " << CurrentTaskNumber;
    CurrentTaskNumber++;
    return CurrentTaskNumber-1;
  }
  
  // void JetScapeTaskSupport::CleanUp()
  // {
  //   VERBOSE(8);
  //   VERBOSE(8)<<"Done ...";
  // }

} // end namespace Jetscape
