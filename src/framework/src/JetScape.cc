// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#include "JetScape.h"
#include "JetScapeXML.h"
#include "JetScapeSignalManager.h"
#include "JetEnergyLossManager.h"
#include "FluidDynamics.h"
#include "JetScapeBanner.h"
#include "InitialState.h"

#include<iostream>

using namespace std;

namespace Jetscape {

  /** Default constructor to create the main task of the JetScape framework. It sets the total number of events to 1.
   */
JetScape::JetScape()
{
  n_events=1;
  VERBOSE(8);
}

  /** This is a destructor for a JetScape.
   */
JetScape::~JetScape()
{
  VERBOSE(8);
  //JetScapeSignalManager::Instance()->Clear();
  //not needed, use weak_ptr in JetScapeSignalManager class (=not owning)
}

void JetScape::Show()
{
  show_jetscape_banner();
}

  /** This function initializes the main task of the JetScape framework. It calls JetScapeTask::InitTaks() function to initialize the modules/tasks of a JetScapeTask.
   */
void JetScape::Init()
{
  Show();
  
  INFO<<BOLDBLACK<<"Intialize JetScape ...";
  
  JetScapeXML::Instance()->OpenXMLFile(GetXMLFileName());
  
  string log_debug = JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement( "debug" )->GetText();
  string log_remark = JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement( "remark" )->GetText();
  int m_vlevel=0;
  
  JetScapeXML::Instance()->GetXMLRoot()->FirstChildElement("vlevel")->QueryIntText(&m_vlevel);
  
  // Some general JetScape interface settings can also be read in from XML file ...
  // Order of task could also be verifed via XML or via sort of TaskList with order parameter
  if ((int) log_debug.find("off")>=0)
    JetScapeLogger::Instance()->SetDebug(false);
  
  if ((int) log_remark.find("on")>=0)
    JetScapeLogger::Instance()->SetRemark(true);

  DEBUG<<"JetScape Debug from XML = "<< log_debug;
  DEBUG<<"JetScape Remark from XML = "<< log_remark;

   if ((int) m_vlevel>0)
     {
       JetScapeLogger::Instance()->SetVerboseLevel(m_vlevel);
       DEBUG<<"JetScape Verbose Level from XML = "<<m_vlevel;
     }

   SetPointers();
   
   // Set up helper. Mostly used for random numbers
   // Needs the XML reader singleton set up
   DEBUG<<"Seeding JetScapeTaskSupport from XML";
   JetScapeTaskSupport::ReadSeedFromXML( );

  // Has to be called explicitly since not really fully recursively (if ever needed)
   // So --> JetScape is "Task Manager" of all modules ...
   
   INFO<<"Found "<<GetNumberOfTasks()<<" Modules Initialize them ... ";
   JetScapeTask::InitTasks();
}

// kind of cluncky, maybe a better way ... ?
// Handle signal/slots in JetScape hence avoid passing pointers to sub tasks ...
void JetScape::SetPointers()
{
  
   // to get hydro pointer for signals, use signal?
  INFO<<"Set Hydro,JetEnergylossManager and IS Pointers for SignalManager to create Signal/Slots";
  
  for (auto it : GetTaskList())
    {
      if (dynamic_pointer_cast<InitialState>(it))
	JetScapeSignalManager::Instance()->SetInitialStatePointer(dynamic_pointer_cast<InitialState>(it));
 
      if (dynamic_pointer_cast<FluidDynamics>(it))
	JetScapeSignalManager::Instance()->SetHydroPointer(dynamic_pointer_cast<FluidDynamics>(it));
  
      if (dynamic_pointer_cast<JetEnergyLossManager>(it))
	JetScapeSignalManager::Instance()->SetJetEnergyLossManagerPointer(dynamic_pointer_cast<JetEnergyLossManager>(it));
      
      if (dynamic_pointer_cast<HardProcess>(it))
	JetScapeSignalManager::Instance()->SetHardProcessPointer(dynamic_pointer_cast<HardProcess>(it));
      
      if (dynamic_pointer_cast<JetScapeWriter>(it) && it->GetActive())
	JetScapeSignalManager::Instance()->SetWriterPointer(dynamic_pointer_cast<JetScapeWriter>(it)); 

      if (dynamic_pointer_cast<PartonPrinter>(it))
        JetScapeSignalManager::Instance()->SetPartonPrinterPointer(dynamic_pointer_cast<PartonPrinter>(it));
    }
}

  /** This function execute the modules/tasks of a JetScapeTask for all the events. It also calls "GetPartons()" function to print parton shower, and  "WriteTasks()" function to store the data in the XML file.
   */
void JetScape::Exec()
{
  INFO<<BOLDBLACK<<"Run JetScape ...";
  INFO<<BOLDBLACK<<"Number of Events = "<<GetNumberOfEvents(); 
  
  // JetScapeTask::ExecuteTasks(); Has to be called explicitly since not really fully recursively (if ever needed)
  // --> JetScape is "Task Manager" of all modules ...

  // Simple way of passing the writer module pointer ...
  weak_ptr<JetScapeWriter> w;
  //weak_ptr<PartonPrinter> p;  

  for (auto it : GetTaskList())
  {
    if (dynamic_pointer_cast<JetScapeWriter>(it))
    {  
      if (it->GetActive())
        w=dynamic_pointer_cast<JetScapeWriter>(it);	           
    }
    //else if(dynamic_pointer_cast<PartonPrinter>(it))
    //{
        //p=dynamic_pointer_cast<PartonPrinter>(it);
    //}
  } 
 
  for (int i=0;i<GetNumberOfEvents();i++)
    {
      INFO<<BOLDBLACK<<"Run Event # = "<<i;
      DEBUG<<"Found "<<GetNumberOfTasks()<<" Modules Execute them ... ";
      
      JetScapeTask::ExecuteTasks();

      //JetScapeTask::GetPartons(p);

      if (w.lock().get())
	JetScapeTask::WriteTasks(w);            

      JetScapeTask::ClearTasks();

      IncrementCurrentEvent();
    }
}

  /**
   */
void JetScape::Finish()
{
  INFO<<BOLDBLACK<<"JetScape finished after "<<GetNumberOfEvents()<<" events!";
  DEBUG<<"More infos wrap up/saving to file/closing file ...";

  // same as in Init() and Exec() ...
  JetScapeTask::FinishTasks(); //dummy so far ...
}

} // end namespace Jetscape
