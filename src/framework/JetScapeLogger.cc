/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * Modular, task-based framework
 * Intial Design: Joern Putschke, Kolja Kauder (Wayne State University)
 * For the full list of contributors see AUTHORS.

 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 * or via email to bugs.jetscape.org@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#include<stddef.h>
#include<fstream>

#include "JetScapeLogger.h"

using namespace std;


// For Linux systems (also Mac)
// only shows max usage throughtout the program

#include <sys/time.h>
#include <sys/resource.h>

long getMemoryUsage() 
{
  struct rusage usage;
  // NOTE: Reported in kB on BSD/Linux, bytes in Mac/Darwin
  // Could try to explicitly catch __linux__ as well
  float mbsize = 1024;
#ifdef __MACH__
  mbsize = 1024 * 1024;
#endif
    
  if(0 == getrusage(RUSAGE_SELF, &usage))
    return usage.ru_maxrss/mbsize;
  else
    return 0;
}

 /*
//For mac os x only ...
// shows current usage
#include <mach/mach.h>
#include <mach/task.h>

int getMemoryUsage()
{
  struct mach_task_basic_info info;
  mach_msg_type_number_t size = MACH_TASK_BASIC_INFO_COUNT;
  kern_return_t kerr = task_info(mach_task_self(),
				 MACH_TASK_BASIC_INFO,
				 (task_info_t)&info,
				 &size);
  if( kerr == KERN_SUCCESS )
    return info.resident_size/1024./1024.;
  else {
    cout<<"Error with task_info(): "<<mach_error_string(kerr)<<endl;
    return -1;}
}
 */
 
// Just some definition of colors for colored terminal log output
// Resetting to default color in LogStreamer!

#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m" //\033[7;30m" bkg     /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */
 
#define CLEAR "\033[2J"  // clear screen escape code 

namespace Jetscape {

std::ostringstream null;     

safe_ostream safe_cout(std::cout);
safe_ostream safe_cerr(std::cerr);
safe_ostream safe_null(null);

JetScapeLogger* JetScapeLogger::m_pInstance = NULL;

JetScapeLogger* JetScapeLogger::Instance()
{
  if (!m_pInstance)
    m_pInstance = new JetScapeLogger();
  
  return m_pInstance;
}

LogStreamer JetScapeLogger::Warn()
{
  string s="[Warning] ";
  //s << __PRETTY_FUNCTION__ <<":"<<__LINE__<<" ";
  //return LogStreamer(std::cout<<s<<__PRETTY_FUNCTION__ <<":"<<__LINE__<<" ");
  return LogStreamer(std::cout<<BOLDRED<<s);
}

LogStreamerThread JetScapeLogger::DebugThread()
{
  if (debug)
    {
      string s="[Debug Thread] ";
      //s <<  __PRETTY_FUNCTION__ <<":"<<__LINE__<<" ";
      s += to_string(getMemoryUsage()); s+="MB ";
      return (LogStreamerThread(safe_cout) << BLUE << s);     
    }
  else
    {
      // check if it is not written in some system log files ...
      //safe_null.setstate(std::ios_base::failbit);
      return LogStreamerThread(safe_null);  
    }
}

LogStreamer JetScapeLogger::Debug()
{
  if (debug)
    {
      string s="[Debug] ";
      //s <<  __PRETTY_FUNCTION__ <<":"<<__LINE__<<" ";
      s += to_string(getMemoryUsage()); s+="MB ";
      return LogStreamer(std::cout<<BLUE<<s);
    }
  else
    {
      null.setstate(std::ios_base::failbit);
      return LogStreamer(null);  
    }
}

LogStreamer JetScapeLogger::Info()
{
  string s="[Info] ";
  // s <<  __PRETTY_FUNCTION__ <<":"<<__LINE__<<" ";
  if (info){
    string s="[Info] ";
    // s <<  __PRETTY_FUNCTION__ <<":"<<__LINE__<<" ";
    s += to_string(getMemoryUsage()); s+="MB ";
    return LogStreamer(std::cout<<s);
  }  else {
    null.setstate(std::ios_base::failbit);
    return LogStreamer(null);
  }

}

LogStreamer JetScapeLogger::InfoNice()
{
  string s="[Info] ";
  // s <<  __PRETTY_FUNCTION__ <<":"<<__LINE__<<" ";
  //s += to_string(getMemoryUsage()); s+="MB ";
 if (info){
   return LogStreamer(std::cout<<s);
 }  else {
   null.setstate(std::ios_base::failbit);
   return LogStreamer(null);
 }
}


LogStreamer JetScapeLogger::Remark()
{
  if (remark)
    {
      string s="[REMARK] ";
      return LogStreamer(std::cout<<BOLDMAGENTA<<s);
    }
  else
    {
      null.setstate(std::ios_base::failbit);
      return LogStreamer(null);  
    }
}

LogStreamer JetScapeLogger::Verbose(unsigned short m_vlevel)
{
  if (m_vlevel<vlevel) // or if (m_vlevel==vlevel)
    {
      string s="[Verbose][";s+= std::to_string(m_vlevel); s+="] ";
      s += to_string(getMemoryUsage()); s+="MB ";
      return LogStreamer(std::cout<<GREEN<<s);
    }
  else
    {
      null.setstate(std::ios_base::failbit);
      return LogStreamer(null);  
    }
}


LogStreamer JetScapeLogger::VerboseShower(unsigned short m_vlevel)
{
  if (m_vlevel<vlevel) // or if (m_vlevel==vlevel)
    {
      string s="[Verbose][";s+= std::to_string(m_vlevel); s+="] ";
      s += to_string(getMemoryUsage()); s+="MB ";
      return LogStreamer(std::cout<<BOLDCYAN<<s);
    }
  else
    {
      null.setstate(std::ios_base::failbit);
      return LogStreamer(null);  
    }
}


LogStreamer JetScapeLogger::VerboseParton(unsigned short m_vlevel,Parton &p)
{
  if (m_vlevel<vlevel) // or if (m_vlevel==vlevel)
    {
      string s="[Verbose][";s+= std::to_string(m_vlevel); s+="] Parton: ";
      return LogStreamer(std::cout<<GREEN<<s<<" "<<p<<endl);
    }    
  else
    {
      null.setstate(std::ios_base::failbit);
      return LogStreamer(null);  
    }
}

LogStreamer JetScapeLogger::VerboseVertex(unsigned short m_vlevel,Vertex &v)
{
  if (m_vlevel<vlevel) // or if (m_vlevel==vlevel)
    {
      string s="[Verbose][";s+= std::to_string(m_vlevel); s+="] Vertex: ";
      return LogStreamer(std::cout<<GREEN<<s<<" "<<v<<endl);
    }    
  else
    {
      null.setstate(std::ios_base::failbit);
      return LogStreamer(null);  
    }
}

} // end namespace Jetscape
