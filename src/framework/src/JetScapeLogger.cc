// JetScape Logger class implementation (meant as singelton)

#include<stddef.h>
#include<fstream>

#include "JetScapeLogger.h"

using namespace std;

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

std::ostringstream null;     

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

LogStreamer JetScapeLogger::Debug()
{
  if (debug)
    {
      string s="[Debug] ";
      //s <<  __PRETTY_FUNCTION__ <<":"<<__LINE__<<" ";
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
  //return LogStreamer(std::cout<<BOLDBLUE<<s);
  return LogStreamer(std::cout<<s);
}


LogStreamer JetScapeLogger::Remark()
{
  if (remark)
    {
      string s="[Remark] ";
      return LogStreamer(std::cout<<CYAN<<s);
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
      //string s="[Verbose] [level = "; s+= std::to_string(vlevel); s+="] ";
      // shorter:
      string s="[Verbose][";s+= std::to_string(m_vlevel); s+="] ";
      return LogStreamer(std::cout<<GREEN<<s);
    }
  else
    {
      null.setstate(std::ios_base::failbit);
      return LogStreamer(null);  
    }
}
