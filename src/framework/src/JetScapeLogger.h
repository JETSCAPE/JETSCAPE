// JetScape Logger class (meant as singelton)

#ifndef JETSCAPELOGGER_H
#define JETSCAPELOGGER_H

#include<iostream>
#include<sstream>

using namespace std;

#define RESET   "\033[0m"
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */

// define nicer macros to be used for logging ...
#define INFO  JetScapeLogger::Instance()->Info()<<" " //<<__PRETTY_FUNCTION__<<" : "
#define INFO_NICE  JetScapeLogger::Instance()->Info()
#define DEBUG JetScapeLogger::Instance()->Debug()<<__PRETTY_FUNCTION__<<" : "
#define REMARK JetScapeLogger::Instance()->Remark()<<__PRETTY_FUNCTION__<<" : "
#define VERBOSE(l) JetScapeLogger::Instance()->Verbose(l)<<__PRETTY_FUNCTION__<<" : "
#define WARN JetScapeLogger::Instance()->Warn()<<__PRETTY_FUNCTION__<<" : "

// --------------------------------
// Just a helper class to make the interface
// consistent with << operator
// In principle simple extension to
// log into a file ... via m_dest
// Think about thread safety ...
// << overload in Parton class not working!? Check ...


class LogStreamer
{
  
  shared_ptr< std::ostringstream > m_collector;
  std::ostream* m_dest;

 public:
  
    LogStreamer( std::ostream& dest )
    {
      m_collector=make_shared<std::ostringstream>();
      m_dest = &dest;
    };
    
    ~LogStreamer()
    {    
      if ( m_collector.unique() && m_dest!=nullptr) {	
	*m_dest << m_collector->str() << RESET << std::endl;
      }
    }
    
    template <typename T>
    LogStreamer& operator<<( T const& value )
    {
        *m_collector << value;
        return *this;
    }
};

// --------------------------------

class JetScapeLogger
{
  
 public:

  static JetScapeLogger* Instance();

  LogStreamer Info();
  LogStreamer Warn();
  LogStreamer Debug();
  LogStreamer Remark();
  //LogStreamer Error(); //to be implemented
  //LogStreamer Fatal(); //to be implemented
  
  LogStreamer Verbose(unsigned short m_vlevel);
  
  void SetDebug(bool m_debug) {debug=m_debug;}
  void SetRemark(bool m_remark) {remark=m_remark;}
  void SetVerboseLevel(unsigned short m_vlevel) {vlevel=m_vlevel;}
  bool GetDebug() {return debug;}  
  bool GetRemark() {return remark;}
  unsigned short GetVerboseLevel() {return vlevel;}
    
 private:

  JetScapeLogger() {debug=true;remark=false;vlevel=0;};
  JetScapeLogger(JetScapeLogger const&) {};
  static JetScapeLogger* m_pInstance;

  bool debug;
  bool remark;
  unsigned short vlevel;
  
};

#endif

