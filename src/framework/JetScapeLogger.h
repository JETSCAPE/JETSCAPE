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

#ifndef JETSCAPELOGGER_H
#define JETSCAPELOGGER_H

#include<iostream>
#include<sstream>
#include <cassert>
#include <iostream>
#include <mutex>
#include <memory>

#include "JetClass.h"

using std::shared_ptr;
using std::make_shared;

// --------------------------------

#define RESET   "\033[0m"
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m" /* Bold White */

// define nicer macros to be used for logging ...
/*
#define INFO  JetScapeLogger::Instance()->Info()<<" " //<<__PRETTY_FUNCTION__<<" : "
#define INFO_NICE  JetScapeLogger::Instance()->InfoNice()
#define JSDEBUG JetScapeLogger::Instance()->Debug()<<__PRETTY_FUNCTION__<<" : "
#define DEBUGTHREAD JetScapeLogger::Instance()->DebugThread()<<__PRETTY_FUNCTION__<<" : "
#define REMARK JetScapeLogger::Instance()->Remark()<<__PRETTY_FUNCTION__<<" : "
#define VERBOSE(l) JetScapeLogger::Instance()->Verbose(l)<<__PRETTY_FUNCTION__<<" : "
#define VERBOSESHOWER(l) JetScapeLogger::Instance()->VerboseShower(l)<<__PRETTY_FUNCTION__<<" : "
#define VERBOSEPARTON(l,p) JetScapeLogger::Instance()->VerboseParton(l,p)<<__PRETTY_FUNCTION__<<" : "
#define VERBOSEPVERTEX(l,v) JetScapeLogger::Instance()->VerboseVertex(l,v)<<__PRETTY_FUNCTION__<<" : "
#define WARN JetScapeLogger::Instance()->Warn()<<__PRETTY_FUNCTION__<<" : "
*/

// define nicer macros to be used for logging and check if they should print stuff ...
// otherwise quite a performance hit ...
#define INFO  JetScapeLogger::Instance()->Info()<<" " //<<__PRETTY_FUNCTION__<<" : "
#define INFO_NICE  JetScapeLogger::Instance()->InfoNice()
#define JSDEBUG if (JetScapeLogger::Instance()->GetDebug()) JetScapeLogger::Instance()->Debug()<<__PRETTY_FUNCTION__<<" : "
#define DEBUGTHREAD if (JetScapeLogger::Instance()->GetDebug()) JetScapeLogger::Instance()->DebugThread()<<__PRETTY_FUNCTION__<<" : "
#define REMARK if (JetScapeLogger::Instance()->GetRemark()) JetScapeLogger::Instance()->Remark()<<__PRETTY_FUNCTION__<<" : "
#define VERBOSE(l) if (l<JetScapeLogger::Instance()->GetVerboseLevel()) JetScapeLogger::Instance()->Verbose(l)<<__PRETTY_FUNCTION__<<" : "
#define VERBOSESHOWER(l) if (l<JetScapeLogger::Instance()->GetVerboseLevel()) JetScapeLogger::Instance()->VerboseShower(l)<<__PRETTY_FUNCTION__<<" : "
#define VERBOSEPARTON(l,p) if (l<JetScapeLogger::Instance()->GetVerboseLevel()) JetScapeLogger::Instance()->VerboseParton(l,p)<<__PRETTY_FUNCTION__<<" : "
#define VERBOSEPVERTEX(l,v) if (l<JetScapeLogger::Instance()->GetVerboseLevel()) JetScapeLogger::Instance()->VerboseVertex(l,v)<<__PRETTY_FUNCTION__<<" : "
#define WARN JetScapeLogger::Instance()->Warn()<<__PRETTY_FUNCTION__<<" : "

namespace Jetscape {

  // Forward declarations. Macro implementation is rather cumbersome and hiccups over include guards
  class Vertex;
  class Parton;
  
// -------------------------------------------
struct safe_ostream {
  struct guarded_impl {
    guarded_impl() = delete;
    guarded_impl(const guarded_impl&) = delete;
    void operator=(const guarded_impl&) = delete;
    guarded_impl(std::ostream& ostream, std::mutex& mutex) : ostream_(ostream), guard_(mutex) {
    }
    ~guarded_impl() {
      ostream_.flush();
    }
    template<typename T> void write(const T& x) {
      ostream_ << x;
    }
    std::ostream& ostream_;
    std::lock_guard<std::mutex> guard_;
  };
  struct impl {
    impl() = delete;
    void operator=(const impl&) = delete;
    impl(std::ostream& ostream, std::mutex& mutex) : unique_impl_(new guarded_impl(ostream, mutex)) {
    }
    impl(const impl& rhs) {
      assert(rhs.unique_impl_.get());
      unique_impl_.swap(rhs.unique_impl_);
    }
    template<typename T> impl& operator<<(const T& x) {
      guarded_impl* p = unique_impl_.get();
      assert(p);
      p->write(x);
      return *this;
    }
    mutable std::unique_ptr<guarded_impl> unique_impl_;
  };
  explicit safe_ostream(std::ostream& ostream) : ostream_(ostream) {
  }
  template<typename T> impl operator<<(const T& x) {
    return impl(ostream_, mutex_) << x;
  }
  std::ostream& ostream_;
  std::mutex mutex_;
};

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
	//*m_dest << m_collector->str() << RESET << std::endl;
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


class LogStreamerThread
{
  
  shared_ptr< std::ostringstream > m_collector;
  safe_ostream* m_dest;

 public:

    LogStreamerThread( safe_ostream& dest )
    {
      m_collector=make_shared<std::ostringstream>();
      m_dest = &dest;
    };
    
    ~LogStreamerThread()
    {    
      if ( m_collector.unique() && m_dest!=nullptr) {	
	*m_dest << m_collector->str() << RESET <<"\n"; 
      }   
    }
    
    template <typename T>
    LogStreamerThread& operator<<( T const& value )
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
  LogStreamer InfoNice();
  LogStreamer Warn();
  LogStreamer Debug();
  LogStreamerThread DebugThread();
  LogStreamer Remark();
  //LogStreamer Error(); //to be implemented
  //LogStreamer Fatal(); //to be implemented
  
  LogStreamer Verbose(unsigned short m_vlevel);
  LogStreamer VerboseShower(unsigned short m_vlevel);
  //Not happy with that fix, still normal << in VERBOSE not working ... follow up.
  LogStreamer VerboseParton(unsigned short m_vlevel,Parton &p);
  LogStreamer VerboseVertex(unsigned short m_vlevel,Vertex &v);
  
  void SetDebug(bool m_debug) {debug=m_debug;}
  void SetRemark(bool m_remark) {remark=m_remark;}
  void SetInfo(bool m_info) {info=m_info;}
  void SetVerboseLevel(unsigned short m_vlevel) {vlevel=m_vlevel;}
  bool GetDebug() {return debug;}  
  bool GetRemark() {return remark;}
  bool GetInfo() {return info;}
  unsigned short GetVerboseLevel() {return vlevel;}
    
 private:

  JetScapeLogger() {info=true;debug=false;remark=false;vlevel=0;};
  JetScapeLogger(JetScapeLogger const&) {};
  static JetScapeLogger* m_pInstance;

  bool debug;
  bool remark;
  bool info;
  unsigned short vlevel;
  
};


} // end namespace Jetscape

#endif
