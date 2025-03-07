/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion
 *collisions
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

#include "JetScapeLogger.h"

#include <fstream>
#include <stddef.h>

using namespace std;

// For Linux systems (also Mac)
// only shows max usage throughtout the program

#include <sys/resource.h>
#include <sys/time.h>

/**
 * @brief Retrieves the memory usage of the current process.
 *
 * This function uses the `getrusage` system call to obtain the memory usage
 * statistics of the current process. The memory usage is reported in kilobytes
 * (kB) on BSD/Linux systems and in bytes on Mac/Darwin systems.
 *
 * @return The maximum resident set size (memory usage) of the current process
 *         in megabytes (MB). If the `getrusage` call fails, the function
 *         returns 0.
 */
long getMemoryUsage() {
  /**
   * Reported in kB on BSD/Linux, bytes in Mac/Darwin
   * Could try to explicitly catch __linux__ as well
   */
  struct rusage usage;
  float mbsize = 1024;
#ifdef __MACH__
  mbsize = 1024 * 1024;
#endif

  if (0 == getrusage(RUSAGE_SELF, &usage))
    return usage.ru_maxrss / mbsize;
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

/** @note Just some definition of colors for colored terminal log output
 * Resetting to default color in LogStreamer!
 */
#define BLACK "\033[30m"              /* Black */
#define RED "\033[31m"                /* Red */
#define GREEN "\033[32m"              //\033[7;30m" bkg     /* Green */
#define YELLOW "\033[33m"             /* Yellow */
#define BLUE "\033[34m"               /* Blue */
#define MAGENTA "\033[35m"            /* Magenta */
#define CYAN "\033[36m"               /* Cyan */
#define WHITE "\033[37m"              /* White */
#define BOLDBLACK "\033[1m\033[30m"   /* Bold Black */
#define BOLDRED "\033[1m\033[31m"     /* Bold Red */
#define BOLDGREEN "\033[1m\033[32m"   /* Bold Green */
#define BOLDYELLOW "\033[1m\033[33m"  /* Bold Yellow */
#define BOLDBLUE "\033[1m\033[34m"    /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m" /* Bold Magenta */
#define BOLDCYAN "\033[1m\033[36m"    /* Bold Cyan */
#define BOLDWHITE "\033[1m\033[37m"   /* Bold White */

#define CLEAR "\033[2J"  // clear screen escape code

namespace Jetscape {

std::ostringstream null;

SafeOstream safe_cout(std::cout);
SafeOstream safe_cerr(std::cerr);
SafeOstream safe_null(null);

JetScapeLogger *JetScapeLogger::m_pInstance = NULL;

/**
 * @brief Singleton instance accessor for JetScapeLogger.
 *
 * This method returns the singleton instance of the JetScapeLogger class.
 * If the instance does not exist, it creates a new one.
 *
 * @return JetScapeLogger* Pointer to the singleton instance of JetScapeLogger.
 */
JetScapeLogger *JetScapeLogger::Instance() {
  if (!m_pInstance)
    m_pInstance = new JetScapeLogger();

  return m_pInstance;
}

/**
 * @brief Logs a warning message.
 *
 * This function creates a LogStreamer object that outputs a warning message
 * to the standard output stream (std::cout). The warning message is prefixed
 * with "[Warning] " and is styled with bold red text.
 *
 * @return A LogStreamer object initialized with the warning message.
 */
LogStreamer JetScapeLogger::Warn() {
  string s = "[Warning] ";
  // s << __PRETTY_FUNCTION__ <<":"<<__LINE__<<" ";
  // return LogStreamer(std::cout<<s<<__PRETTY_FUNCTION__ <<":"<<__LINE__<<" ");
  return LogStreamer(std::cout << BOLDRED << s);
}

/**
 * @brief Creates a debug log streamer thread.
 *
 * This function checks if debugging is enabled and creates a log streamer
 * thread with debug information. If debugging is not enabled, it returns a log
 * streamer thread that writes to a null stream.
 *
 * @return LogStreamerThread A log streamer thread with debug information if
 * debugging is enabled, otherwise a log streamer thread that writes to a null
 * stream.
 */
LogStreamerThread JetScapeLogger::DebugThread() {
  if (debug) {
    string s = "[Debug Thread] ";
    // s <<  __PRETTY_FUNCTION__ <<":"<<__LINE__<<" ";
    s += to_string(getMemoryUsage());
    s += "MB ";
    return (LogStreamerThread(safe_cout) << BLUE << s);
  } else {
    // check if it is not written in some system log files ...
    // safe_null.setstate(std::ios_base::failbit);
    return LogStreamerThread(safe_null);
  }
}

/**
 * @brief Logs a debug message if debugging is enabled.
 *
 * This function checks if debugging is enabled and, if so, constructs a debug
 * message that includes the current memory usage. The message is prefixed with
 * "[Debug]" and is output to the standard output stream with a blue color.
 * If debugging is not enabled, the function returns a LogStreamer with a
 * failed state.
 *
 * @return LogStreamer An object that streams the debug message to the standard
 * output if debugging is enabled, or a failed state stream if debugging is
 * disabled.
 */
LogStreamer JetScapeLogger::Debug() {
  if (debug) {
    string s = "[Debug] ";
    // s <<  __PRETTY_FUNCTION__ <<":"<<__LINE__<<" ";
    s += to_string(getMemoryUsage());
    s += "MB ";
    return LogStreamer(std::cout << BLUE << s);
  } else {
    null.setstate(std::ios_base::failbit);
    return LogStreamer(null);
  }
}

/**
 * @brief Logs an informational message.
 *
 * This function constructs an informational log message, optionally including
 * memory usage information, and returns a LogStreamer object for logging.
 *
 * @return LogStreamer object for logging the informational message. If logging
 *         is disabled, returns a LogStreamer object with a failed state.
 */
LogStreamer JetScapeLogger::Info() {
  string s = "[Info] ";
  // s <<  __PRETTY_FUNCTION__ <<":"<<__LINE__<<" ";
  if (info) {
    string s = "[Info] ";
    // s <<  __PRETTY_FUNCTION__ <<":"<<__LINE__<<" ";
    s += to_string(getMemoryUsage());
    s += "MB ";
    return LogStreamer(std::cout << s);
  } else {
    null.setstate(std::ios_base::failbit);
    return LogStreamer(null);
  }
}

/**
 * @brief Logs an informational message with a nice format.
 *
 * This function creates a log message prefixed with "[Info] ".
 * If the `info` flag is set to true, it returns a LogStreamer object
 * that streams to `std::cout`. Otherwise, it returns a LogStreamer
 * object that streams to a null stream with a fail state.
 *
 * @return LogStreamer object for logging the informational message.
 */
LogStreamer JetScapeLogger::InfoNice() {
  string s = "[Info] ";
  // s <<  __PRETTY_FUNCTION__ <<":"<<__LINE__<<" ";
  // s += to_string(getMemoryUsage()); s+="MB ";
  if (info) {
    return LogStreamer(std::cout << s);
  } else {
    null.setstate(std::ios_base::failbit);
    return LogStreamer(null);
  }
}

/**
 * @brief Logs a remark message.
 *
 * This function returns a LogStreamer object that logs a remark message.
 * If the remark flag is set to true, it prefixes the message with "[REMARK] "
 * and outputs it to the standard output stream with bold magenta formatting.
 * If the remark flag is false, it sets the null stream to a failed state and
 * returns a LogStreamer object associated with the null stream.
 *
 * @return LogStreamer object for logging the remark message.
 */
LogStreamer JetScapeLogger::Remark() {
  if (remark) {
    string s = "[REMARK] ";
    return LogStreamer(std::cout << BOLDMAGENTA << s);
  } else {
    null.setstate(std::ios_base::failbit);
    return LogStreamer(null);
  }
}

/**
 * @brief Logs a verbose message with a specified verbosity level.
 *
 * This function checks if the provided verbosity level is less than the current
 * verbosity level. If it is, it constructs a verbose log message that includes
 * the verbosity level and current memory usage, and returns a LogStreamer
 * object that streams the message to the standard output. If the provided
 * verbosity level is not less than the current verbosity level, it returns a
 * LogStreamer object that is set to a fail state.
 *
 * @param m_vlevel The verbosity level of the message to be logged.
 * @return LogStreamer object that streams the verbose message if the verbosity
 * level is appropriate, otherwise a LogStreamer object in a fail state.
 */
LogStreamer JetScapeLogger::Verbose(unsigned short m_vlevel) {
  if (m_vlevel < vlevel)  // or if (m_vlevel==vlevel)
  {
    string s = "[Verbose][";
    s += std::to_string(m_vlevel);
    s += "] ";
    s += to_string(getMemoryUsage());
    s += "MB ";
    return LogStreamer(std::cout << GREEN << s);
  } else {
    null.setstate(std::ios_base::failbit);
    return LogStreamer(null);
  }
}

/**
 * @brief Logs verbose shower information if the specified verbosity level is
 * less than the current verbosity level.
 *
 * This function creates a log message with a verbosity level and memory usage
 * information. If the specified verbosity level is less than the current
 * verbosity level, it constructs a log message string and returns a LogStreamer
 * object that streams the message to standard output. Otherwise, it returns a
 * LogStreamer object that is set to a failed state.
 *
 * @param m_vlevel The verbosity level to compare against the current verbosity
 * level.
 * @return LogStreamer object that streams the log message if the verbosity
 * level is less, otherwise a LogStreamer object in a failed state.
 */
LogStreamer JetScapeLogger::VerboseShower(unsigned short m_vlevel) {
  if (m_vlevel < vlevel)  // or if (m_vlevel==vlevel)
  {
    string s = "[Verbose][";
    s += std::to_string(m_vlevel);
    s += "] ";
    s += to_string(getMemoryUsage());
    s += "MB ";
    return LogStreamer(std::cout << BOLDCYAN << s);
  } else {
    null.setstate(std::ios_base::failbit);
    return LogStreamer(null);
  }
}

/**
 * @brief Logs verbose information about a Parton object.
 *
 * This function logs detailed information about a given Parton object if the
 * specified verbosity level is less than the current verbosity level of the
 * logger.
 *
 * @param m_vlevel The verbosity level for this log entry.
 * @param p The Parton object to be logged.
 * @return A LogStreamer object that streams the log output.
 */
LogStreamer JetScapeLogger::VerboseParton(unsigned short m_vlevel, Parton &p) {
  if (m_vlevel < vlevel)  // or if (m_vlevel==vlevel)
  {
    string s = "[Verbose][";
    s += std::to_string(m_vlevel);
    s += "] Parton: ";
    return LogStreamer(std::cout << GREEN << s << " " << p << endl);
  } else {
    null.setstate(std::ios_base::failbit);
    return LogStreamer(null);
  }
}

/**
 * @brief Logs a verbose message for a given vertex if the verbosity level is
 * met.
 *
 * @param m_vlevel The verbosity level of the message.
 * @param v The vertex to be logged.
 * @return LogStreamer object that streams the log message if the verbosity
 * level is met, otherwise returns a LogStreamer with a failed state.
 */
LogStreamer JetScapeLogger::VerboseVertex(unsigned short m_vlevel, Vertex &v) {
  if (m_vlevel < vlevel)  // or if (m_vlevel==vlevel)
  {
    string s = "[Verbose][";
    s += std::to_string(m_vlevel);
    s += "] Vertex: ";
    return LogStreamer(std::cout << GREEN << s << " " << v << endl);
  } else {
    null.setstate(std::ios_base::failbit);
    return LogStreamer(null);
  }
}

}  // end namespace Jetscape
