/* This software is distributed under the GNU Lesser General Public License */
//==========================================================================
//
//   debug.cpp
//
//==========================================================================
// $Id: debug.cpp,v 1.10 2001/11/07 13:58:09 pick Exp $

#include <GTL/debug.h>

#include <fstream>
#include <cstdarg>
#include <cstdio>

#ifdef __GTL_MSVCC
#   ifdef _DEBUG
#	ifndef SEARCH_MEMORY_LEAKS_ENABLED
#	error SEARCH NOT ENABLED
#	endif
#	define new DEBUG_NEW
#	undef THIS_FILE
	static char THIS_FILE[] = __FILE__;
#   endif   // _DEBUG
#endif	// __GTL_MSVCC

__GTL_BEGIN_NAMESPACE

ostream* GTL_debug::GTLerr = 0;

void GTL_debug::debug_message (const char* message, ...) 
{
#ifdef _DEBUG
    va_list arg_list;
    va_start(arg_list, message);

    char buf[1024];
    vsprintf(buf, message, arg_list);
    if (GTLerr) {
	os() << buf;
    }
#endif
}

void GTL_debug::init_debug () 
{
    if (!GTLerr) {
#ifdef __GTL_MSVCC
	GTLerr = new ofstream ("ERRLOG.txt", ios::out | ios::app);
#else
	GTLerr = &cerr;
#endif
    }
}

void GTL_debug::close_debug () 
{
    if (GTLerr) {
#ifdef __GTL_MSVCC 
	((ofstream*) GTLerr)->close();
	delete GTLerr;
	GTLerr = 0;
#endif
    }
}

__GTL_END_NAMESPACE

//--------------------------------------------------------------------------
//   end of file
//--------------------------------------------------------------------------
