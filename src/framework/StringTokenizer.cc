/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * For the full list of contributors see AUTHORS.
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 * or via email to bugs.jetscape.org@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/
// StringTokenizer
// General purpose string tokenizer (C++ string version)
// based on https://github.com/ViDA-NYU/birdvis/blob/master/Tokenizer.cpp

#include "StringTokenizer.h"
#include <iostream>
using std::cout;
using std::endl;

namespace Jetscape {

StringTokenizer::~StringTokenizer()
{
}

bool StringTokenizer::isGraphEntry() const
{
  if (buffer.length()==0)      return false;
  if (buffer.find("[")<1)      return true;
  return false;
}

bool StringTokenizer::isNodeEntry() const
{
  if (buffer.length()==0)      return false;
  if (buffer.find("] V") <100 )  return true;
  return false;
}

bool StringTokenizer::isNodeZero() const
{
  if ( !isGraphEntry() )     return false;
  if ( !isNodeEntry() )      return false;
  if ( isHadronEntry() )     return false;

  if (buffer.find("[0]")<3)	return true;
  
  return false;
}

bool StringTokenizer::isEdgeEntry() const
{
  if (buffer.length()==0)      return false;
  if ( !isGraphEntry() )       return false;
  if (buffer.find("]=>[")<100)   return true;

  return false;
}

bool StringTokenizer::isCommentEntry() const
{
  if (buffer.length()==0)      return false;
  if (buffer.find("#")<1)      return true;
  return false;
}

bool StringTokenizer::isEventEntry() const
{
  if (buffer.length()==0)      return false;
  if (buffer.find("Event")<100)// && !isCommentEntry())
    return true;

  return false;
}

bool StringTokenizer::isHadronEntry() const
{
  if (buffer.length()==0)      return false;
  if (buffer.find("] H") <100)  return true;
  return false;
}

///////////////////////////////////////////////////////////////////////////////
// reset string buffer, delimiter and the currsor position
///////////////////////////////////////////////////////////////////////////////
void StringTokenizer::set(const std::string& str, const std::string& delimiter)
{
    this->buffer = str;
    this->delimiter = delimiter;
    this->currPos = buffer.begin();
}

void StringTokenizer::setString(const std::string& str)
{
    this->buffer = str;
    this->currPos = buffer.begin();
}

void StringTokenizer::setDelimiter(const std::string& delimiter)
{
    this->delimiter = delimiter;
    this->currPos = buffer.begin();
}



///////////////////////////////////////////////////////////////////////////////
// return the next token
// If cannot find a token anymore, return "".
///////////////////////////////////////////////////////////////////////////////
std::string StringTokenizer::next()
{
  //if(buffer.size() <= 0) return "";           // skip if buffer is empty

  if(buffer.size() <= 0) {this->currPos = buffer.end(); return "";}
    
    token.clear();                              // reset token string

    this->skipDelimiter();                      // skip leading delimiters

    // append each char to token string until it meets delimiter
    while(currPos != buffer.end() && !isDelimiter(*currPos))
    {
        token += *currPos;
        ++currPos;
    }
    return token;
}



///////////////////////////////////////////////////////////////////////////////
// skip ang leading delimiters
///////////////////////////////////////////////////////////////////////////////
void StringTokenizer::skipDelimiter()
{
    while(currPos != buffer.end() && isDelimiter(*currPos))
        ++currPos;
}



///////////////////////////////////////////////////////////////////////////////
// return true if the current character is delimiter
///////////////////////////////////////////////////////////////////////////////
bool StringTokenizer::isDelimiter(char c)
{
    return (delimiter.find(c) != std::string::npos);
}



///////////////////////////////////////////////////////////////////////////////
// split the input string into multiple tokens
// This function scans tokens from the current cursor position.
///////////////////////////////////////////////////////////////////////////////
std::vector<std::string> StringTokenizer::split()
{
    std::vector<std::string> tokens;
    std::string token;
    while((token = this->next()) != "")
    {
        tokens.push_back(token);
    }

    return tokens;

}

} // end namespace Jetscape
