// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

///////////////////////////////////////////////////////////////////////////////
// StringTokenizer
// =========
// General purpose string tokenizer (C++ string version)
///////////////////////////////////////////////////////////////////////////////

#include "StringTokenizer.h"

namespace Jetscape {

///////////////////////////////////////////////////////////////////////////////
// destructor
///////////////////////////////////////////////////////////////////////////////
StringTokenizer::~StringTokenizer()
{
}

bool StringTokenizer::isGraphEntry() const
{
  if (buffer.length()>0)
    if (buffer.find("[")<1)
      return true;
    else
      return false;
  else
    return false;
}

bool StringTokenizer::isNodeEntry() const
{
  return isGraphEntry();
}

bool StringTokenizer::isNodeZero() const
{
  if (isGraphEntry() && !isEdgeEntry())
    {
      if (buffer.find("[0]")<3)
	return true;
      else
	return false;
    }
  else
    return 0;
      
}

bool StringTokenizer::isEdgeEntry() const
{
  if (buffer.length()>0)
    if (isGraphEntry())
      if (buffer.find("-->")<5)
	return true;
      else
	return false;
    else
      return false;
  else
    return false;
}

bool StringTokenizer::isCommentEntry() const
{
  if (buffer.length()>0)
    if (buffer.find("#")<1)
      return true;
    else
      return false;
  else
    return false;
}

bool StringTokenizer::isEventEntry() const
{
  if (buffer.length()>0)
    if (buffer.find("Event")<100)// && !isCommentEntry())
      return true;
    else
      return false;
  else
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
