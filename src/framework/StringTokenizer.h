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

///////////////////////////////////////////////////////////////////////////////
// StringTokenizer
// =========
// General purpose string tokenizer (C++ string version)
///////////////////////////////////////////////////////////////////////////////

#ifndef STRINGTOKENIZER_H
#define STRINGTOKENIZER_H

#include <string>
#include <vector>

namespace Jetscape {

// default delimiter string (space, tab, newline, carriage return, form feed and =,>,[,].
const std::string DEFAULT_DELIMITER = " \t\v\n\r\f=>[]";

class StringTokenizer
{
public:
    // ctor/dtor
    StringTokenizer() {};
    StringTokenizer(const std::string& str, const std::string& delimiter=DEFAULT_DELIMITER);
    ~StringTokenizer();

    // set string and delimiter
    void set(const std::string& str, const std::string& delimiter=DEFAULT_DELIMITER);
    void setString(const std::string& str);             // set source string only
    void setDelimiter(const std::string& delimiter);    // set delimiter string only
    
    std::string next();                                 // return the next token, return "" if it ends

    std::vector<std::string> split();                   // return array of tokens from current cursor

    bool done() const { return currPos == buffer.end();}

    // Specific to potential JetScape Ascii format ...
    bool isGraphEntry() const;
    bool isNodeEntry() const;
    bool isNodeZero() const;
    bool isEdgeEntry() const;
    bool isCommentEntry() const;
    bool isEventEntry() const;
    bool isHadronEntry() const;
    
private:
    
    void skipDelimiter();                               // ignore leading delimiters
    bool isDelimiter(char c);                           // check if the current char is delimiter

    std::string buffer;                                 // input string
    std::string token;                                  // output string
    std::string delimiter;                              // delimiter string
    std::string::const_iterator currPos;                // string iterator pointing the current position

};

} // end namespace Jetscape

#endif // STRINGTOKENIZER_H


///////////////////////////////////////////////////////////////////////////////
// Usage of Tokenizer Class: Example program
///////////////////////////////////////////////////////////////////////////////

/*
// testing Tokenizer class

#include "StringTokenizer.h"
#include <string>
#include <iostream>


using std::string;
using std::cout;
using std::endl;

int main(int argc, char* argv[])
{
    // instanciate Tokenizer class
    StringTokenizer str;
    string token;
    int counter = 0;
    
    string m_str2="[0]-->[1] 100. 0 0 75.";
    cout<<m_str2<<endl;
    str.set(m_str2);
    //str.setDelimiter(" ->[]");

    cout<<str.isGraphEntry()<<endl;
    cout<<str.isEdgeEntry()<<endl;

    // Two ways ...
    //while((token = str.next()) != "")
    while (!str.done())
      {
	token = str.next();
        ++counter;
        cout << counter << ": " << token << endl;
      }

    string m_str3="#[-->Just --> a Comment []..";
    cout<<m_str3<<endl;    
    str.set(m_str3);
    
    cout<<str.isGraphEntry()<<endl;
    cout<<str.isEdgeEntry()<<endl;
    cout<<str.isCommentEntry()<<endl;
    
    return 0;
}

*/
