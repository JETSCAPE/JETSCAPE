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
// StringTokenizer
// General purpose string tokenizer (C++ string version)
// based on https://github.com/ViDA-NYU/birdvis/blob/master/Tokenizer.cpp

#ifndef STRINGTOKENIZER_H
#define STRINGTOKENIZER_H

#include <string>
#include <vector>

namespace Jetscape {

/** 
 * @brief Default delimiter string.
 * 
 * Contains space, tab, vertical tab, newline, carriage return, form feed, '=', '>', '[', and ']'.
 */
const std::string DEFAULT_DELIMITER = " \t\v\n\r\f=>[]";

/**
 * @class StringTokenizer
 * @brief A utility class for tokenizing strings.
 *
 * This class provides methods to split a string into tokens based on a given delimiter set.
 */
class StringTokenizer {
 public:
  /** @brief Default constructor. */
  StringTokenizer(){};

  /** 
   * @brief Constructor to initialize with a string and an optional delimiter.
   * @param str The input string to tokenize.
   * @param delimiter The delimiter characters (default: DEFAULT_DELIMITER).
   */
  StringTokenizer(const std::string &str,
                  const std::string &delimiter = DEFAULT_DELIMITER);
  
  /** @brief Destructor. */
  ~StringTokenizer();

  /** 
   * @brief Set a new string and delimiter for tokenization.
   * @param str The new input string.
   * @param delimiter The new delimiter string.
   */
  void set(const std::string &str,
           const std::string &delimiter = DEFAULT_DELIMITER);
  
  /** 
   * @brief Set a new input string only.
   * @param str The new input string.
   */
  void setString(const std::string &str);

  /** 
   * @brief Set a new delimiter string only.
   * @param delimiter The new delimiter string.
   */
  void setDelimiter(const std::string &delimiter);

  /** 
   * @brief Get the next token from the string.
   * @return The next token as a string, or an empty string if no more tokens exist.
   */
  std::string next();

  /** 
   * @brief Split the remaining string into tokens.
   * @return A vector containing all tokens from the current cursor.
   */
  std::vector<std::string> split();

  /** 
   * @brief Check if all tokens have been processed.
   * @return True if no more tokens exist, false otherwise.
   */
  bool done() const { return currPos == buffer.end(); }

  /// @name Special format detection functions
  /// Functions to detect specific formats in JetScape ASCII format.
  ///@{

  /** @brief Check if the entry is a graph entry. */
  bool isGraphEntry() const;
  /** @brief Check if the entry is a node entry. */
  bool isNodeEntry() const;
  /** @brief Check if the entry is a node zero entry. */
  bool isNodeZero() const;
  /** @brief Check if the entry is an edge entry. */
  bool isEdgeEntry() const;
  /** @brief Check if the entry is a comment entry. */
  bool isCommentEntry() const;
  /** @brief Check if the entry is an event entry. */
  bool isEventEntry() const;
  /** @brief Check if the entry is a hadron entry. */
  bool isHadronEntry() const;

  ///@}

 private:
  /** @brief Skip leading delimiters. */
  void skipDelimiter();

  /** 
   * @brief Check if a character is a delimiter.
   * @param c The character to check.
   * @return True if the character is a delimiter, false otherwise.
   */
  bool isDelimiter(char c);

  std::string buffer; ///< Input string to be tokenized.
  std::string token; ///< Current token.
  std::string delimiter; ///< Delimiter characters.
  std::string::const_iterator 
    currPos; ///< Iterator pointing to the current position.
};

}  // end namespace Jetscape

#endif  // STRINGTOKENIZER_H

/**
 * @example
 * Example usage of the StringTokenizer class:
 * @code
 * #include <iostream>
 * #include <string>
 * #include "StringTokenizer.h"
 *
 * using std::string;
 * using std::cout;
 * using std::endl;
 *
 * int main(int argc, char* argv[]) {
 *     // Create a StringTokenizer object
 *     StringTokenizer str;
 *     string token;
 *     int counter = 0;
 *
 *     string m_str2="[0]-->[1] 100. 0 0 75.";
 *     cout<<m_str2<<endl;
 *     str.set(m_str2);
 *     cout<<str.isGraphEntry()<<endl;
 *     cout<<str.isEdgeEntry()<<endl;
 *
 *     while (!str.done()) {
 *         token = str.next();
 *         ++counter;
 *         cout << counter << ": " << token << endl;
 *     }
 *
 *     string m_str3="# [-->Just --> a Comment []..";
 *     cout<<m_str3<<endl;
 *     str.set(m_str3);
 *
 *     cout<<str.isGraphEntry()<<endl;
 *     cout<<str.isEdgeEntry()<<endl;
 *     cout<<str.isCommentEntry()<<endl;
 *
 *     return 0;
 * }
 * @endcode
 */