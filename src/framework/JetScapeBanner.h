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

#ifndef JETSCAPE_BANNER_H
#define JETSCAPE_BANNER_H

#include "JetScapeLogger.h"

namespace Jetscape {

/**
 * @brief Displays the JETSCAPE banner with ASCII art and license information.
 *
 * This function prints a stylized ASCII banner with the JETSCAPE logo and key
 * details about the framework, including citation information and the licensing terms.
 * The output is handled using the `INFO_NICE` logging macro.
 *
 * @note The function does not take any parameters or return values. It only prints
 *       information to the standard output using the logger.
 */
void ShowJetscapeBanner() {
  INFO_NICE
      << "*--------------------------------------------------------------*";
  INFO_NICE
      << "|                                                              |";
  INFO_NICE << "|                  /" << (char)92
            << "                                          |";
  INFO_NICE << "|                 /  " << (char)92
            << "                                         |";
  INFO_NICE << "|                / |  " << (char)92
            << "                                        |";
  INFO_NICE << "|               /  |   " << (char)92 << "              /"
            << (char)92 << "                       |";
  INFO_NICE << "|              /   |    " << (char)92 << "          /"
            << (char)92 << "/ "
            << " " << (char)92 << "                      |";
  INFO_NICE << "|             / " << (char)92 << "  |  /  " << (char)92
            << "      /" << (char)92 << "/"
            << "   |  " << (char)92 << "/" << (char)92
            << "                   |";
  INFO_NICE << "|            /   " << (char)92 << " | /    " << (char)92
            << "    /    % | %   " << (char)92 << "                  |";
  INFO_NICE << "|         __/     " << (char)92 << "|/      " << (char)92
            << "__/"
            << "      %|%     " << (char)92 << "/" << (char)92
            << "__             |";
  // INFO_NICE<<"|            /   "<<(char) 92<<" | /    "<<(char) 92<<"    / *
  // * | *   "<<(char) 92<<"                  |"; INFO_NICE<<"|          _/
  // "<<(char) 92<<"|/      "<<(char) 92<<"__/"<<"   *  *|*     "<<(char)
  // 92<<"/"<<(char) 92<<"__             |";
  INFO_NICE
      << "|                                                              |";
  INFO_NICE
      << "|                          JETSCAPE                            |";
  INFO_NICE
      << "|                                                              |";
  INFO_NICE
      << "|     The Jet Energy-loss Tomography with a Statistically      |";
  INFO_NICE
      << "|       and Computationally Advanced Program Envelope          |";
  INFO_NICE
      << "|                     http://jetscape.org                      |";
  INFO_NICE
      << "|                                                              |";
  INFO_NICE
      << "| Please cite arXiv:1903.07706 if you use this package for     |";
  INFO_NICE
      << "| scientific work.                                             |";
  INFO_NICE
      << "|                                                              |";
  INFO_NICE
      << "| JETSCAPE is provided without warranty under the terms        |";
  INFO_NICE
      << "| of the GNU GPLv3. It uses xxx code(s).                       |";
  INFO_NICE
      << "| See COPYING file for details.                                |";
  INFO_NICE
      << "|                                                              |";
  INFO_NICE
      << "*--------------------------------------------------------------*";
}

}  // end namespace Jetscape

#endif
