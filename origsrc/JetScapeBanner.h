// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#ifndef JETSCAPE_BANNER_H
#define JETSCAPE_BANNER_H

#include "JetScapeLogger.h"

namespace Jetscape {

void show_jetscape_banner()
{
  INFO_NICE<<"*--------------------------------------------------------------*";
  INFO_NICE<<"|                                                              |";
  INFO_NICE<<"|                  /"<<(char) 92<<"                                          |";
  INFO_NICE<<"|                 /  "<<(char) 92<<"                                         |";
  INFO_NICE<<"|                / |  "<<(char) 92<<"                                        |";
  INFO_NICE<<"|               /  |   "<<(char) 92<<"              /"<<(char) 92<<"                       |";
  INFO_NICE<<"|              /   |    "<<(char) 92<<"          /"<<(char) 92<<"/ "<<" "<<(char) 92<<"                      |";
  INFO_NICE<<"|             / "<<(char) 92<<"  |  /  "<<(char) 92<<"      /"<<(char) 92<<"/"<<"   |  "<<(char) 92<<"/"<<(char) 92<<"                   |";
  INFO_NICE<<"|            /   "<<(char) 92<<" | /    "<<(char) 92<<"    /    % | %   "<<(char) 92<<"                  |";
  INFO_NICE<<"|         __/     "<<(char) 92<<"|/      "<<(char) 92<<"__/"<<"      %|%     "<<(char) 92<<"/"<<(char) 92<<"__             |";
  //INFO_NICE<<"|            /   "<<(char) 92<<" | /    "<<(char) 92<<"    / *  * | *   "<<(char) 92<<"                  |";
  //INFO_NICE<<"|          _/     "<<(char) 92<<"|/      "<<(char) 92<<"__/"<<"   *  *|*     "<<(char) 92<<"/"<<(char) 92<<"__             |";
  INFO_NICE<<"|                                                              |";
  INFO_NICE<<"|                 JETSCAPE beta release 0.1                    |";
  INFO_NICE<<"|                                                              |";
  INFO_NICE<<"|     The Jet Energy-loss Tomography with a Statistically      |";
  INFO_NICE<<"|       and Computationally Advanced Program Envelope          |";
  INFO_NICE<<"|                http://jetscape.wayne.edu                     |";
  INFO_NICE<<"|                                                              |";
  INFO_NICE<<"| Please cite xxx if you use this package for scientific work. |";
  INFO_NICE<<"|                                                              |";
  INFO_NICE<<"| JETSCAPE is provided without warranty under the terms        |";
  INFO_NICE<<"| of the GNU GPLv3. It uses xxx code(s).                       |";
  INFO_NICE<<"| See COPYING file for details.                                |";
  INFO_NICE<<"|                                                              |";
  INFO_NICE<<"*--------------------------------------------------------------*";
}

} // end namespace Jetscape

#endif
