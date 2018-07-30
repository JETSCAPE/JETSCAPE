// main34.cc is a part of the PYTHIA event generator.
// Copyright (C) 2017 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Author: Philip Ilten, December 2015.

// An example where the hard process (p p -> mu+ mu-) is automatically
// produced externally with MadGraph 5, read in, and the remainder of
// the event is then produced by Pythia (MPI, showers, hadronization,
// and decays). A comparison is made between events produced with
// Pythia at LO, MadGraph 5 at LO, and aMC@NLO at NLO.

// For this example to run, MadGraph 5 must be installed and the
// command "exe" (set by default as "mg5_aMC") must be available via
// the command line. Additionally, GZIP support must be enabled via
// the "--with-gzip" configuration option(s). Note that this example has
// only been tested with MadGraph 5 version 2.3.3; due to rapid
// MadGraph development, this example may not work with other
// versions. For more details on the LHAMadgraph class see the
// comments of Pythia8Plugins/LHAMadgraph.h.

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/LHAMadgraph.h"

using namespace Pythia8;

//==========================================================================

// A simple method to run Pythia, analyze the events, and fill a histogram.

void run(Pythia* pythia, Hist& hist, int nEvent) {
  pythia->readString("Random:setSeed = on");
  pythia->readString("Random:seed = 1");
  pythia->init();
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia->next()) continue;
    int iMu1(0), iMu2(0);
    for (int i = 0; i < pythia->event.size(); ++i) {
      if (!iMu1 && pythia->event[i].id() == 13)  iMu1 = i;
      if (!iMu2 && pythia->event[i].id() == -13) iMu2 = i;
      if (iMu1 && iMu2) {
        iMu1 = pythia->event[iMu1].iBotCopyId();
        iMu2 = pythia->event[iMu2].iBotCopyId();
        hist.fill((pythia->event[iMu1].p() + pythia->event[iMu2].p()).pT());
        break;
      }
    }
  }
  pythia->stat();
}

//==========================================================================

int main() {

  // The name of the MadGraph5_aMC@NLO executable.
  // You must prepend this string with the path to the executable
  // on your local installation, or otherwise make it available.
  string exe("mg5_aMC");

  // Create the histograms.
  Hist pyPtZ("Pythia dN/dpTZ", 100, 0., 100.);
  Hist mgPtZ("MadGraph dN/dpTZ", 100, 0., 100.);
  Hist amPtZ("aMC@NLO dN/dpTZ", 100, 0., 100.);

  // Produce leading-order events with Pythia.
  Pythia* pythia = new Pythia();
  pythia->readString("Beams:eCM = 13000.");
  pythia->readString("WeakSingleBoson:ffbar2gmZ = on");
  pythia->readString("23:onMode = off");
  pythia->readString("23:onIfMatch = -13 13");
  pythia->readString("PhaseSpace:mHatMin = 80.");
  run(pythia, pyPtZ, 1000);
  delete pythia;

  // Produce leading-order events with MadGraph 5.
  pythia = new Pythia();
  LHAupMadgraph madgraph(pythia, true, "madgraphrun", exe);
  madgraph.readString("generate p p > mu+ mu-");
  // Note the need for a blank character before "set".
  madgraph.readString(" set ebeam1 6500");
  madgraph.readString(" set ebeam2 6500");
  madgraph.readString(" set mmll 80");
  pythia->setLHAupPtr(&madgraph);
  run(pythia, mgPtZ, 1000);
  delete pythia;

  // Produce next-to-leading-order events with aMC@NLO.
  pythia = new Pythia();
  LHAupMadgraph amcatnlo(pythia, true, "amcatnlorun", exe);
  amcatnlo.readString("generate p p > mu+ mu- [QCD]");
  // Note the need for a blank character before "set".
  amcatnlo.readString(" set ebeam1 6500");
  amcatnlo.readString(" set ebeam2 6500");
  amcatnlo.readString(" set mll 80");
  pythia->setLHAupPtr(&amcatnlo);
  run(pythia, amPtZ, 1000);
  delete pythia;

  // Print the histograms.
  cout << pyPtZ;
  cout << mgPtZ;
  cout << amPtZ;
  return 0;
}
