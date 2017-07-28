// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: LongGang Pang (2017)
//                (UC Berkeley and LBNL)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...
#include "../JetScape.h"
#include "../TrentoInitial.h"
#include "../FluidDynamics.h"
#include "gtest/gtest.h"

using namespace Jetscape;


TEST(InitialFromXMLTest, TEST_XML){
  auto jetscape = make_shared<JetScape>("./jetscape_init.xml",3);
  jetscape->SetId("primary");
  auto ini = make_shared<TrentoInitial>();

  jetscape->Add(ini);

  // jetscape->Init() must come before make_shred<JetScapeInitial>() 
  // such that jetscape_init.xml file is read in to memory in prior
  jetscape->Init();

  // Run JetScape with all task/modules as specified ...
  jetscape->Exec();

  // Most thinkgs done in write and clear ...
  jetscape->Finish();
  
}

