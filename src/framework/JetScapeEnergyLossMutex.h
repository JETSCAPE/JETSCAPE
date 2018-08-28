 #ifndef JETSCAPEENERGYLOSSMUTEX_H
 #define JETSCAPEENERGYLOSSMUTEX_H

 #include <vector>
 #include <memory>
 #include "JetScapeTask.h"
 #include "tinyxml2.h"

 using namespace std; 
 using std::shared_ptr;
 
 
 namespace Jetscape 
 {
 
 class JetScapeEnergyLossMutex
 {
   public:
     JetScapeEnergyLossMutex();
     ~JetScapeEnergyLossMutex();
 
     bool CheckEnergyLossModules(vector<shared_ptr<JetScapeTask>> jLossModules);
 
 };
 
 } // end namespace Jetscape
 
 #endif
