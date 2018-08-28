 #include <vector>
 #include <string>
 #include<iostream>
 #include <memory>

 #include "JetScapeTask.h"
 #include "JetScapeEnergyLossMutex.h"
 
 using namespace std;
 using std::shared_ptr;
 
 namespace Jetscape 
 {
 
 JetScapeEnergyLossMutex::JetScapeEnergyLossMutex()
 {
 }
 
 JetScapeEnergyLossMutex::~JetScapeEnergyLossMutex()
 {
 
 }

 bool JetScapeEnergyLossMutex::CheckEnergyLossModules(vector<shared_ptr<JetScapeTask>> jLossModules)
 {
   bool isMartini = false;
   bool isLbt = false;
   bool isAdscft = false;

   int eLossModuleNum = jLossModules.size();
   if(eLossModuleNum <= 1)  return true;

   for(auto module : jLossModules)
   {
     string name = module->GetId();
     if(!name.compare("Martini")) isMartini = true;
     if(!name.compare("LBT")) isLbt = true;
     if(!name.compare("AdSCFT")) isAdscft = true; 
   }

   if(isMartini && isLbt) return false;
   if(isMartini && isAdscft) return false;
   if(isLbt && isAdscft) return false;
   return true;
 }

 } // end namespace Jetscape
