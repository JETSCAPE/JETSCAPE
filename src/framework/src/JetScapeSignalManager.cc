// -----------------------------------------
// JetScape (modular/task) based framework
// Intial Design: Joern Putschke (2017)
//                (Wayne State University)
// -----------------------------------------
// License and Doxygen-like Documentation to be added ...

#include "JetScapeSignalManager.h"
#include "JetScapeLogger.h"
#include <stdlib.h>

using namespace std;

JetScapeSignalManager* JetScapeSignalManager::m_pInstance = NULL;

JetScapeSignalManager* JetScapeSignalManager::Instance()
{
  if (!m_pInstance)
    {
      INFO<<"Created JetScapeSignalManager Instance";
      m_pInstance = new JetScapeSignalManager();
    }
  
  return m_pInstance;
}

void JetScapeSignalManager::ConnectGetHardPartonListSignal(shared_ptr<JetEnergyLossManager> jm)
{
  if (!jm->GetGetHardPartonListConnected() && GetHardProcessPointer().lock().get())
    {
      jm->GetHardPartonList.connect(GetHardProcessPointer().lock().get(),&HardProcess::GetHardPartonList);
      jm->SetGetHardPartonListConnected(true);
    }
}

void JetScapeSignalManager::ConnectJetSignal(shared_ptr<JetEnergyLoss> j) 
{  
  if (!j->GetJetSignalConnected())
    {
      j->jetSignal.connect(GetHydroPointer().lock().get(),&FluidDynamics::UpdateEnergyDeposit);
      j->SetJetSignalConnected(true);
      jet_signal_map.emplace(num_jet_signals,(weak_ptr<JetEnergyLoss>) j);
      
      num_jet_signals++;
    }
}

void JetScapeSignalManager::ConnectEdensitySignal(shared_ptr<JetEnergyLoss> j) 
{  
  if (!j->GetEdensitySignalConnected())
    {
      j->edensitySignal.connect(GetHydroPointer().lock().get(),&FluidDynamics::GetEnergyDensity);
      j->SetEdensitySignalConnected(true);
      edensity_signal_map.emplace(num_edensity_signals,(weak_ptr<JetEnergyLoss>) j);
      
      num_edensity_signals++;
    }
}

void JetScapeSignalManager::ConnectAddJetSourceSignal(shared_ptr<JetEnergyLoss> j)
{
  // to be done ...
}

void JetScapeSignalManager::ConnectGetTemperatureSignal(shared_ptr<JetEnergyLoss> j)
{
  // to be done ... (obsolete with HydroCell) ...
}


void JetScapeSignalManager::ConnectGetHydroCellSignal(shared_ptr<JetEnergyLoss> j)
{
  if (!j->GetGetHydroCellSignalConnected())
    {
      j->GetHydroCellSignal.connect(GetHydroPointer().lock().get(),&FluidDynamics::GetHydroCell);
      j->SetGetHydroCellSignalConnected(true);
      GetHydroCellSignal_map.emplace(num_GetHydroCellSignals,(weak_ptr<JetEnergyLoss>) j);
      
      num_GetHydroCellSignals++;
    }
}

void JetScapeSignalManager::ConnectSentInPartonsSignal(shared_ptr<JetEnergyLoss> j,shared_ptr<JetEnergyLoss> j2)
{
  if (!j2->GetSentInPartonsConnected())
    {
      j->SentInPartons.connect(j2.get(), &JetEnergyLoss::DoEnergyLoss);
      j2->SetSentInPartonsConnected(true);
      SentInPartons_map.emplace(num_SentInPartons,(weak_ptr<JetEnergyLoss>) j2);

      num_SentInPartons++;
    }
}
			   
void JetScapeSignalManager::CleanUp()
{
  VERBOSE(8);

  // hmmm wrong caintainer .. should have used vectore with struct instead of map!!!!
  
  if (jloss.lock())
    { 
      int nEnd=GetHydroCellSignal_map.size();
      int nStart=jloss.lock()->GetTaskAt(0)->GetNumberOfTasks();
      
      for (int i=nStart;i<nEnd;i++)
	{
	  jet_signal_map.erase(i);
	  num_jet_signals--;
	  
	  edensity_signal_map.erase(i);
	  num_edensity_signals--;
	  
	  GetHydroCellSignal_map.erase(i);
	  num_GetHydroCellSignals--;

	  SentInPartons_map.erase(i);
	  num_SentInPartons--;
	}
    }
  else
    {
      jet_signal_map.clear();edensity_signal_map.clear();GetHydroCellSignal_map.clear(),SentInPartons_map.clear();
      // think better here how to handle the clean of when the instance goes out of scope ...!???
    }

  PrintGetHydroCellSignalMap();
  PrintSentInPartonsSignalMap();
  
  VERBOSE(8)<<"Done ...";
}

void JetScapeSignalManager::PrintJetSignalMap()
{
  for (auto& x: jet_signal_map)
    VERBOSE(8) << "[" << x.first << ':' << x.second.lock().get() << ']'<<" "<<x.second.lock()->GetId();
}

void JetScapeSignalManager::PrintEdensitySignalMap()
{
  for (auto& x: edensity_signal_map)
    VERBOSE(8) << "[" << x.first << ':' << x.second.lock().get() << ']'<<" "<<x.second.lock()->GetId();
}

void JetScapeSignalManager::PrintGetHydroCellSignalMap()
{
   for (auto& x: GetHydroCellSignal_map)
     //if (x.second.lock().get())
       VERBOSE(8) << "[" << x.first << ':' << x.second.lock().get() << ']'<<" "<<x.second.lock()->GetId();
}

void JetScapeSignalManager::PrintSentInPartonsSignalMap()
{
   for (auto& x: SentInPartons_map)
     //if (x.second.lock().get())
       VERBOSE(8) << "[" << x.first << ':' << x.second.lock().get() << ']'<<" "<<x.second.lock()->GetId();
}

/*
void JetScapeSignalManager::Clear()
{
  // if use of shared pointers ...
  //hydro=nullptr;
  //jloss=nullptr;
  // ...
}
*/
