// JetScape SignalManager init reader class implementation (meant as singelton)

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

/*
void JetScapeSignalManager::ConnectJetSignal(shared_ptr<JetScapeModuleBase> j) //, shared_ptr<JetScapeModuleBase> h)
{
  if (!dynamic_pointer_cast<JetEnergyLoss>(j)->GetJetSignalConnected())
    {
      dynamic_pointer_cast<JetEnergyLoss>(j)->jetSignal.connect(GetHydroPointer().get(),&FluidDynamics::UpdateEnergyDeposit);
      dynamic_pointer_cast<JetEnergyLoss>(j)->SetJetSignalConnected(true);
    }
}
*/

void JetScapeSignalManager::ConnectJetSignal(shared_ptr<JetEnergyLoss> j) //, shared_ptr<JetScapeModuleBase> h)
{  
  if (!j->GetJetSignalConnected())
    {
      j->jetSignal.connect(GetHydroPointer().lock().get(),&FluidDynamics::UpdateEnergyDeposit);
      j->SetJetSignalConnected(true);
      jet_signal_map.emplace(num_jet_signals,(weak_ptr<JetEnergyLoss>) j);
      
      num_jet_signals++;
    }
}

void JetScapeSignalManager::ConnectEdensitySignal(shared_ptr<JetEnergyLoss> j) //, shared_ptr<JetScapeModuleBase> h)
{  
  if (!j->GetEdensitySignalConnected())
    {
      j->edensitySignal.connect(GetHydroPointer().lock().get(),&FluidDynamics::GetEnergyDensity);
      j->SetEdensitySignalConnected(true);
      edensity_signal_map.emplace(num_edensity_signals,(weak_ptr<JetEnergyLoss>) j);
      
      num_edensity_signals++;
    }
}

void JetScapeSignalManager::CleanUp()
{
  //some cleanup ...
  // not sure what actually happens to the "connect"!???
  // needed only if task gets destroyed after connected ...
  
  map<int,weak_ptr<JetEnergyLoss>>::iterator it;
  for (it=jet_signal_map.begin(); it!=jet_signal_map.end();++it)
    {
      if (!it->second.lock().get())
	{
	  //VERBOSE(8) << "[" << it->first << ':' << it->second.lock().get() << ']';	  
	  jet_signal_map.erase(it);
	  num_jet_signals--;
	}
    }

  map<int,weak_ptr<JetEnergyLoss>>::iterator it2;
  for (it2=edensity_signal_map.begin(); it2!=edensity_signal_map.end();++it2)
    {
      if (!it->second.lock().get())
	{
	  //VERBOSE(8) << "[" << it->first << ':' << it->second.lock().get() << ']';	  
	  jet_signal_map.erase(it2);
	  num_edensity_signals--;
	}
    }
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

/*
void JetScapeSignalManager::Clear()
{
  // if use of shared pointers ...
  //hydro=nullptr;
  //jloss=nullptr;
}
*/
/*
*/
