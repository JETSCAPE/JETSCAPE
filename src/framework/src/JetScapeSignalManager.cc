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

void JetScapeSignalManager::ConnectGetHardPartonListSignal(shared_ptr<JetEnergyLossManager> jm)
{
  if (!jm->GetGetHardPartonListConnected() && GetHardProcessPointer().lock().get())
    {
      jm->GetHardPartonList.connect(GetHardProcessPointer().lock().get(),&HardProcess::GetHardPartonList);
      jm->SetGetHardPartonListConnected(true);

      // map and signal count later ...
    }

}

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

void JetScapeSignalManager::CleanUp()
{
  VERBOSE(8);
  //some cleanup ...
  // not sure what actually happens to the "connect"!???
  // needed only if task gets destroyed after connected ...
  // smarter !!! seperate private function ...

  //PrintJetSignalMap();
  //map<int,weak_ptr<JetEnergyLoss>>::iterator it;
  //map<int,weak_ptr<JetEnergyLoss>>::reverse_iterator itEnd=jet_signal_map.rend();

  //cout<<jet_signal_map.size()<<endl;
  //for (it=jet_signal_map.rbegin(); it!=itEnd;++it)
  //for (it=jet_signal_map.rbegin(); it!=jet_signal_map.rend();++it)
  /*
  for (auto it = jet_signal_map.cbegin(); it != jet_signal_map.cend(); ++it)
    {
      cout<<"---"<<endl;
      cout<<it->first<<" "<<it->second.lock().get()<<endl;
      cout<<num_jet_signals<<endl;
      cout<< "a) [" << it->first << ':' << it->second.lock().get() << ']'<<endl;
      
      if (!it->second.lock().get())
	{
	  VERBOSE(8) << "a) [" << it->first << ':' << it->second.lock().get() << ']';	  
	  //jet_signal_map.erase(it->first);
	  //jet_signal_map.erase(it); //not working !????
	  num_jet_signals--;
	}
      cout<<num_jet_signals<<endl;
      cout<<"---"<<endl;
    }
    cout<<"After .."<<endl;
  */

  // hmmm wrong caintainer .. should have used vectore with struct instead of map!!!!
  //cout<<GetHydroCellSignal_map.size()<<endl;
  //int n=GetHydroCellSignal_map.size();
  //for (int i=1;i<n;i++)
  if (jloss.lock())
    {
      //cout<<jloss.lock()->GetTaskAt(0)->GetTaskList().size()<<endl;
  
      int nEnd=GetHydroCellSignal_map.size();
      // better acess !???
      int nStart=jloss.lock()->GetTaskAt(0)->GetNumberOfTasks(); //GetTaskList().size();
      
      for (int i=nStart;i<nEnd;i++)
	{
	  jet_signal_map.erase(i);
	  num_jet_signals--;
	  
	  edensity_signal_map.erase(i);
	  num_edensity_signals--;
	  
	  //cout<<i<<endl;
	  GetHydroCellSignal_map.erase(i);
	  num_GetHydroCellSignals--;
	}
    }
  else
    {
      jet_signal_map.clear();edensity_signal_map.clear();GetHydroCellSignal_map.clear();
    }
  // think better here how to handle the clean of when the instance goes out of scope ...!???

  //cout<<GetHydroCellSignal_map.size()<<endl;
  PrintGetHydroCellSignalMap();
  
  /*
  cout<<num_jet_signals<<endl;
  PrintJetSignalMap();
  cout<<"after print"<<endl;
  */
  
  /*
  map<int,weak_ptr<JetEnergyLoss>>::iterator it2;
  for (it2=edensity_signal_map.begin(); it2!=edensity_signal_map.end();++it2)
    {
      if (!it2->second.lock().get())
	{
	  VERBOSE(8) << "b) [" << it->first << ':' << it->second.lock().get() << ']';	  
	  jet_signal_map.erase(it2);
	  num_edensity_signals--;
	}
    }

  map<int,weak_ptr<JetEnergyLoss>>::iterator it3;
  for (it3=GetHydroCellSignal_map.begin(); it3!=GetHydroCellSignal_map.end();++it3)
    {
      if (!it3->second.lock().get())
	{
	  VERBOSE(8) << "c) [" << it->first << ':' << it->second.lock().get() << ']';	  
	  GetHydroCellSignal_map.erase(it3);
	  num_GetHydroCellSignals--;
	}
    }
  */
  
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
