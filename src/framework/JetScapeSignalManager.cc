/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
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

#include "JetScapeSignalManager.h"
#include "JetScapeLogger.h"
#include <stdlib.h>

using namespace std;

namespace Jetscape {

JetScapeSignalManager* JetScapeSignalManager::m_pInstance = NULL;

JetScapeSignalManager* JetScapeSignalManager::Instance()
{
  if (!m_pInstance)
    {
      JSINFO<<"Created JetScapeSignalManager Instance";
      m_pInstance = new JetScapeSignalManager();
    }
  
  return m_pInstance;
}

void JetScapeSignalManager::ConnectGetHardPartonListSignal(shared_ptr<JetEnergyLossManager> jm){
  if (!jm->GetGetHardPartonListConnected()){
    auto hpp = GetHardProcessPointer().lock();
    if ( hpp ) {
      jm->GetHardPartonList.connect(hpp.get(),&HardProcess::GetHardPartonList);
      jm->SetGetHardPartonListConnected(true);
    }
  }
}

void JetScapeSignalManager::ConnectGetFinalPartonListSignal(shared_ptr<HadronizationManager> hm) {
  if ( !hm->GetGetFinalPartonListConnected() ){
      auto elp = GetEnergyLossPointer().lock();
      if ( elp ) {
          hm->GetFinalPartonList.connect(elp.get(),&JetEnergyLoss::SendFinalStatePartons);
          hm->SetGetFinalPartonListConnected(true);
    }
  }

  if ( !hm->GetGetHadronListConnected() ){
    auto hpp = GetHardProcessPointer().lock();
    if ( hpp ) {
      hm->GetHadronList.connect(hpp.get(),&HardProcess::GetHadronList);
      hm->SetGetHadronListConnected(true);
    }
  }
}

void JetScapeSignalManager::ConnectJetSignal(shared_ptr<JetEnergyLoss> j) 
{  
  if (!j->GetJetSignalConnected()) {
    auto hp = GetHydroPointer().lock();
    if ( hp ){
      j->jetSignal.connect( hp.get(),&FluidDynamics::UpdateEnergyDeposit);
      j->SetJetSignalConnected(true);
      jet_signal_map.emplace(num_jet_signals,(weak_ptr<JetEnergyLoss>) j);      
      num_jet_signals++;
    }
  }
}

void JetScapeSignalManager::ConnectEdensitySignal(shared_ptr<JetEnergyLoss> j) 
{  
  if (!j->GetEdensitySignalConnected()){
    auto hp = GetHydroPointer().lock();
    if ( hp ){      
      j->edensitySignal.connect(hp.get(),&FluidDynamics::GetEnergyDensity);
      j->SetEdensitySignalConnected(true);
      edensity_signal_map.emplace(num_edensity_signals,(weak_ptr<JetEnergyLoss>) j);      
      num_edensity_signals++;
    }
  }
}

void JetScapeSignalManager::ConnectGetHydroCellSignal(shared_ptr<JetEnergyLoss> j)
{
  if (!j->GetGetHydroCellSignalConnected()) {
    auto hp = GetHydroPointer().lock();
    if ( hp ){
      j->GetHydroCellSignal.connect(hp.get(),&FluidDynamics::GetHydroCell);
      j->SetGetHydroCellSignalConnected(true);
      GetHydroCellSignal_map.emplace(num_GetHydroCellSignals,(weak_ptr<JetEnergyLoss>) j);      
      num_GetHydroCellSignals++;
    }
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

void JetScapeSignalManager::ConnectTransformPartonsSignal(shared_ptr<Hadronization> h,shared_ptr<Hadronization> h2)
{
  if (!h2->GetTransformPartonsConnected())
    {
      h->TransformPartons.connect(h2.get(), &Hadronization::DoHadronization);
      h2->SetTransformPartonsConnected(true);
      TransformPartons_map.emplace(num_TransformPartons,(weak_ptr<Hadronization>) h2);

      num_TransformPartons++;
    }
}
			   
void JetScapeSignalManager::CleanUp()
{
  VERBOSE(8);

  // hmmm wrong caintainer .. should have used vectore with struct instead of map!!!!

  auto loss = jloss.lock();
  if ( loss ) { 
    int nEnd=SentInPartons_map.size();
    int nStart=loss->GetTaskAt(0)->GetNumberOfTasks();
    
    for (int i=nStart;i<nEnd;i++){
      jet_signal_map.erase(i);
      num_jet_signals--;
      
      edensity_signal_map.erase(i);
      num_edensity_signals--;
      
      GetHydroCellSignal_map.erase(i);
      num_GetHydroCellSignals--;
      
      SentInPartons_map.erase(i);
      num_SentInPartons--;  
      
      TransformPartons_map.erase(i);
      num_TransformPartons--;
    }
  } else {
    jet_signal_map.clear();edensity_signal_map.clear();GetHydroCellSignal_map.clear(),SentInPartons_map.clear();
    TransformPartons_map.clear();
    // think better here how to handle the clean of when the instance goes out of scope ...!???
  }

  PrintGetHydroCellSignalMap();
  PrintSentInPartonsSignalMap();
  PrintTransformPartonsSignalMap();  

  VERBOSE(8)<<"Done ...";
}
  
void JetScapeSignalManager::PrintJetSignalMap()
{
  for (auto& x: jet_signal_map){
    auto xs = x.second.lock();
    if ( xs ){
      VERBOSE(8) << "[" << x.first << ':' << xs.get() << ']'<<" "<< xs->GetId();
    }
  }
}

void JetScapeSignalManager::PrintEdensitySignalMap()
{
  for (auto& x: edensity_signal_map){
    auto xs = x.second.lock();
    if ( xs ){
      VERBOSE(8) << "[" << x.first << ':' << xs.get() << ']'<<" "<<xs->GetId();
    }
  }
}

void JetScapeSignalManager::PrintGetHydroCellSignalMap()
{
  for (auto& x: GetHydroCellSignal_map){
    auto xs = x.second.lock();
    if ( xs ){
      VERBOSE(8) << "[" << x.first << ':' << xs.get() << ']'<<" "<<xs->GetId();
    }
  }
}

void JetScapeSignalManager::PrintSentInPartonsSignalMap()
{
  for (auto& x: SentInPartons_map){
    auto xs = x.second.lock();
    if ( xs ){
      VERBOSE(8) << "[" << x.first << ':' << xs.get() << ']'<<" "<<xs->GetId();
    }
  }
}

void JetScapeSignalManager::PrintTransformPartonsSignalMap()
{
  for (auto& x: TransformPartons_map){
    auto xs = x.second.lock();
    if ( xs ){      
      VERBOSE(8) << "[" << x.first << ':' << xs.get() << ']'<<" "<<xs->GetId();
    }
  }
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

} // end namespace Jetscape
