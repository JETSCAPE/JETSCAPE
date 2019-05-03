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

#ifndef JETENERGYLOSS_H
#define JETENERGYLOSS_H

#include "JetScapeModuleBase.h"
#include "FluidDynamics.h"
#include "FluidCellInfo.h"
#include "JetClass.h"
#include "JetScapeWriter.h"
#include "PartonShower.h"
#include "PartonPrinter.h"
#include "MakeUniqueHelper.h"
#include "LiquefierBase.h"
#include <vector>
#include <random>

namespace Jetscape {


class JetEnergyLoss : public JetScapeModuleBase, public std::enable_shared_from_this<JetEnergyLoss> 
{
  
 public:
  /** Default constructor. It sets the value of qhat, deltaT and maxT to -99.99, 0.0 and 0.0, respectively. Standard signal slot flags are set to false.
   */
  JetEnergyLoss();
  
  /** A copy constructor for Jet Energy Loss Physics Task.
      @param j A pointer of type JetEnergyLoss class.
   */
  JetEnergyLoss(const JetEnergyLoss &j);

  /** Destructor
   */
  virtual ~JetEnergyLoss();

  /** Deep copy.
      @return Null pointer.
   */
  virtual shared_ptr<JetEnergyLoss> Clone() const {return nullptr;}

  /** It reads the input parameters from a XML file under the tag <Eloss>.
      Sets the Parton class @a inP and PartonShower class @a pShower to null.
      Also initializes the tasks attached to the JetEnergyLoss module.
  */  
  virtual void Init();

  /** It calls DoShower() for all shower-initiating partons.
      To avoid abuse, this can NOT be overwritten. Eloss happens on a parton-by-parton level,
      Exec() should only be executed once per event.
   */
  virtual void Exec() final; // prevents eloss modules from overwrting and missusing

  /** Write output information for each tasks/subtasks attached to the JetEnergyLoss module using JetScapeWriter functionality.
      @param w A pointer of type JetScapeWriter.
  */
  virtual void WriteTask(weak_ptr<JetScapeWriter> w); 

  /** Reset the parton shower information.
  */
  virtual void Clear();

  /** Default function to perform the energy loss for partons at time "time". It should be overridden by different energy loss tasks.
      @param deltaT Step-size.
      @param time Current time.
      @param Q2 Current virtuality of the parton.
      @param pIn Vector of current partons.
      @param pOut Vector of partons at time "time+deltaT".
   */
  virtual void DoEnergyLoss(double deltaT, double time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut) {};
  
  //! Core signal to receive information from the medium
  sigslot::signal5<double, double, double, double, std::unique_ptr<FluidCellInfo>&,multi_threaded_local> GetHydroCellSignal;

  /** For future development. A signal to connect the JetEnergyLoss object to the function UpdateEnergyDeposit() of the FluidDynamics class.
   */
  sigslot::signal2<int, double,multi_threaded_local> jetSignal;

  /** For future development. A signal to connect the JetEnergyLoss object to the function GetEnergyDensity() of the FluidDynamics class. 
   */
  sigslot::signal2<int, double&,multi_threaded_local> edensitySignal;

  /** A signal to connect the JetEnergyLoss object to the function DoEnergyLoss() function.
      Send all a list of shower-initiating partons to all attached eloss modules.
      They in turn decide whether they are responsible or not.
      @TODO Rename...
   */
  sigslot::signal5<double, double, double, vector<Parton>&, vector<Parton>&, multi_threaded_local> SentInPartons;

  /** Sets the value of qhat to "m_qhat". 
      @param m_qhat Jet quenching parameter q-hat.
   */
  void SetQhat(double m_qhat) {qhat=m_qhat;}

  /** @return The current value of qhat.
   */
  const double GetQhat() const {return qhat;}
  
  /** It adds a initiating parton @a p to create the parton shower in an energy loss task.
      @param p A pointer of type parton class.
   */
  void AddShowerInitiatingParton(shared_ptr<Parton> p) {inP=p;}

  /** @return The parton which initiated the parton shower.
   */
  shared_ptr<Parton> GetShowerInitiatingParton() {return inP;}  

  void PrintShowerInitiatingParton();

  /** @return The time-step "deltaT" used by energy loss task. 
   */
  double GetDeltaT() {return deltaT;}

  /** @return The maximum time limit for parton shower.
   */
  double GetMaxT() {return maxT;}

  /** @return The current shower.
   */
  shared_ptr<PartonShower> GetShower() {return pShower;}

  // old test signals ======================
  //! TODO: Remove
  /** Set the flag m_jetSignalConnected to true, if JetEnergyLoss had sent a signal to the function UpdateEnergyDeposit() of the class FluidDynamics.
      @param m_jetSignalConnected A boolean flag.
   */
  void SetJetSignalConnected(bool m_jetSignalConnected) {jetSignalConnected=m_jetSignalConnected;}

  /**  @return A boolean flag. Its status indicates whether JetEnergyLoss had sent a signal to the function UpdateEnergyDeposit() of the class FluidDynamics. 
   */
  const bool GetJetSignalConnected() const {return jetSignalConnected;}  
  
  /** Set the flag m_edensitySignalConnected to true, if JetEnergyLoss had sent a signal to the function GetEnergyDensity() of the class FluidDynamics.
      @param m_edensitySignalConnected A boolean flag. 
   */
  void SetEdensitySignalConnected(bool m_edensitySignalConnected) {edensitySignalConnected=m_edensitySignalConnected;}

  /** 
     @return A boolean flag. Its status indicates whether JetEnergyLoss had sent a signal to the function GetEnergyDensity() of the class FluidDynamics.
   */
  const bool GetEdensitySignalConnected() const {return edensitySignalConnected;}
    
  /** Set the flag m_GetHydroCellSignalConnected to true, if JetEnergyLoss had sent a signal to the function GetHydroCell() of the class FluidDynamics.
      @param m_GetHydroCellSignalConnected A boolean flag.
   */
  void SetGetHydroCellSignalConnected(bool m_GetHydroCellSignalConnected) {GetHydroCellSignalConnected=m_GetHydroCellSignalConnected;}

  /** 
     @return A boolean flag. Its status indicates whether JetEnergyLoss had sent a signal to the function GetHydroCell() of the class FluidDynamics.
   */
  const bool GetGetHydroCellSignalConnected() {return GetHydroCellSignalConnected;}

  /** Set the flag m_SentInPartonsConnected to true, if JetEnergyLoss had sent a signal to the function DoEnergyLoss().
      @param m_SentInPartonsConnected A boolean flag.
   */
  void SetSentInPartonsConnected(bool m_SentInPartonsConnected) {SentInPartonsConnected=m_SentInPartonsConnected;}

  /** 
      @return A boolean flag. Its status indicates whether JetEnergyLoss had sent a signal to the function DoEnergyLoss().
   */
  const bool GetSentInPartonsConnected() {return SentInPartonsConnected;}
    
    void add_a_liqueifier(std::shared_ptr<LiquefierBase> new_liqueifier) {
        liquefier_ptr = new_liqueifier;
    }
  
    // The Slot method to send the vector of Hadronization module
    void SendFinalStatePartons(vector<vector<shared_ptr<Parton>>>& fPartons) {
        fPartons = final_Partons;
    }

 protected:
    std::weak_ptr<LiquefierBase> liquefier_ptr;

  void GetFinalPartonsForEachShower(shared_ptr<PartonShower> shower);


 private:

  double deltaT;
  double maxT;
  
  double qhat;
  shared_ptr<Parton> inP;
  shared_ptr<PartonShower> pShower;	  

  bool GetHydroCellSignalConnected;
  bool SentInPartonsConnected;

  /** This function executes the shower process for the partons produced from the hard scaterring.                                                                         
  */
  void DoShower();

  node vStart;
  node vEnd;

  //old test signals 
  bool jetSignalConnected;
  bool edensitySignalConnected;
  
  // Vector of final state partons for each shower as a vector
  vector<vector<shared_ptr<Parton>>> final_Partons; 

};

} // end namespace Jetscape

#endif
