/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * For the full list of contributors see AUTHORS.
 * Report issues at https://github.com/amajumder/JETSCAPE-COMP/issues
 * or via email to bugs.jetscape.org@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#ifndef JETENERGYLOSS_H
#define JETENERGYLOSS_H

#include "JetScapeModuleBase.h"
#include "FluidDynamics.h"
#include "JetClass.h"
#include "JetScapeWriter.h"
#include "PartonShower.h"
#include "PartonPrinter.h"
#include "MakeUniqueHelper.h"
#include <vector>
#include <random>

namespace Jetscape {


class JetEnergyLoss : public JetScapeModuleBase, public std::enable_shared_from_this<JetEnergyLoss> 
{
  
 public:
  /** Default constructor. It sets the value of qhat, deltaT and maxT to -99.99, 0.0 and 0.0, respectively. Standard signal slot flags are set to false.
*/
  JetEnergyLoss();

  /** Standard constructor. Default value of qhat is set to -99.99. Standard signal slot flags are set to false.
      @param m_name is a name of the control XML file which contains the input parameters under the tag <Eloss>.
   */
  JetEnergyLoss(string m_name) : JetScapeModuleBase (m_name)
  {qhat=-99.99;SetId("JetEnergyLoss");jetSignalConnected=false;
    edensitySignalConnected=false; AddJetSourceSignalConnected=false;
    GetTemperatureSignalConnected=false; GetHydroCellSignalConnected=false;}

  /** A copy constructor for Jet Energy Loss Physics Task.
      @param j A pointer of type JetEnergyLoss class.
   */
  JetEnergyLoss(const JetEnergyLoss &j);

  /** Destructor for the Jet Energy Loss Physics Task.
   */
  virtual ~JetEnergyLoss();

  /** Clones the JetEnergyLoss class. It can be overridden by other modules.
      @return Null pointer.
   */
  virtual shared_ptr<JetEnergyLoss> Clone() const {return nullptr;}  // const = 0;

  /** It reads the input parameters from a XML file under the tag <Eloss>. It sets the Parton class @a inP and PartonShower class @a pShower to null. It also initializes the tasks attached to the JetEnergyLoss module. It can be overridden by other modules such as Martini, Matter and AdSCFT.
  */  
  virtual void Init();

  /** It calls  DoShower() to execute the parton shower in the medium.
   */
  virtual void Exec() final; // prevents eloss modules from overwrting and missusing

  /** It writes the output information for each tasks/subtasks attached to the JetEnergyLoss module using JetScapeWriter functionality.
      @param w A pointer of type JetScapeWriter.
  */
  virtual void WriteTask(weak_ptr<JetScapeWriter> w); 

  /** It resets the parton shower information associated with the energy loss tasks.
  */
  virtual void Clear();

  //virtual void DoEnergyLoss(double deltaT, double Q2, const vector<Parton>& pIn, vector<Parton>& pOut) {};

  /** Default function to perform the energy loss for partons at time "time". It should be overridden by different energy loss tasks.
      @param deltaT Step-size.
      @param time Current time.
      @param Q2 Current virtuality of the parton.
      @param pIn Vector of current partons.
      @param pOut Vector of partons at time "time+deltaT".
   */
  virtual void DoEnergyLoss(double deltaT, double time, double Q2, vector<Parton>& pIn, vector<Parton>& pOut) {};
  
  // test only ...
  /** A signal to connect the JetEnergyLoss object to the function UpdateEnergyDeposit() of the FluidDynamics class.
   */
  sigslot::signal2<int, double,multi_threaded_local> jetSignal;

  /** A signal to connect the JetEnergyLoss object to the function GetEnergyDensity() of the FluidDynamics class. 
   */
  sigslot::signal2<int, double&,multi_threaded_local> edensitySignal;

  // for brick/gubser test ...  
  sigslot::signal5<double, double, double, double, JetSource,multi_threaded_local> AddJetSourceSignal;
  sigslot::signal5<double, double, double, double, double&,multi_threaded_local> GetTemperatureSignal;
  // sigslot::signal5<double, double, double, double, FluidCellInfo*,multi_threaded_local> GetHydroCellSignal;
  sigslot::signal5<double, double, double, double, std::unique_ptr<FluidCellInfo>&,multi_threaded_local> GetHydroCellSignal;

  // signal to all energy loss modules ... get intial list and delta T ... (think more !???)
  // test first ...
  // deltaT , criteria , list
  //sigslot::signal4<double, double, const vector<Parton>&, vector<Parton>&, multi_threaded_local> SentInPartons;
  /** A signal to connect the JetEnergyLoss object to the function DoEnergyLoss() function.
   */
  sigslot::signal5<double, double, double, vector<Parton>&, vector<Parton>&, multi_threaded_local> SentInPartons;

  sigslot::signal1<vector<Parton>&, multi_threaded_local> GetOutPartons; // probably not needed ... do in SentInPartons with return ...

  /** Sets the value of qhat to "m_qhat". 
      @param m_qhat Jet quenching parameter q-hat.
   */
  void SetQhat(double m_qhat) {qhat=m_qhat;}

  /** @return The current value of qhat.
   */
  const double GetQhat() const {return qhat;}
  
  // old test signals ...
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

  // signals for jetscape (brick test so far, might be extended)

  void SetAddJetSourceSignalConnected(bool m_AddJetSourceSignalConnected) {AddJetSourceSignalConnected=m_AddJetSourceSignalConnected;}
  
  const bool GetAddJetSourceSignalConnected() {return AddJetSourceSignalConnected;}
  
  void SetGetTemperatureSignalConnected(bool m_GetTemperatureSignalConnected) {GetTemperatureSignalConnected=m_GetTemperatureSignalConnected;}
  
  const bool GetGetTemperatureSignalConnected() {return GetTemperatureSignalConnected;}
  
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

  void SetGetOutPartonsConnected(bool m_GetOutPartonsConnected) {GetOutPartonsConnected=m_GetOutPartonsConnected;}
  
  const bool GetGetOutPartonsConnected() {return GetOutPartonsConnected;}

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

 private:

  double deltaT;
  double maxT; // quick fix here ...
  
  double qhat;
  //old test signals 
  bool jetSignalConnected;
  bool edensitySignalConnected;
  
  bool AddJetSourceSignalConnected;
  bool GetTemperatureSignalConnected; //probably not needed if everything via HydroCell
  bool GetHydroCellSignalConnected;
  bool SentInPartonsConnected;
  bool GetOutPartonsConnected;
  
  shared_ptr<Parton> inP;
  //unique_ptr<PartonShower> pShower;
  shared_ptr<PartonShower> pShower;	  

  node vStart;
  node vEnd;

  /** This function executes the shower process for the partons produced from the har\
d scaterring.                                                                         
  */
  void DoShower();

};

} // end namespace Jetscape

#endif
