/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2017
 *
 * For the full list of contributors see AUTHORS.
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 * or via email to bugs.jetscape.org@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

#ifndef HARDPROCESS_H
#define HARDPROCESS_H

#include "InitialState.h"
#include "JetScapeModuleBase.h"
#include "tinyxml2.h"
#include "JetClass.h"
#include <vector>

namespace Jetscape {

  /**
     @class 
     Interface for the hard process. 
   */
class HardProcess : public JetScapeModuleBase 
{
  
 public:

  /** Default constructor to create a Hard Process Physics task. Sets the task ID as "HardProcess".
  */    
  HardProcess();

  /** Standard constructor. Sets the task ID as "HardProcess".
      @param m_name  is a name of the control XML file which contains the input parameters relevant to the hard process under the tag <Hard>.
   */
  HardProcess(string m_name) : JetScapeModuleBase (m_name)
    {SetId("HardProcess");}

  /** Destructor for the Hard Process Physics task.
   */
  virtual ~HardProcess(); 

  /** It reads the input parameters relevant to the hard scattering from the XML file under the name tag <Hard>. Uses JetScapeSingnalManager Instance to retrieve the Initial State Physics information. Calls InitTask(); This explicit call can be used for actual initialization of modules such as @a PythiaGun if attached as a @a polymorphic class. It also initializes the tasks within the current module.
    @sa Read about @a polymorphism in C++.
  */
  virtual void Init();

  /** Calls JetScapeTask::ExecuteTasks() for recursive execution of tasks attached to HardProcess module. It can be overridden by the attached module.
   */
  virtual void Exec();

  /** Erases the hard partons stored in the vector @a hp_list of the hard process module. It can be overridden by the attached module.
  */
  virtual void Clear();

  /** It writes the output information obtained from the HardProcess Task into a file.
      @param w is a pointer of type JetScapeWrite class.
  */
  virtual void WriteTask(weak_ptr<JetScapeWriter> w);

  /** Collect header information for writer modules
      @param w is a pointer of type JetScapeWrite class.
  */
  virtual void CollectHeader( weak_ptr<JetScapeWriter> w );
    
  /** @return A pointer to the XML elements. Such XML elements are the input parameters stored in the XML file under the tag <Hard>.
   */
  tinyxml2::XMLElement* GetHardXML() {return fd;}

  // connect the InitialState module with hard process
  /** A pointer of type InitialState class. 
   */
  std::shared_ptr<InitialState> ini;

  /** 
      @return The number of hard partons.
   */
  int GetNHardPartons() {return hp_list.size();}

  /** @return A pointer to the Parton class for ith hard parton.
      @param i Index of a vector of the hard parton.
   */
  shared_ptr<Parton> GetPartonAt(int i) {return hp_list[i];}

  /** @return A vector of the Parton class. These parton classes correspond to the hard partons.
   */
  vector<shared_ptr<Parton>>& GetPartonList() {return hp_list;}

  /** It adds a parton class pointer p into an existing vector of hard Parton class, and increases the vector size by 1.
      @param p Parton class pointer for a hard parton.
   */
  void AddParton(shared_ptr<Parton> p) {hp_list.push_back(p);}
  
  // Slots ...
  /** This function stores the vector of hard partons into a vector plist.
      @param plist A output vector of Parton class.
   */
  void GetHardPartonList(vector<shared_ptr<Parton>> &plist) {plist=hp_list;}

  /** Generated cross section.
      To be overwritten by implementations that have such information.
  */
  virtual double GetSigmaGen(){ return 1; };
  /** Generated cross section error.
      To be overwritten by implementations that have such information.
  */
  virtual double GetSigmaErr(){ return 0; };

  /** Generated weight.
      This is in addition to sigmaGen, e.g. coming from dynamic oversampling.
      To be overwritten by implementations that have such information.
  */
  virtual double GetEventWeight(){ return 1; };
    
 private:

  tinyxml2::XMLElement *fd;

  // Think of always using unique_ptr for any vector in jetscape framework !???
  // To be discussed ...
  vector<shared_ptr<Parton>> hp_list;

 
};

} // end namespace Jetscape

#endif
