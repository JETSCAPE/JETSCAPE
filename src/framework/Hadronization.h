/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion
 *collisions
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

#ifndef HADRONIZATION_H
#define HADRONIZATION_H

#include "JetClass.h"
#include "JetScapeModuleBase.h"
#include "JetScapeWriter.h"
#include "PartonShower.h"
#include "SurfaceFinder.h"
//#include "MakeUniqueHelper.h"
#include <random>
#include <vector>

namespace Jetscape {

/**
 * @class Hadronization
 * @brief Class responsible for hadronization in the JETSCAPE framework.
 *
 * The Hadronization class handles the transformation of partons into hadrons,
 * managing internal data structures and communication with other framework 
 * components.
 */
class Hadronization : public JetScapeModuleBase,
                      public std::enable_shared_from_this<Hadronization> {
 public:
  /**
   * @brief Constructor.
   */
  Hadronization();

  /**
   * @brief Destructor.
   */
  virtual ~Hadronization();

  /**
   * @brief Clone function (returns nullptr by default).
   * @return nullptr.
   */
  virtual shared_ptr<Hadronization> Clone() const { return nullptr; }

  /**
   * @brief Initialize the hadronization module.
   */
  virtual void Init();

  /**
   * @brief Execute the hadronization process.
   */
  virtual void Exec();

  /**
   * @brief Perform the hadronization process.
   * 
   * @param pIn Input partons.
   * @param hOut Output hadrons.
   * @param pOut Remaining partons after hadronization.
   */
  virtual void DoHadronization(vector<vector<shared_ptr<Parton>>> &pIn,
                               vector<shared_ptr<Hadron>> &hOut,
                               vector<shared_ptr<Parton>> &pOut){};
  
  /**
   * @brief Write hadronization data to an output stream.
   * @param w Writer object.
   */
  virtual void WriteTask(weak_ptr<JetScapeWriter> w);

  /**
   * @brief Clear internal data structures.
   */
  virtual void Clear();

  /**
   * @brief Retrieve the list of hadrons.
   * @param signal Vector to store the hadrons.
   */
  void GetHadrons(vector<shared_ptr<Hadron>> &signal) { signal = outHadrons; }

  /// Signal for parton transformation.
  sigslot::signal3<vector<vector<shared_ptr<Parton>>> &,
                   vector<shared_ptr<Hadron>> &, vector<shared_ptr<Parton>> &,
                   multi_threaded_local>
      TransformPartons;

  /// Signal for obtaining the hydrodynamic hypersurface.
  sigslot::signal2<Jetscape::real, std::vector<SurfaceCellInfo> &,
                   multi_threaded_local>
      GetHydroHyperSurface;

  /// Signal for obtaining hydrodynamic cell information.
  sigslot::signal5<double, double, double, double,
                   std::unique_ptr<FluidCellInfo> &, multi_threaded_local>
      GetHydroCellSignal;

  /**
   * @brief Get the list of hadrons.
   * @return Vector of shared pointers to hadrons.
   */
  vector<shared_ptr<Hadron>> GetHadrons() { return outHadrons; }

  /**
   * @brief Get the list of remaining partons.
   * @return Vector of shared pointers to partons.
   */
  vector<shared_ptr<Parton>> GetOutPartons() { return outPartons; }

  /**
   * @brief Add input partons for hadronization.
   * @param ip Input partons.
   */
  void AddInPartons(vector<vector<shared_ptr<Parton>>> ip) { inPartons = ip; }

  /**
   * @brief Set whether the TransformPartons signal is connected.
   * @param m_TransformPartonsConnected Boolean flag.
   */
  void SetTransformPartonsConnected(bool m_TransformPartonsConnected) {
    TransformPartonsConnected = m_TransformPartonsConnected;
  }

  /**
   * @brief Check if the TransformPartons signal is connected.
   * @return True if connected, otherwise false.
   */
  const bool GetTransformPartonsConnected() {
    return TransformPartonsConnected;
  }

  /**
   * @brief Set whether the GetHydroHyperSurface signal is connected.
   * @param m_GetHydroHyperSurfaceConnected Boolean flag.
   */
  void SetGetHydroHyperSurfaceConnected(bool m_GetHydroHyperSurfaceConnected) {
    HydroHyperSurfaceConnected_ = m_GetHydroHyperSurfaceConnected;
  }

  /**
   * @brief Check if the GetHydroHyperSurface signal is connected.
   * @return True if connected, otherwise false.
   */
  const bool GetGetHydroHyperSurfaceConnected() {
    return HydroHyperSurfaceConnected_;
  }

  /**
   * @brief Set whether the GetHydroCellSignal is connected.
   * @param m_GetHydroCellSignalConnected Boolean flag.
   */
  void SetGetHydroCellSignalConnected(bool m_GetHydroCellSignalConnected) {
    GetHydroCellSignalConnected_ = m_GetHydroCellSignalConnected;
  }

  /**
   * @brief Check if the GetHydroCellSignal signal is connected.
   * @return True if connected, otherwise false.
   */
  const bool GetGetHydroCellSignalConnected() {
    return GetHydroCellSignalConnected_;
  }

  /**
   * @brief Add input hadrons.
   * @param ih Input hadrons.
   */
  void AddInHadrons(vector<shared_ptr<Hadron>> ih) { outHadrons = ih; }

  /**
   * @brief Clear the hadron output vector.
   * 
   * Explanation for necessity is provided in HadronizationManager.h.
   */
  void DeleteHadrons();

  /**
   * @brief Remove hadrons with a positive status flag.
   */
  void DeleteRealHadrons();

 private:
  /// Input partons.
  vector<vector<shared_ptr<Parton>>> inPartons;
  
  /// Output hadrons.
  vector<shared_ptr<Hadron>> outHadrons;

  /// Remaining partons after hadronization.
  vector<shared_ptr<Parton>> outPartons;

  /// Performs the hadronization process.
  void DoHadronize();

  /// Status flags for signal connections.
  bool TransformPartonsConnected;
  bool HydroHyperSurfaceConnected_;
  bool GetHydroCellSignalConnected_;
};

}  // namespace Jetscape

#endif
