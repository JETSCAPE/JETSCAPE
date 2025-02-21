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
// This is a general basic class for hydrodynamics

#ifndef SURFACECELLINFO_H
#define SURFACECELLINFO_H

#include "RealType.h"

namespace Jetscape {
/**
 * @class SurfaceCellInfo
 * @brief A class representing a single cell on a hyper-surface in the simulation of heavy-ion collisions.
 *
 * This class holds data related to the hydrodynamic properties and conditions of a surface cell 
 * in a heavy-ion collision simulation. The information includes variables for energy density, 
 * entropy density, temperature, chemical potentials, shear stress tensor, and other physical properties 
 * relevant to the hydrodynamics of the system.
 */
class SurfaceCellInfo {
 public:
  // data structure for outputing hyper-surface information
  /** 
   * @brief Tau (proper time).
   * @details The proper time of the cell in the surface.
   */
  Jetscape::real tau;
  /** 
   * @brief X coordinate.
   * @details The x-coordinate of the cell in the surface.
   */
  Jetscape::real x;
  /** 
   * @brief Y coordinate.
   * @details The y-coordinate of the cell in the surface.
   */
  Jetscape::real y;
  /** 
   * @brief Eta (spatial rapidity).
   * @details The eta coordinate (spatial rapidity) of the cell.
   */
  Jetscape::real eta;
  /** 
   * @brief Surface vector.
   * @details This is a 4-dimensional vector representing the surface element. 
   *          It is used in the calculation of the flow dynamics.
   */
  Jetscape::real d3sigma_mu[4];    
   /** 
   * @brief Energy density.
   * @details Local energy density in the cell, measured in GeV/fm^3.
   */
  Jetscape::real energy_density;   
   /** 
   * @brief Entropy density.
   * @details Local entropy density in the cell, measured in 1/fm^3.
   */
  Jetscape::real entropy_density;  
  /** 
   * @brief Temperature.
   * @details Local temperature of the cell, measured in GeV.
   */
  Jetscape::real temperature;      
  /** 
   * @brief Thermal pressure.
   * @details Thermal pressure within the cell, measured in GeV/fm^3.
   */
  Jetscape::real pressure;         
  /** 
   * @brief Baryon density.
   * @details Net baryon density in the cell, measured in 1/fm^3.
   */
  Jetscape::real baryon_density;   
  /** 
   * @brief QGP fraction.
   * @details The fraction of the cell in the QGP+HRG phase, indicating the presence of quark-gluon plasma.
   */
  Jetscape::real qgp_fraction;     //!< Fraction of quark gluon plasma assuming
                                   //!< medium is in QGP+HRG phase.
  /** 
   * @brief Baryon chemical potential.
   * @details Net baryon chemical potential in the cell, measured in GeV.
   */                   
  Jetscape::real mu_B;             
  /** 
   * @brief Charge chemical potential.
   * @details Net charge chemical potential in the cell, measured in GeV.
   */
  Jetscape::real mu_Q;             
  /** 
   * @brief Strangeness chemical potential.
   * @details Net strangeness chemical potential in the cell, measured in GeV.
   */
  Jetscape::real mu_S;     
   /** 
   * @brief Flow velocity.
   * @details The flow velocity of the cell in 4 dimensions, used in hydrodynamic simulations.
   */
  Jetscape::real umu[4];   
  /** 
   * @brief Shear stress tensor.
   * @details Shear stress tensor in the cell, representing the viscosity of the medium. 
   *          Measured in GeV/fm^3.
   */
  Jetscape::real pi[10];   
  /** 
   * @brief Bulk viscous pressure.
   * @details Bulk viscous pressure in the cell, measured in GeV/fm^3.
   */
  Jetscape::real bulk_Pi;  

   /** 
   * @brief Default constructor.
   * @details Initializes the object with default values.
   */
  SurfaceCellInfo() = default;

 /** 
   * @brief Destructor.
   * @details Cleans up any allocated resources (if any) when the object is destroyed.
   */
  ~SurfaceCellInfo(){};
};

}  // namespace Jetscape

#endif  // SURFACECELLINFO_H
