#ifndef __CL_IDEAL__
#define __CL_IDEAL__
/*!< Use cpp exceptionn to handel errors */
#define __CL_ENABLE_EXCEPTIONS 
// System includes
//#include <CL/cl.hpp>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <ctime>
#include <algorithm>
#include <map>

#include <random>

#include "Config.h"
#include "opencl_backend.h"

namespace clvisc {

typedef struct
{
    int block_size; // num of threads along one dim on gpu
    int nx;  // number of grids along x
    int ny;  // number of grids along y
    int nz;  // number of grids along etas
    double tau0; // initial thermalization time
    double dt;   // time step
    float dx;   // x step 
    float dy;   // y step 
    float dz;   // etas step
    float etaos_xmin; // parameterized eta/s (T for minimum etaos)
    float etaos_ymin; // parameterized eta/s (minumum etaos)
    float etaos_left_slop; // parameterized eta/s (left slop)
    float etaos_right_slop; // parameterized eta/s (left slop)
    std::string result_directory;
} Config;

/*! \class CLIdeal c++ wrapper for JetScape 
 *  */
class CLIdeal
{
    private:
    double tau_;
    Config cfg_;
    OpenclBackend backend_;
    CompileOption opts_;

    std::vector<cl_real4> h_ev_;
    cl::Buffer d_ev_[3];
    cl::Buffer d_src_;

    // image2d_t for eos_table
    cl::Image2D eos_table_;

    // submax and d_submax is used to compute the maximum
    // energy density of the fluctuating QGP
    std::vector<cl_real> submax_;
    cl::Buffer d_submax_;

    // stores the maximum energy density history
    std::vector<cl_real> max_ed_history_;

	cl::Kernel kernel_kt_src_christoffel_;
	cl::Kernel kernel_kt_src_alongx_;
	cl::Kernel kernel_kt_src_alongy_;
	cl::Kernel kernel_kt_src_alongz_;
	cl::Kernel kernel_update_ev_;
	cl::Kernel kernel_reduction_;

	//cl::Buffer d_pi_;

    void read_eos_table_(std::string fname);

    void initialize_gpu_buffer_();
    
    public:

	CLIdeal(const Config & cfg, std::string device_type, int device_id);

    // read initial energy density from external vector
    void read_ini_ed(const std::vector<cl_real> & ed);

    // read initial ed, vx, vy, vz vector
    void read_ini(const std::vector<cl_real> & ed, 
                  const std::vector<cl_real> & vx, 
                  const std::vector<cl_real> & vy, 
                  const std::vector<cl_real> & vz);

    // step update for Runge-Kutta method, step = {1, 2}
    void step_update(int step);

    void evolve();

	~CLIdeal();

};

} // end namespace clvisc

#endif 


/*! \mainpage 
 *  \section  */

/*! 
 *  \example 
*
 *  \code
 *
 *


 * \endcode
*/


