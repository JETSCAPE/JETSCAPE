#ifndef Meidum_Reader_H
#define Meidum_Reader_H

#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <boost/multi_array.hpp>
#include "H5Cpp.h"
#include "lorentz.h"
#include "predefine.h"

// N is the dimension of the grid, we load hydro frame by frame,
// therefore, no need to load the time all at once
template <size_t N> 
class Medium{
private:
	double _dx,_dy,_dtau,_tau0, _xl, _xh, _yl, _yh;
	int _iXL, _iXH, _iYL, _iYH; 
	std::vector<size_t> _shape;
	std::vector<double> _xl_limits, _xh_limits, _x_steps;
	const std::string _filename;
	int _frame_count, _number_of_frames;
	size_t _power_rank;
	H5::H5File _file;
	H5::Group _event;
	boost::multi_array<double, N> _buffer;
    boost::multi_array<double, N+1> _Temp;
    boost::multi_array<double, N+1> _Vx;
	boost::multi_array<double, N+1> _Vy;
public:
	Medium(std::string filename);
	~Medium(){_file.close();};
	void init();
	bool load_next();
	void interpolate(fourvec x, double & T, double & vx, double & vy, double & vz);
	double get_tauL() {return _tau0+_dtau*(_frame_count-1);};
	double get_tauH() {return _tau0+_dtau*_frame_count;};
	double get_hydro_time_step() {return _dtau;};
};

struct ele{
	int i,j;
	double val;
};

class TransverPositionSampler{
private:
	double _dx, _dy, _x_min, _y_min, _x_max, _y_max;
	size_t _Nx, _Ny;
	H5::H5File _file;
	std::string _datasetname;
	H5::Group _event;
	boost::multi_array<double, 2> _TAB;
	std::vector<ele> PDF;
	double _TAB_SUM, Ncoll;

public:
	TransverPositionSampler(std::string filename, int iev);
	void SampleXY(double & x, double & y);
	double get_Ncoll(void) const { return Ncoll;}
};

#endif
