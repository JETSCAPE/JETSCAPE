#include <sstream>
#include <iomanip>
#include <algorithm>
#include "simpleLogger.h"
#include "Medium_Reader.h"
#include "random.h"

template <size_t N>
Medium<N>::Medium(std::string filename):
_filename(filename),
_file(filename, H5F_ACC_RDONLY),
_event(_file.openGroup("/Event"))
{
	init();
}

template <size_t N>
void Medium<N>::init(){
	_frame_count = 0;
	hdf5_read_scalar_attr(_event, "DX", _dx);
	hdf5_read_scalar_attr(_event, "DY", _dy);
	hdf5_read_scalar_attr(_event, "Tau0", _tau0);
	hdf5_read_scalar_attr(_event, "dTau", _dtau);
	hdf5_read_scalar_attr(_event, "XL", _iXL);
	hdf5_read_scalar_attr(_event, "XH", _iXH);
	hdf5_read_scalar_attr(_event, "YL", _iYL);
	hdf5_read_scalar_attr(_event, "YH", _iYH);
	_dx *=  fmc_to_GeV_m1;
	_dy *=  fmc_to_GeV_m1;
	_tau0 *=  fmc_to_GeV_m1;
	_dtau *=  fmc_to_GeV_m1;

	_xl = _iXL * _dx;
	_xh = _iXH * _dx;
	_yl = _iYL * _dy;
	_yh = _iYH * _dy;

    hsize_t  num_obj;
    H5Gget_num_objs(_event.getId(), &num_obj);
	_number_of_frames = static_cast<int>(num_obj);
	_shape.clear();
	_shape.resize(N+1); // Nt=2, Nx = 2*iXH+1, Ny = 2*iYH+1
	_shape[0] = 2; _shape[1] = 2*_iXH+1; _shape[2] =  2*_iYH+1;
	std::vector<size_t> spatial_shape(N);
	for (auto i=0; i<N; ++i) spatial_shape[i]=_shape[i+1];	
	_buffer.resize(spatial_shape);
	_Temp.resize(_shape);
	_Vx.resize(_shape);
	_Vy.resize(_shape);
	_power_rank = std::pow(2, N+1);

	_xl_limits.resize(N+1); _xh_limits.resize(N+1); _x_steps.resize(N+1);
	_xl_limits[0] = _tau0; _xh_limits[0] = _tau0 + (_number_of_frames-1)*_dtau;
	_x_steps[0] = _dtau;
	_xl_limits[1] = _xl; _xh_limits[1] = _xh; _x_steps[1] = _dx;
	_xl_limits[2] = _yl; _xh_limits[2] = _yh; _x_steps[2] = _dy;
}

template <size_t N>
bool Medium<N>::load_next(){
	if (_frame_count >= _number_of_frames-1) return false; // time to stop
	
	hsize_t spatial_dims[N];
	H5::DataSet dataset;
	for (auto i=0; i<N; ++i) spatial_dims[i]=_shape[i+1];	
	H5::DataSpace spatial_dataspace(N, spatial_dims);

	// load data
	for (int iTau=0; iTau<2; ++iTau){
		int istart = _buffer.num_elements()*iTau;
		std::stringstream FrameNumber;
		FrameNumber << std::setw(4) << std::setfill('0') << _frame_count+iTau;
		std::string FrameName = "/Event/Frame_"+FrameNumber.str();
		LOG_INFO << "loading " << FrameName;
		
		// Temp
		dataset = _file.openDataSet(FrameName+"/Temp");
		dataset.read(_buffer.data(), H5::PredType::NATIVE_DOUBLE,
					 spatial_dataspace, dataset.getSpace());	
		for(int i=0; i<_buffer.num_elements(); ++i) 
			_Temp.data()[i+istart] = _buffer.data()[i];

		// Vx
		dataset = _file.openDataSet(FrameName+"/Vx");
		dataset.read(_buffer.data(), H5::PredType::NATIVE_DOUBLE,
					 spatial_dataspace, dataset.getSpace());
		for(int i=0; i<_buffer.num_elements(); ++i) 
			_Vx.data()[i+istart] = _buffer.data()[i];
		// Vy
		dataset = _file.openDataSet(FrameName+"/Vy");
		dataset.read(_buffer.data(), H5::PredType::NATIVE_DOUBLE,
					 spatial_dataspace, dataset.getSpace());
		for(int i=0; i<_buffer.num_elements(); ++i) 
			_Vy.data()[i+istart] = _buffer.data()[i];
	}

	//regulate v
	for(int i=0; i< _Vx.num_elements(); ++i){
		double vx = _Vx.data()[i];
		double vy = _Vy.data()[i];
		double vz = 0.0; // at mid-rapidity
		double vabs = std::sqrt(vx*vx + vy*vy + vz*vz);
		if (vabs > 1.-1e-6) {
			double rescale = (1.-1e-6)/vabs;
			vx *= rescale;
			vy *= rescale;	
			vz *= rescale;	
		}
		double gamma = 1.0/std::sqrt(1.- vx*vx - vy*vy - vz*vz);
		_Vx.data()[i] = vx;
		_Vy.data()[i] = vy;
	}
	
	_frame_count ++;
	return true;
}

template <size_t N>
void Medium<N>::interpolate(fourvec x, double & T, double & vx, double & vy, double & vz){
	double tau = std::sqrt(x.t()*x.t() - x.z()*x.z());
	if ( tau < get_tauL()*.99 || tau > get_tauH()*1.01)
		LOG_WARNING << "tau = " << tau << " GeV^{-1} outof considered range ["
				 << get_tauL() << ", " << get_tauH() << "]";
	for(int i=1; i<N+1; ++i) {
		if ( x.a[i] < _xl_limits[i] || x.a[i] > _xh_limits[i] ){
			T = 0.;
			vx = 0.;
			vy = 0.;
			vz = 0.;
			return;
		}
	}

	std::vector<size_t> start_index(N+1);
	std::vector<double> w(N+1);
	for(int i=0; i<N+1; ++i) {
		if (i==0){ // time component
			w[0] = (tau-get_tauL())/_dtau;
			start_index[0] = 0;
		}
		else { // spatial component
			double var = (x.a[i]-_xl_limits[i])/_x_steps[i];
			var = std::min(std::max(var, 0.), _shape[i]-2.);
			size_t n = size_t(floor(var));
			double rx = var-n;
			w[i] = rx;
			start_index[i] = n;
		}
	}
	std::vector<size_t> index(N+1);
	T = 0.;
	vx = 0.;
	vy = 0.;
	vz = x.z()/x.t();
	for(int i=0; i<_power_rank; ++i) {
		double W = 1.0;
		for (int j=0; j<N+1; ++j) {
		    index[j] = start_index[j] + ((i & ( 1 << j )) >> j);
		    W *= (index[j]==start_index[j])?(1.-w[j]):w[j];
		}
		T = T + _Temp(index)*W;
		vx = vx + _Vx(index)*W;
		vy = vy + _Vy(index)*W;
	}
	
	double gammaZ = 1./std::sqrt(1. - vz*vz);
	vx /= gammaZ;
	vy /= gammaZ;
	
	return;
}


////////////////////////// Trento ////////////////////////////////////
TransverPositionSampler::TransverPositionSampler(std::string filename, int iev):
_file(filename, H5F_ACC_RDONLY),
_datasetname("/event_"+std::to_string(iev)),
_event(_file.openGroup(_datasetname))
{
	hdf5_read_scalar_attr(_event, "Nx", _Nx);
	hdf5_read_scalar_attr(_event, "Ny", _Ny);
	hdf5_read_scalar_attr(_event, "dxy", _dx);
	hdf5_read_scalar_attr(_event, "dxy", _dy);
	_dx *=  fmc_to_GeV_m1;
	_dy *=  fmc_to_GeV_m1;
	_x_min = -0.5*_Nx*_dx;
	_y_min = -0.5*_Ny*_dy;
	_x_max = 0.5*_Nx*_dx;
	_y_max = 0.5*_Ny*_dy;
	std::vector<size_t> shape = {_Nx, _Ny};
	_TAB.resize(shape);
	hsize_t dims[2]; dims[0] = _Nx; dims[1] = _Ny;
	H5::DataSet dataset;	
	H5::DataSpace dataspace(2, dims);
	
	dataset = _file.openDataSet(_datasetname+"/Ncoll_density");
	dataset.read(_TAB.data(), H5::PredType::NATIVE_DOUBLE, 
				dataspace, dataset.getSpace());
	//////////// Init PDF ///////////////////
	PDF.clear();
	_TAB_SUM = 0.0;
	for(int i=0;i<_Nx;i++){
		for(int j=0;j<_Ny;j++){
			if(_TAB[i][j] > 0.0){ // skip zeros
				_TAB_SUM += _TAB[i][j];
				ele e1;	
				e1.val = _TAB_SUM; e1.i = i; e1.j = j;
				PDF.push_back(e1);
			}
		}
	}
	//normalize
	for(int i=0;i<PDF.size();i++) PDF[i].val /= _TAB_SUM;
	// Ncoll
	Ncoll = _TAB_SUM * _dx * _dy;
}

void TransverPositionSampler::SampleXY(double & x, double & y){
	double P = Srandom::init_dis(Srandom::gen);
	int low=0, high=PDF.size()-1, mid;
	if (P<PDF[low].val) mid = 0;
	else if (P>=PDF[high].val) mid = PDF.size()-1;
	else{
		do{
			mid = floor((low+high)/2.0);
			if(PDF[low].val <= P && P < PDF[mid].val) high = mid;
			else low = mid;
		}while (low < high-1);
	}
	int ix = PDF[mid].i;
	int iy = PDF[mid].j;
	y = _x_min + _dx*(ix - 0.5 + Srandom::init_dis(Srandom::gen) );
	x = _y_min + _dy*(iy - 0.5 + Srandom::init_dis(Srandom::gen) );
}



template class Medium<2>; // 2+1D hydro
//template class Medium<3>; // 3+1D hydro
