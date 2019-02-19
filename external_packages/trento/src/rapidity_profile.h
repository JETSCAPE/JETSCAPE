// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// TRENTO3D: Three-dimensional extension of TRENTO by Weiyao Ke
// MIT License

#ifndef ETA_H
#define ETA_H

#include <cmath>
#include <vector>
#include <gsl/gsl_fft_complex.h>

/// GSL fft of complex array, z_i = x_i + j*y_i
/// data are stored in one real array A that A[2i] = x_i and A[2i+1] = y_i
#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

double const sqrt2 = std::sqrt(2);
constexpr double TINY = 1e-12;
constexpr int relative_skew_switch = 1;
constexpr int absolute_skew_switch = 2;
/// the mean of the rapidity distribution as function of Ta(x,y),
/// Tb(x,y) and exp(y_beam).
/// For the case both Ta and Tb are tiny small, mean = 0.0 is returned.
/// For the case only Ta(Tb) = 0, it returns +(-)y_beam,
/// For the case y_beam >> 1, mean ~ 0.5*log(Ta/Tb)
double inline mean_function(double ta, double tb, double exp_ybeam){
  if (ta < TINY && tb < TINY) return 0.;
  return 0.5 * std::log((ta*exp_ybeam + tb/exp_ybeam) / std::max(ta/exp_ybeam + tb*exp_ybeam, TINY));
}

/// The coefficient of normalized width of the rapidity distribution
/// Currently, it is simply a constant independent of Ta(x,y) and Tb(x,y)
double inline std_function(double ta, double tb){
  (void)ta;
  (void)tb;
  return 1.;
}

/// The normalized skewness as function of Ta(x,y) and Tb(x,y)

double inline skew_function(double ta, double tb, int type_switch){
  if (type_switch == relative_skew_switch)
    return (ta - tb)/std::max(ta + tb, TINY);
  else if(type_switch == absolute_skew_switch)
      return ta-tb;
  else return 0.;
}

/// A fast pseudorapidity to rapidity transformer using pretabulated values
/// It returns the transfoamtion y(eta) and the jacobian dy/deta(eta)
class fast_eta2y {
private:
  double etamax_;
  double deta_;
  std::size_t neta_;
  std::vector<double> y_;
  std::vector<double> dydeta_;
public:
  fast_eta2y(double J, double etamax, double deta)
      : etamax_(etamax),
        deta_(deta),
        neta_(std::ceil(2.*etamax_/(deta_+1e-15))+1),
        y_(neta_, 0.),
        dydeta_(neta_, 0.) {

    for (std::size_t ieta = 0; ieta < neta_; ++ieta) {
      double eta = -etamax_ + ieta*deta_;
      double Jsh = J*std::sinh(eta);
      double sq = std::sqrt(1. + Jsh*Jsh);
      y_[ieta] = std::log(sq + Jsh);
      dydeta_[ieta] = J*std::cosh(eta)/sq;
    }
  }

  double rapidity(double eta){
    double steps = (eta + etamax_)/deta_;
    double xi = std::fmod(steps, 1.);
    std::size_t index = std::floor(steps);
    return y_[index]*(1. - xi) + y_[index+1]*xi;
  }

  double Jacobian(double eta){
    double steps = (eta + etamax_)/deta_;
    double xi = std::fmod(steps, 1.);
    std::size_t index = std::floor(steps);
    return dydeta_[index]*(1. - xi) + dydeta_[index+1]*xi;
  }

};


/// A class for cumulant generating function inversion
/// This class handles inverse fourier transform the cumulant generating function
/// A direct transformation F^{-1} exp(i*m*k-(std*k)^2/2-skew*(std*k)^3/6) results
/// in an ill-behaved function. Instead, skew is replaced by skew*exp(-(std*k)^2/2).
/// In this way higher order cumulant are included to regulated the function.
/// The reconstruction range is taken to be [-3.33, 3.33]*std, which does not 
/// seem to be very large, however, by trail and error, this range elimiates 
/// oscillations at large rapidity and leaves the details close to mid-rapidity 
/// unchanged. We used GSL FFT library (more specialized library would be FFTW, 
/// e.g.), with 256 points. The transformed results are stored and interpolated. 
class cumulant_generating{
private:
  size_t const N;
  double * data, * dsdy;
  double eta_max;
  double deta;
  double center;

public:
  cumulant_generating(): N(256), data(new double[2*N]), dsdy(new double[2*N]){};

  /// This function set the mean, std and skew of the profile and use FFT to
  /// transform cumulant generating function at zero mean.
  void calculate_dsdy(double mean, double std, double skew){
    double k1, k2, k3, amp, arg;
    // adaptive eta_max = 3.33*std;
    center = mean;
    eta_max = std*3.33;
    deta = 2.*eta_max/(N-1.);
    double fftmean = eta_max/std;
      for(size_t i=0;i<N;i++){
          k1 = M_PI*(i-N/2.0)/eta_max*std;
      k2 = k1*k1;
      k3 = k2*k1;

          amp = std::exp(-k2/2.0);
          arg = fftmean*k1+skew/6.0*k3*amp;

      REAL(data,i) = amp*std::cos(arg);
          IMAG(data,i) = amp*std::sin(arg);
      }
       gsl_fft_complex_radix2_forward(data, 1, N);

      for(size_t i=0;i<N;i++){
          dsdy[i] = REAL(data,i)*(2.0*static_cast<double>(i%2 == 0)-1.0);
      }
  }

  /// When interpolating the funtion, the mean is put back by simply shifting 
  /// the function by y = y - mean + dy/2, the last term is correcting for
  /// interpolating bin edge instead of bin center
  double interp_dsdy(double y){
    y = y-center+deta/2.;
    if (y < -eta_max || y >= eta_max) return 0.0;
    double xy = (y+eta_max)/deta;
    size_t iy = std::floor(xy);
    double ry = xy-iy;
    return dsdy[iy]*(1.-ry) + dsdy[iy+1]*ry;
  }
};
#endif
