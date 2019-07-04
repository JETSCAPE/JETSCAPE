#include "lorentz.h"
fourvec measure_perp(fourvec p, fourvec k){
	double kdotp = (k.x()*p.x() + k.y()*p.y() + k.z()*p.z());
	return k - p*(kdotp/p.pabs2()); // k_perp  vector
}

double dot(const fourvec& A, const fourvec& B) {
	return A.t()*B.t() - A.x()*B.x() - A.y()*B.y() - A.z()*B.z();
}
