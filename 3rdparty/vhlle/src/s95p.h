
// initial state for hydrodynamic evolution from the table
// in the format used in VISH2+1 hydro code

namespace s95p {
// energy density -> entropy density usgin s95 EoS
double s95p_s(double e);
// entropy density -> energy density usgin s95 EoS
double s95p_e(double s);

// load the table for the initial entropy density profile
// factor: factor to scale the entropy density
void loadSongIC(char* filename, double factor);
// get the energy density at a given point in the transverse plane
double getSongEps(double x, double y);
}
