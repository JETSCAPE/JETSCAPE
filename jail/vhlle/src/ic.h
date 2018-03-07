
class EoS;
class Fluid;

// this class takes care of the initial conditions for hydrodynamic evolution
class IC {
 public:
  IC(char* icInputFile, double s0ScaleFactor);
  ~IC(void);
  // setIC: initializes entire hydro grid at a given initial proper time tau
  void setIC(Fluid *f, EoS *eos, double tau);
};
