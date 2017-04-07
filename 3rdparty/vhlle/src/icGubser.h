
class EoS;
class Fluid;

// initial state from Gubser analytical solution (ideal hydro)
class ICGubser {
 public:
  ICGubser(void);
  ~ICGubser(void);
  // setIC: initializes entire hydro grid at a given initial proper time tau
  void setIC(Fluid *f, EoS *eos, double tau);
};
