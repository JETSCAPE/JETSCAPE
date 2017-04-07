class EoS ;

// this class contains the information about the transport coefficients
// of the fluid: eta/s, zeta/s and the corresponding relaxation times,
// taupi (\tau_\pi) and tauPi (\tau_\Pi)
class TransportCoeff
{
 double etaS, zetaS, taupi, tauPi ;
 EoS *eos ; // EoS instance is needed optionally for zeta/s parametrization, which depends on the speed of sound
 public:
 TransportCoeff(double _etaS, double _zetaS, EoS *_eos) ;
 ~TransportCoeff() {} ;
 // returns (optionally temperature dependent) eta/s and zeta/s
 void getEta(double e, double T, double &_etaS, double &_zetaS) ;
 // returns shear and bulk relaxation times
 void getTau(double T, double &_taupi, double &_tauPi) ;
 // isViscous tells whether the fluid is viscous or inviscid
 inline bool isViscous() { if(etaS>0. || zetaS>0.) return true ; else return false ; }
} ;
