#include "trancoeff.h"
#include "eos.h"
#include "inc.h"

TransportCoeff::TransportCoeff(double _etaS, double _zetaS, EoS *_eos)
{
 etaS = _etaS ;
 zetaS = _zetaS ;
 eos = _eos ;
}

void TransportCoeff::getEta(double e, double T, double &_etaS, double &_zetaS)
{
 _etaS = etaS ;
 _zetaS=zetaS*(1./3.-eos->cs2(e))/(exp((0.16-T)/0.001)+1.) ;
}


void TransportCoeff::getTau(double T, double &_taupi, double &_tauPi)
{
	if(T>0.) _taupi=std::max(3.0/5.068*etaS/T,0.003) ; else _taupi=0.1 ;
	//if(T>0.) _tauPi=std::max(3.0/5.068*(1./4./C_PI)/T,0.005) ; else _tauPi=0.1 ;
  _tauPi = 0.1 ; // for I-S analytical solution
}
