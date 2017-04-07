#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <ctime>
#include <sstream>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>

#include "eos.h"
#include "eo3.h"
#include "eo1.h"
#include "eoChiral.h"
#include "eoHadron.h"
#include "rmn.h"

using namespace std ;


int main(int argc, char **argv)
{
  TApplication theApp("App", &argc, argv);
  // pointers to all the main objects
  EoS *eos ;

  char * eosfile = "eos/Laine_nf3.dat" ;
  int ncols = 3, nrows = 286 ;
  //eos = new EoSs(eosfile,ncols) ;
  eos = new EoSChiral() ;
  EoS* eosH = new EoSHadron("eos/eosHadron3D.dat") ;
  
  double e=0.5, p, nb=0.2, nq=0., ns=0., vx=0.8, vy=0., vz=0. ;
  double Q [7] ;
  transformCV(e, eos->p(e,nb,nq,ns), nb, nq, ns, vx, vy, vz, Q) ;
  transformPV(eosH, Q, e, p, nb, nq, ns, vx, vy, vz) ;
  cout<<"HelloWorld;\n" ;
  transformCV(e, eosH->p(e,nb,nq,ns), nb, nq, ns, vx, vy, vz, Q) ;
  cout<<"HelloWorld;\n" ;

  const int N = 500 ;
  double x[N], y[N] ;
  for(int i=0; i<N; i++){
   double Tt, mubt, muqt, must, pt ;
   double e = i*1.0/N ;
   double nb = i*0.2/N ;
   double nq = i*0.1/N ;
   eosH->eos(0.5, nb, 0., 0., Tt, mubt, muqt, must, pt) ;
   x[i] = nb;
   y[i] = mubt ;
  }
  TGraph *g = new TGraph(N,x,y) ;
  g->SetMarkerStyle(22) ;
  g->SetMarkerSize(0.8) ;
  g->Draw("AP") ;
  theApp.Run() ;

  delete eos ;
  delete eosH ;
}
