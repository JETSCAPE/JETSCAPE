#include "ctrl.h"

#include "inc.h"
#include <ctime>
#include <cmath>
#include <iomanip>
//#include <TMath.h>
#include "hdo.h"
#include "ickw.h"
#include "conv.h"
#include "eo3.h"
#include "eo1.h"
#include "trancoeff.h"
#include <unistd.h>

#include <execinfo.h>
#include <signal.h>

ofstream ofile ; // cout + cerr output

EoS *eos ;
TransportCoeff *trcoeff ;
Fluid *f ;
IC_KW *ic_kw ;
Hydro* h ;
time_t start, end ;

double tau0, tau_max, dtau ;
int maxstep ;
//string directory, outkw ;

void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  cout<<"Error: signal "<<sig<<endl ;
  backtrace_symbols_fd(array, size, 2);
  exit(1);
}

void redirecterrorshlle_(char* filename)
{
    signal(SIGSEGV, handler); // handler for segfault
    ofile.open(filename) ;
    cout.rdbuf(ofile.rdbuf()) ;
    cerr.rdbuf(ofile.rdbuf()) ;
}

void initeoshlle_(char *filename, int* ncols)
{
    eos = new EoSs(filename,*ncols) ; // << CE ; set BAG_ASYMPT !!
}


void initeoshlle3f_(char *filename, double *B, double *volex0, double *delta0, double *aaa, double *bbb)
{
    eos = new EoS3f(filename, *B, *volex0, *delta0, *aaa, *bbb) ;
}

void initeoshlle1f_(char *filename)
{
    eos = new EoS1f(filename) ;
}


void inittrcoeff_(double *etaS, double *zetaS)
{
 trcoeff = new TransportCoeff(*etaS, *zetaS, eos) ;
}

void eoshlle_(double *e, double *nb, double *nq, double *ns, double *T, double *mub, double *muq, double *mus, double *p)
{
	eos->eos(*e, *nb, *nq, *ns, *T, *mub, *muq, *mus, *p) ;
}


void eosrangeshlle_(double *emax, double *e0, double *nmax, double *n0, int *ne, int *nn)
{
	eos->eosranges(*emax, *e0, *nmax, *n0, *ne, *nn) ;
}


void initfluidhlle_(int* nx, int* ny, int* nz, double* minx, double* maxx, double* miny, double* maxy, double* minz, double* maxz)
{
    f = new Fluid(eos, trcoeff, *nx, *ny, *nz, *minx, *maxx, *miny, *maxy, *minz, *maxz) ;
}



void initichlle_(char *filename, double *tau0)
{
    ic_kw = new IC_KW(filename) ; //---------for Klaus
//    ic_kw->setIC(f, eos, *tau0) ; //!!! z-symmetry is included in ic_kw.cpp
}


void icgethlle3f_(double* x, double* y, double* eta, double* e, double* nB, double* nQ, double* nS, double* vx, double* vy, double* vz)
{
	double nu, nd, ns ;
    ic_kw->getICs(*x, *y, *eta, *e, *vx, *vy, *vz, nu, nd, ns) ;
//	if(*e<1.) *e = 0. ;
	*nB = 1./3.*(nu + nd + ns) ;
	*nQ = 1./3.*(2.*nu - nd - ns) ;
	*nS = - ns ;
}


void icgethlle_(double* x, double* y, double* eta, double* e, double* nB, double* nQ, double* nS, double* vx, double* vy, double* vz)
{
	double nu, nd, ns ;
    ic_kw->getICs(*x, *y, *eta, *e, *vx, *vy, *vz, nu, nd, ns) ;
	*nB = *nQ = *nS = 0. ;
}


void icsethlle_(int* ix, int* iy, int* iz, double* tau0, double* e, double* nb, double* nq, double* ns, double* vx, double* vy, double* vz)
{
Cell *c = f->getCell(*ix-1,*iy-1,*iz-1) ;
if((*vx)*(*vx)+(*vy)*(*vy)+(*vz)*(*vz)>1.){
//	cerr << "setIC : " << ix << "  " << iy << "  " << iz << "  e = " << e << "  v^2 = " << (*vx)*(*vx)+(*vy)*(*vy)+vz*vz << endl ;
	double factor = sqrt((*vx)*(*vx)+(*vy)*(*vy)+(*vz)*(*vz)) ;
	(*vx) = (*vx)*0.99/factor ;
	(*vy) = (*vy)*0.99/factor ;
	(*vz) = (*vz)*0.99/factor ;
  }
  if(isinf(*e) || isnan(*e) || isinf(*nb) || isnan(*nb) || isinf(*nq) || isnan(*nq) || isinf(*ns) || isnan(*ns) ||
  isinf(*vx) || isnan(*vx) || isinf(*vy) || isnan(*vy) || isinf(*vz) || isnan(*vz))
  {cout<<"fluid: bad init,"<<setw(14)<<"e"<<setw(14)<<"nb"<<setw(14)<<"nq"<<setw(14)<<"ns"<<setw(14)<<"vx"<<setw(14)<<"vy"<<setw(14)<<"vz"<<endl ;
   cout<<"----------------"<<setw(14)<<e<<setw(14)<<nb<<setw(14)<<nq<<setw(14)<<ns<<setw(14)<<vx<<setw(14)<<vy<<setw(14)<<vz<<endl ;
   cout<<"cell  "<<ix<<"  "<<iy<<"  "<<iz<<endl ;
   exit(1) ; }
  c->setPrimVar(eos, *tau0, *e, *nb, *nq, *ns, (*vx), (*vy), (*vz)) ;
  c->saveQprev() ;
  if(*e>0.) c->setAllM(1.) ;
}


void inithydrohlle_(double* _tau0, double* _tau_max, double* _dtau)
{
    tau0 = *_tau0 ;
    tau_max = *_tau_max ;
    dtau = *_dtau ;
    h = new Hydro(f, eos, trcoeff, tau0, dtau) ;
    maxstep = ceil((tau_max-tau0)/dtau) ;
    
    start = 0;
    time(&start);
    h->setNSvalues() ; // initialize viscous terms
}


void dtauhlle_(double* dtau)
{
  h->setDtau(*dtau) ;
}


void initoutputhlle_(char* dir)
{
  f->initOutput(dir, maxstep, tau0, 2) ;
  f->outputPDirections(h->getTau());
  f->outputTransverseAverages(h->getTau()) ;
}


int getmaxstephlle_(void)
{
    return maxstep ;
}


void makestephlle_(int *i)
{
  h->performStep() ;
	f->outputPDirections(h->getTau());
  f->outputTransverseAverages(h->getTau()) ;
}


void getvalueshlle_(int* ix, int* iy, int* iz, double* e, double *p, double *nb, double *nq, double *ns, double* vx, double* vy, double* vz, double* viscCorrCutFlag)
{
    Cell *c = f->getCell(*ix-1, *iy-1, *iz-1) ;
    c->getPrimVar(eos, h->getTau(), *e, *p, *nb, *nq, *ns, *vx, *vy, *vz) ;
    *viscCorrCutFlag = c->getViscCorrCutFlag() ;
}


void getvflaghlle_(int* ix, int* iy, int* iz, double* viscCorrCutFlag)
{
  *viscCorrCutFlag = f->getCell(*ix-1, *iy-1, *iz-1)->getViscCorrCutFlag() ;
}


void getvischlle_(int* ix, int* iy, int* iz, double *pi, double *Pi)
{
 Cell* c=f->getCell(*ix-1, *iy-1, *iz-1) ;
 // fortran and C have reverse array alignment, (a,b) --> [b,a]
 // but since pi[mu][nu] is symmetric, this is not important
 for(int i=0; i<4; i++)
 for(int j=0; j<4; j++){
  *(pi+4*i+j) = c->getpi(i,j) ;
 }
 *Pi=c->getPi(); 
 //if(c->Pi!=0.0) cout <<"nonzero Pi: " << c->Pi << endl ;
}


double getxhlle_(int *ix)
{
    return f->getX(*ix-1) ;
}


double getyhlle_(int *iy)
{
    return f->getY(*iy-1) ;
}


double getzhlle_(int *iz)
{
    return f->getZ(*iz-1) ;
}


void destroyhlle_(void)
{
//    delete ic ; // do we need IC instance?
    delete f ;
    delete h ;
	end=0 ;
    time(&end); float diff = difftime(end, start);
//    cout<<"Execution time = "<<diff<< " [sec]" << endl;
    ofile.close() ; // close output and error file
}


void destroyeoshlle_(void)
{
	delete eos ;
}


double gettimehlle_(void)
{
    end=0 ;
    time(&end); float diff = difftime(end, start);
    return diff ;
}

double getenergyhlle_(void)
{
    double ene = 0., e, p, nb, nq, ns, vx, vy, vz, _vx, _vy, _vz ;
    double tau = h->getTau() ;
    for(int ix=0; ix<f->getNX(); ix++)
    for(int iy=0; iy<f->getNY(); iy++)
    for(int iz=0; iz<f->getNZ(); iz++){
	f->getCell(ix, iy, iz)->getPrimVar(eos, tau, e, p, nb, nq, ns, _vx, _vy, _vz) ;
	double eta = f->getZ(iz) ;
	vx = _vx/(cosh(eta) + _vz*sinh(eta)) ;
	vy = _vy/(cosh(eta) + _vz*sinh(eta)) ;
	vz = (_vz*cosh(eta)+sinh(eta))/(cosh(eta) + _vz*sinh(eta)) ;
	ene += (e+p)/(1.-vx*vx-vy*vy-vz*vz) - p ;
	}
    return e ;
}
