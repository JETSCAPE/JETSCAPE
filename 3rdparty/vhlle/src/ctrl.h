
extern "C" void redirecterrorshlle_(char* filename) ;
                // initeoshlle: simple p(e) table
extern "C" void initeoshlle_(char *filename, int* ncols) ;
                // initeoshlle3f: p(e,nb,nq,ns) table
extern "C" void initeoshlle3f_(char *filename,double *B, double *volex0, double *delta0, double *aaa, double *bbb) ;
extern "C" void initeoshlle1f_(char *filename) ;
extern "C" void inittrcoeff_(double *etaS, double *zetaS) ;
extern "C" void eoshlle_(double *e, double *nb, double *nq, double *ns, double *T, double *mub, double *muq, double *mus, double *p) ;
extern "C" void eosrangeshlle_(double *emax, double *e0, double *nmax, double *n0, int *ne, int *nn) ;
extern "C" void initfluidhlle_(int* nx, int* ny, int* nz, double* minx, double* maxx, double* miny, double* maxy, double* minz, double* maxz) ;
extern "C" void initichlle_(char *filename, double *tau0) ;
extern "C" void icgethlle3f_(double* x, double* y, double* eta, double* e, double* nb, double* nq, double* ns, double* vx, double* vy, double* vz) ;
extern "C" void icgethlle_(double* x, double* y, double* eta, double* e, double* nb, double* nq, double* ns, double* vx, double* vy, double* vz) ;
              //icsethlle: ix=1..nx, iy=1..ny, iz=1..nz
extern "C" void icsethlle_(int* ix, int* iy, int* iz, double *tau0, double* e, double* nb, double* nq, double* ns, double* vx, double* vy, double* vz) ;
extern "C" void inithydrohlle_(double* _tau0, double* _tau_max, double* _dtau) ;
extern "C" void initoutputhlle_(char* dir) ;
extern "C" int getmaxstephlle_(void) ;
extern "C" void dtauhlle_(double* dtau) ;
             // makestephlle: i=0..maxstep-1
extern "C" void makestephlle_(int *i) ;
extern "C" void destroyhlle_(void) ;
extern "C" void destroyeoshlle_(void) ;

           //getvalueshlle: ix=1..nx, iy=1..ny, iz=1..nz
           //               viscCorrCutFlag=1.0 if viscous corrections are ok
extern "C" void getvalueshlle_(int* ix, int* iy, int* iz, double* e, double *p, double *nb, double *nq, double *ns, double* vx, double* vy, double* vz, double* viscCorrCutFlag) ;
           //getvflaghlle: get viscCorrCutFlag only
extern "C" void getvflaghlle_(int* ix, int* iy, int* iz, double* viscCorrCutFlag) ;
          // getvischlle: ix=1..nx, iy=1..ny, iz=1..nz,
          //              pi is 4*4 array, index order=(tau,x,y,eta)
extern "C" void getvischlle_(int* ix, int* iy, int* iz, double *pi, double *Pi) ;

          // ix=1..nx
extern "C" double getxhlle_(int *ix) ;
          // iy=1..ny
extern "C" double getyhlle_(int *iy) ;
          // iz=1..nz
extern "C" double getzhlle_(int *iz) ;

extern "C" double gettimehlle_(void) ;
extern "C" double getenergyhlle_(void) ;

