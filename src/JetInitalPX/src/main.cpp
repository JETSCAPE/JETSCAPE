#include <iostream>
#include <complex>
#include <fstream>
#include <math.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

#include "initialPX.hpp"

using namespace std;

int main(){  					

     srand(time(NULL));

     initialPX initialInfo;
     
     for(int i=1; i<=10; i++) {
         initialInfo.next();
         cout << initialInfo.id() << "  " << initialInfo.px() << "  " << initialInfo.py() << "  " << initialInfo.pz() << "  " << initialInfo.e() << "  " << initialInfo.m() << "  " << initialInfo.x() << "  " << initialInfo.y() << "  " << initialInfo.z() << "  " << initialInfo.t() << endl;
     }

}


