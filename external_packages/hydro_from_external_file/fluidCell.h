#ifndef fluidCell_H_
#define fluidCell_H_

<<<<<<< HEAD
<<<<<<< HEAD

=======
>>>>>>> a8fdc27b03dd460fc82996bb1aa469ebf9cbe306
=======
>>>>>>> a8fdc27b03dd460fc82996bb1aa469ebf9cbe306
struct fluidCell {
   double ed, sd, temperature, pressure;
   double vx, vy, vz;
   double pi[4][4];
   double bulkPi;
};

<<<<<<< HEAD
<<<<<<< HEAD

=======
>>>>>>> a8fdc27b03dd460fc82996bb1aa469ebf9cbe306
=======
>>>>>>> a8fdc27b03dd460fc82996bb1aa469ebf9cbe306
struct fluidCell_2D {
    double temperature;
    double ux, uy, ueta;
    // the shear stress tensor are already divided by e+P
    double pi00, pi01, pi02;
    double pi11, pi12;
    double pi22;
    double pi33;
    double bulkPi;
};

struct fluidCell_3D {
    double temperature;
    double vx, vy, vz;
    // the shear stress tensor are already divided by e+P
    double pi00, pi01, pi02, pi03;
    double pi11, pi12, pi13;
    double pi22, pi23;
    double pi33;
    double bulkPi;
};

struct fluidCell_3D_new {
    int itau, ix, iy, ieta;
    double temperature;
    double ux, uy, ueta;
    // the shear stress tensor are already divided by e+P
    double pi11, pi12, pi13;
    double pi22, pi23;
    double bulkPi;
};

#endif  // fluidCell_H_
