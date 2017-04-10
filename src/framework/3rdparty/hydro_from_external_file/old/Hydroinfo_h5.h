#ifndef Hydroinfo_h5_H
#define Hydroinfo_h5_H

#include<fstream>
#include<sstream>
#include<string>

#include "hdf5.h"

using namespace std;

struct fluidCell {
   double ed, sd, temperature, pressure;
   double vx, vy, vz;
   double pi[4][4];
   double bulkPi;
};

class HydroinfoH5
{
   private:
      int readinFlag;
      int outputFlag;

      string filename;

      int Visflag;  // flag to determine whether to read evolutions for viscous variables
      int Buffersize;
      hid_t H5file_id, H5groupEventid;

      int grid_XL, grid_XH, grid_YL, grid_YH;
      int grid_Framenum;
      double grid_X0, grid_Y0;
      double grid_Xmax, grid_Ymax;
      double grid_Tau0, grid_dTau, grid_dx, grid_dy;
      double grid_Taumax;

      int grid_LSX, grid_LSY, grid_LST;
      int LST_cur;

      int dimensionX, dimensionY;
      double ***ed, ***sd, ***vx, ***vy, ***Temperature, ***Pressure;
      double ***pi00, ***pi01, ***pi02, ***pi03, ***pi11, ***pi12, ***pi13, ***pi22, ***pi23, ***pi33;
      double ***BulkPi;

   public:
      HydroinfoH5();
      HydroinfoH5(string filename_in, int bufferSize_in, int Visflag_in);
      HydroinfoH5(int XL_in, int XH_in, double DX_in, int LSX_in, int YL_in, int YH_in, double DY_in, int LSY_in, double Tau0_in, double dTau_in, double LST_in, int Visflag_in, string filename_in);

      ~HydroinfoH5();
     
      // functions to write into hdf5 file
      void setHydroFiles(int XL_in, int XH_in, double DX_in, int LSX_in, int YL_in, int YH_in, double DY_in, int LSY_in, double Tau0_in, double dTau_in, double LST_in, int Visflag_in, string filename_in);
      void writeGroupattribute(hid_t H5groupEventid);
      void addGroupattributeInt(hid_t H5groupEventid, string attName, int attValue);
      void addGroupattributeDouble(hid_t H5groupEventid, string attName, double attValue);
      void writeHydroBlock(int Time_id, double **ed_in, double **sd_in, double **p_in, double **Temp_in, double **Vx_in, double **Vy_in, double **Pi00_in, double **Pi01_in, double **Pi02_in, double **Pi03_in, double **Pi11_in, double **Pi12_in, double **Pi13_in, double **Pi22_in, double **Pi23_in, double **Pi33_in, double ** BulkPi_in);
      void CSH5dumpBlockdata(hid_t group_id, const hsize_t * dims, string DatasetName, double** Dataset);

      // functions to read hdf5 files
      void readHydroinfoH5(string filename, int bufferSize_in, int Visflag_in);
      void readHydrogridInfo();
      void printHydrogridInfo();
      int readH5Attribute_int(hid_t id, string attributeName);
      double readH5Attribute_double(hid_t id, string attributeName);

      void readHydroinfoBuffered_total();
      void readHydroinfoSingleframe(int frameIdx);
      void readH5Dataset_double(hid_t id, string datasetName, double** dset_data);

      int getNumberofFrames() {return((int)grid_Framenum);};
      double getHydrogridDX() {return(grid_dx);};
      double getHydrogridDY() {return(grid_dy);};
      double getHydrogridDTau() {return(grid_dTau);};
      double getHydrogridTau0() {return(grid_Tau0);};
      double getHydrogridTaumax() {return(grid_Taumax);};
      double getHydrogridNX() {return(grid_XH - grid_XL + 1);};
      double getHydrogridNY() {return(grid_YH - grid_YL + 1);};
      double getHydrogridX0() {return(grid_X0);};
      double getHydrogridY0() {return(grid_Y0);};
      double getHydrogridXmax() {return(grid_Xmax);};
      double getHydrogridYmax() {return(grid_Ymax);};
      void getHydroinfoOnlattice(int frameIdx, int xIdx, int yIdx, fluidCell* fluidCellptr);
      void getHydroinfo(double tau, double x, double y, fluidCell* fluidCellptr);
      void setZero_fluidCell(fluidCell* fluidCellptr);

      double cubeInterpShell(int idx_x, int idx_y, int idx_z, double x, double y, double z, double ***dataset);
      double cubeInterp(double x, double y, double z, double A000, double A100, double A010, double A110, double A001, double A101, double A011, double A111);

};

#endif
