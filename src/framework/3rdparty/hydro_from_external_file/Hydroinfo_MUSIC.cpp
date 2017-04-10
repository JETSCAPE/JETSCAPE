// Hydroinfo_MUSIC.cpp is a part of the MARTINI event generator.
// Copyright (C) 2009-2010 Bjoern Schenke.
// MARTINI is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains routines to read in hydro data from files and functions
// that return interpolated data at a given space-time point

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <vector>
#include <string>

#include "./Hydroinfo_MUSIC.h"

using namespace std;

Hydroinfo_MUSIC::Hydroinfo_MUSIC() {
    hbarC = 0.19733;
    lattice_2D = new vector<fluidCell_2D>;
    lattice_3D = new vector<fluidCell_3D>;
    lattice_3D_new = new vector<fluidCell_3D_new>;
}

Hydroinfo_MUSIC::~Hydroinfo_MUSIC() {
    if (boost_invariant) {
        lattice_2D->clear();
    } else {
        lattice_3D->clear();
        lattice_3D_new->clear();
    }
    delete lattice_2D;
    delete lattice_3D;
    delete lattice_3D_new;
}

void Hydroinfo_MUSIC::readHydroData(int whichHydro, int nskip_tau_in) {
    // all hydro data is stored in tau steps (not t)
    // evolution is converted to tau when accessing the hydro data
    lattice_2D->clear();
    lattice_3D->clear();
    lattice_3D_new->clear();

    // read in setups of the hydro simulation
    ostringstream config_file;
    config_file << "results/music_input";
    ifstream configuration;
    configuration.open(config_file.str().c_str(), ios::in);
    if (!configuration) {
        cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
             << "Unable to open file: " << config_file.str() << endl;
        exit(1);
    }
    string temp1;
    string temp_name;
    while (!configuration.eof()) {
        getline(configuration, temp1);
        stringstream ss(temp1);
        ss >> temp_name;

        // read in grid information
        if (temp_name == "Initial_time_tau_0") {
            ss >> hydroTau0;
        } else if (temp_name == "Delta_Tau") {
            ss >> hydroDtau;
        } else if (temp_name == "X_grid_size_in_fm") {
            double temp;
            ss >> temp;
            hydroXmax = temp/2.;
        } else if (temp_name == "Grid_size_in_x") {
            ss >> ixmax;
        } else if (temp_name == "Eta_grid_size") {
            double temp;
            ss >> temp;
            hydro_eta_max = temp/2.;
        } else if (temp_name == "Grid_size_in_eta") {
            ss >> ietamax;
        } else if (temp_name == "output_evolution_every_N_timesteps") {
            ss >> nskip_tau;
        } else if (temp_name == "output_evolution_every_N_x") {
            ss >> nskip_x;
        } else if (temp_name == "output_evolution_every_N_eta") {
            ss >> nskip_eta;
        }
        // read in additioinal information
        if (temp_name == "Include_Rhob_Yes_1_No_0") {
            ss >> turn_on_rhob;
        } else if (temp_name == "Include_Shear_Visc_Yes_1_No_0") {
            ss >> turn_on_shear;
        } else if (temp_name == "Include_Bulk_Visc_Yes_1_No_0") {
            ss >> turn_on_bulk;
        } else if (temp_name == "turn_on_baryon_diffusion") {
            ss >> turn_on_diff;
        }
    }
    configuration.close();

    hydroDx = 2.*hydroXmax/(ixmax - 1.);
    hydroDeta = 2.*hydro_eta_max/(static_cast<double>(ietamax));

    hydroWhichHydro = whichHydro;
    use_tau_eta_coordinate = 1;

    if (use_tau_eta_coordinate == 0) {
        cout << "Hydroinfo_MUSIC:: Warning hydro grid is set to "
             << "cartesian coordinates, please make sure this is correct!"
             << endl;
    }

    if (whichHydro != 6 && whichHydro != 8 && whichHydro != 9 && whichHydro !=10) {
        cout << "Hydroinfo_MUSIC:: This option is obsolete! whichHydro = "
             << whichHydro << endl;
        exit(1);
    } else if (whichHydro == 6) {
        // 3+1D MUSIC hydro (Schenke, Jeon, Gale)
        cout << "Using 3+1D Jeon Schenke hydro reading data ..." << endl;
        boost_invariant = false;

        ixmax = static_cast<int>(2.*hydroXmax/hydroDx + 0.001);
        ietamax = static_cast<int>(2.*hydro_eta_max/hydroDeta + 0.001);

        // read in temperature, QGP fraction , flow velocity
        // The name of the evolution file: evolution_name
        string evolution_name = "results/evolution_xyeta.dat";
        string evolution_name_Wmunu =
            "results/evolution_Wmunu_over_epsilon_plus_P_xyeta.dat";
        string evolution_name_Pi = "results/evolution_bulk_pressure_xyeta.dat";
        cout << "Evolution file name = " << evolution_name << endl;
        ifstream fin;
        fin.open(evolution_name.c_str(), ios::in);
        if (!fin) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "Unable to open file: " << evolution_name << endl;
            exit(1);
        }
        ifstream fin1;
        fin1.open(evolution_name_Wmunu.c_str(), ios::in);
        if (!fin) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "Unable to open file: " << evolution_name_Wmunu << endl;
            exit(1);
        }
        ifstream fin2;
        fin2.open(evolution_name_Pi.c_str(), ios::in);
        if (!fin) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "Unable to open file: " << evolution_name_Pi << endl;
            exit(1);
        }

        double T, vx, vy, vz, QGPfrac;
        fluidCell_3D newCell;
        int ik = 0;
        while (!fin.eof()) {
            ik++;
            fin >> T;
            fin >> QGPfrac;
            fin >> vx;
            fin >> vy;
            fin >> vz;
            newCell.temperature = T;
            newCell.vx = vx;
            newCell.vy = vy;
            newCell.vz = vz;

            lattice_3D->push_back(newCell);
            if (ik%50000 == 0)
                cout << "o" << flush;
        }
        cout << ik << endl;
        fin.close();
        fin1.close();
        fin2.close();
    } else if (whichHydro == 8) {
        // event-by-event (2+1)-d MUSIC hydro from JF
        // there are two slices in medium in eta_s
        // one at eta_s = -15. and the other at eta_s = 0.0
        // only the medium at middle rapidity will be kept in the memory
        boost_invariant = true;
        cout << "Reading event-by-event hydro evolution data from JF ..."
             << endl;

        ixmax = static_cast<int>(2.*hydroXmax/hydroDx + 0.001);
        ietamax = 1;
        nskip_tau = nskip_tau_in;

        hydroDx *= nskip_x;
        hydroDtau *= nskip_tau;
        hydroDeta *= nskip_eta;

        int n_eta = 2;  // there are two slices in eta_s
        // number of fluid cell in the transverse plane
        int num_fluid_cell_trans = ixmax*ixmax;

        // read in hydro evolution
        string evolution_name = "results/evolution_xyeta.dat";
        string evolution_name_Wmunu =
                "results/evolution_Wmunu_over_epsilon_plus_P_xyeta.dat";
        string evolution_name_Pi = "results/evolution_bulk_pressure_xyeta.dat";

        std::FILE *fin;
        string evolution_file_name = evolution_name;
        cout << "Evolution file name = " << evolution_file_name << endl;
        fin = std::fopen(evolution_file_name.c_str(), "rb");

        if (fin == NULL) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "Unable to open file: " << evolution_file_name << endl;
            exit(1);
        }

        std::FILE *fin1 = NULL;
        if (turn_on_shear == 1) {
            fin1 = std::fopen(evolution_name_Wmunu.c_str(), "rb");
            if (fin1 == NULL) {
                cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                     << "Unable to open file: " << evolution_name_Wmunu << endl;
                exit(1);
            }
        }

        std::FILE *fin2 = NULL;
        if (turn_on_bulk == 1) {
            fin2 = std::fopen(evolution_name_Pi.c_str(), "rb");
            if (fin2 == NULL) {
                cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                     << "Unable to open file: " << evolution_name_Pi << endl;
                exit(1);
            }
        }

        int ik = 0;
        fluidCell_2D newCell;
        float T, QGPfrac, vx, vy, vz;
        double ux, uy, ueta;
        float pi00 = 0.0;
        float pi01 = 0.0;
        float pi02 = 0.0;
        float pi03 = 0.0;
        float pi11 = 0.0;
        float pi12 = 0.0;
        float pi13 = 0.0;
        float pi22 = 0.0;
        float pi23 = 0.0;
        float pi33 = 0.0;;
        float bulkPi = 0.0;
        float e_plus_P = 1e-15;
        float cs2 = 0.0;
        int size = sizeof(float);
        while (true) {
            int status = 0;
            status = std::fread(&T, size, 1, fin);
            status += std::fread(&QGPfrac, size, 1, fin);
            status += std::fread(&vx, size, 1, fin);
            status += std::fread(&vy, size, 1, fin);
            status += std::fread(&vz, size, 1, fin);

            if (status != 5) {  // this is the end of file
                break;
            }

            int status_pi = 0;
            if (turn_on_shear == 1) {
                status_pi = std::fread(&pi00, size, 1, fin1);
                status_pi += std::fread(&pi01, size, 1, fin1);
                status_pi += std::fread(&pi02, size, 1, fin1);
                status_pi += std::fread(&pi03, size, 1, fin1);
                status_pi += std::fread(&pi11, size, 1, fin1);
                status_pi += std::fread(&pi12, size, 1, fin1);
                status_pi += std::fread(&pi13, size, 1, fin1);
                status_pi += std::fread(&pi22, size, 1, fin1);
                status_pi += std::fread(&pi23, size, 1, fin1);
                status_pi += std::fread(&pi33, size, 1, fin1);
            
                if (status_pi != 10) {
                    cout << "Error:Hydroinfo_MUSIC::readHydroData: "
                         << "Wmunu file does not have the same number of "
                         << "fluid cells as the ideal file!" << endl;
                    exit(1);
                }
            }

            int status_bulkPi = 0;
            if (turn_on_bulk == 1) {
                status_bulkPi = std::fread(&bulkPi, size, 1, fin2);
                status_bulkPi += std::fread(&e_plus_P, size, 1, fin2);
                status_bulkPi += std::fread(&cs2, size, 1, fin2);
                
                if (status_bulkPi != 3) {
                    cout << "Error:Hydroinfo_MUSIC::readHydroData: "
                         << "bulkPi file does not have the same number of "
                         << "fluid cells as the ideal file!" << endl;
                    exit(1);
                }
            }

            int ieta_idx = static_cast<int>(ik/num_fluid_cell_trans) % n_eta;
            int itau_idx = static_cast<int>(ik/(num_fluid_cell_trans*n_eta));
            ik++;
            if (itau_idx%nskip_tau != 0)  // skip in tau
                continue;

            // print out tau information
            double tau_local = hydroTau0 + itau_idx*hydroDtau/nskip_tau;
            if ((ik-1)%(num_fluid_cell_trans*n_eta) == 0) {
                cout << "read in tau frame: " << itau_idx
                     << " tau_local = " << setprecision(3) << tau_local
                     << " fm ..."<< endl;
            }

            if (ieta_idx == (n_eta-1)) {
                // store the hydro medium at eta_s = 0.0
                double v2 = vx*vx + vy*vy + vz*vz;
                if (v2 > 1.0) {
                    cerr << "[Hydroinfo_MUSIC::readHydroData:] Error: "
                         << "v > 1! vx = " << vx << ", vy = " << vy
                         << ", vz = " << vz << ", T = " << T << endl;
                    if (T > 0.01) {
                        exit(1);
                    } else {
                        v2 = 0.0;
                    }
                }
                double gamma = 1./sqrt(1. - v2);
                ux = gamma*vx;
                uy = gamma*vy;
                ueta = gamma*vz;  // assuming eta = 0

                newCell.temperature = T;
                // convert vx and vy to longitudinal co-moving frame
                newCell.ux = ux;
                newCell.uy = uy;
                newCell.ueta = ueta;

                // pi^\mu\nu tensor
                newCell.pi00 = pi00;
                newCell.pi01 = pi01;
                newCell.pi02 = pi02;
                newCell.pi11 = pi11;
                newCell.pi12 = pi12;
                newCell.pi22 = pi22;
                newCell.pi33 = pi33;

                // bulk pressure
                if (T > 0.18) {
                    // QGP phase prefactor is divided out here
                    newCell.bulkPi = bulkPi/(15.*(1./3. - cs2)*e_plus_P);
                } else {
                    newCell.bulkPi = bulkPi;   // [1/fm^4]
                }
                lattice_2D->push_back(newCell);
            }
        }
        std::fclose(fin);
        if (turn_on_shear == 1) {
            std::fclose(fin1);
        }
        if (turn_on_bulk == 1) {
            std::fclose(fin2);
        }
        cout << endl;
        cout << "number of fluid cells: " << lattice_2D->size() << endl;
    } else if (whichHydro == 9) {
        // event-by-event (2+1)-d MUSIC hydro
        // the output medium is at middle rapidity
        boost_invariant = true;
        cout << "Reading event-by-event hydro evolution data "
             << "from (2+1)D MUSIC ..." << endl;

        ietamax = 1;

        hydroDx *= nskip_x;
        hydroDtau *= nskip_tau;
        hydroDeta *= nskip_eta;

        nskip_tau = nskip_tau_in;
        ixmax = static_cast<int>(2.*hydroXmax/hydroDx + 0.001);

        int n_eta = 1;
        // number of fluid cell in the transverse plane
        int num_fluid_cell_trans = ixmax*ixmax;

        // read in hydro evolution
        string evolution_name = "results/evolution_xyeta.dat";
        string evolution_name_Wmunu =
                "results/evolution_Wmunu_over_epsilon_plus_P_xyeta.dat";
        string evolution_name_Pi = "results/evolution_bulk_pressure_xyeta.dat";

        std::FILE *fin;
        string evolution_file_name = evolution_name;
        cout << "Evolution file name = " << evolution_file_name << endl;
        fin = std::fopen(evolution_file_name.c_str(), "rb");

        if (fin == NULL) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "Unable to open file: " << evolution_file_name << endl;
            exit(1);
        }

        std::FILE *fin1 = NULL;
        if (turn_on_shear == 1) {
            fin1 = std::fopen(evolution_name_Wmunu.c_str(), "rb");
            if (fin1 == NULL) {
                cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                     << "Unable to open file: " << evolution_name_Wmunu << endl;
                exit(1);
            }
        }

        std::FILE *fin2 = NULL;
        if (turn_on_bulk == 1) {
            fin2 = std::fopen(evolution_name_Pi.c_str(), "rb");
            if (fin2 == NULL) {
                cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                     << "Unable to open file: " << evolution_name_Pi << endl;
                exit(1);
            }
        }

        int ik = 0;
        fluidCell_2D newCell;
        double T, QGPfrac, ux, uy, ueta;
        double vx, vy, vz;
        float pi00 = 0.0;
        float pi01 = 0.0;
        float pi02 = 0.0;
        float pi03 = 0.0;
        float pi11 = 0.0;
        float pi12 = 0.0;
        float pi13 = 0.0;
        float pi22 = 0.0;
        float pi23 = 0.0;
        float pi33 = 0.0;;
        float bulkPi = 0.0;
        float e_plus_P = 1e-15;
        float cs2 = 0.0;
        int size = sizeof(double);
        while (true) {
            int status = 0;
            status = std::fread(&T, size, 1, fin);
            status += std::fread(&QGPfrac, size, 1, fin);
            status += std::fread(&vx, size, 1, fin);
            status += std::fread(&vy, size, 1, fin);
            status += std::fread(&vz, size, 1, fin);
            if (status != 5) {  // this is the end of file
                break;
            }

            double v2 = vx*vx + vy*vy + vz*vz;
            if (v2 > 1.) {
                cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                     << "v > 1! vx = " << vx << ", vy = " << vy
                     << ", vz = " << vz << endl;
                exit(1);
            }
            double gamma = 1./sqrt(1. - v2);
            ux = vx*gamma;
            uy = vy*gamma;
            ueta = vz*gamma;  // assuming at the eta = 0
            
            int status_pi = 0;
            if (turn_on_shear == 1) {
                status_pi = std::fread(&pi00, size, 1, fin1);
                status_pi += std::fread(&pi01, size, 1, fin1);
                status_pi += std::fread(&pi02, size, 1, fin1);
                status_pi += std::fread(&pi03, size, 1, fin1);
                status_pi += std::fread(&pi11, size, 1, fin1);
                status_pi += std::fread(&pi12, size, 1, fin1);
                status_pi += std::fread(&pi13, size, 1, fin1);
                status_pi += std::fread(&pi22, size, 1, fin1);
                status_pi += std::fread(&pi23, size, 1, fin1);
                status_pi += std::fread(&pi33, size, 1, fin1);
                
                if (status_pi != 10) {
                    cout << "Error:Hydroinfo_MUSIC::readHydroData: "
                         << "Wmunu file does not have the same number of "
                         << "fluid cells as the ideal file!" << endl;
                    exit(1);
                }
            }

            int status_bulkPi = 0;
            if (turn_on_bulk == 1) {
                status_bulkPi = std::fread(&bulkPi, size, 1, fin2);
                status_bulkPi += std::fread(&e_plus_P, size, 1, fin2);
                status_bulkPi += std::fread(&cs2, size, 1, fin2);
                
                if (status_bulkPi != 3) {
                    cout << "Error:Hydroinfo_MUSIC::readHydroData: "
                         << "bulkPi file does not have the same number of "
                         << "fluid cells as the ideal file!" << endl;
                    exit(1);
                }
            }

            int ieta_idx = static_cast<int>(ik/num_fluid_cell_trans) % n_eta;
            int itau_idx = static_cast<int>(ik/(num_fluid_cell_trans*n_eta));
            ik++;
            if (itau_idx%nskip_tau != 0)  // skip in tau
                continue;

            // print out tau information
            double tau_local = hydroTau0 + itau_idx*hydroDtau/nskip_tau;
            if ((ik-1)%(num_fluid_cell_trans*n_eta) == 0) {
                cout << "read in tau frame: " << itau_idx
                     << " tau_local = " << setprecision(3) << tau_local
                     << " fm ..."<< endl;
            }

            if (ieta_idx == (n_eta-1)) {
                newCell.temperature = T;
                newCell.ux = ux;
                newCell.uy = uy;
                newCell.ueta = ueta;

                // pi^\mu\nu tensor
                newCell.pi00 = pi00;
                newCell.pi01 = pi01;
                newCell.pi02 = pi02;
                newCell.pi11 = pi11;
                newCell.pi12 = pi12;
                newCell.pi22 = pi22;
                newCell.pi33 = pi33;

                // bulk pressure
                if (T > 0.18) {
                    // QGP phase prefactor is divided out here
                    newCell.bulkPi = bulkPi/(15.*(1./3. - cs2)*e_plus_P);
                } else {
                    newCell.bulkPi = bulkPi;   // [1/fm^4]
                }
                lattice_2D->push_back(newCell);
            }
        }
        std::fclose(fin);
        if (turn_on_shear == 1) {
            std::fclose(fin1);
        }
        if (turn_on_bulk == 1) {
            std::fclose(fin2);
        }
        cout << endl;
        cout << "number of fluid cells: " << lattice_2D->size() << endl;
    } else if (whichHydro == 10) {
        // new 3+1D MUSIC hydro (Schenke, Jeon, Gale, Shen)
        cout << "Using 3+1D new MUSIC hydro reading data ..." << endl;
        boost_invariant = false;

        ixmax = static_cast<int>(2.*hydroXmax/hydroDx + 0.001);
        ietamax = static_cast<int>(2.*hydro_eta_max/hydroDeta + 0.001);

        // read in temperature, QGP fraction , flow velocity
        // The name of the evolution file: evolution_name
        string evolution_name = "results/evolution_all_xyeta.dat";
        cout << "Evolution file name = " << evolution_name << endl;
        std::FILE *fin;
        fin = std::fopen(evolution_name.c_str(), "rb");
        if (fin == NULL) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "Unable to open file: " << evolution_name << endl;
            exit(1);
        }
        
        int idx[4];
        int itau_max = 0;
        double ideal_variables[4];
        fluidCell_3D_new newCell;
        int ik = 0;
        while (true) {
            int status = 0;
            status = std::fread(&idx, sizeof(int), 4, fin);
            if (status == 0) break;
            
            status = std::fread(&ideal_variables, sizeof(double), 4, fin);
            if (status == 0) {
                cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                     << "file format is read in wrong" << endl;
                exit(1);
            }

            double muB_local = 0.0;
            if (turn_on_rhob == 1) {
                status = std::fread(&muB_local, sizeof(double), 1, fin);
                if (status == 0) {
                    cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                         << "file format is read in wrong" << endl;
                    exit(1);
                }
            }
            
            double Wmunu[5] = {0., 0., 0., 0., 0.};
            if (turn_on_shear == 1) {
                status = std::fread(&Wmunu, sizeof(double), 5, fin);
                if (status == 0) {
                    cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                         << "file format is read in wrong" << endl;
                    exit(1);
                }
            }

            double pi11 = Wmunu[0];
            double pi12 = Wmunu[1];
            double pi13 = Wmunu[2];
            double pi22 = Wmunu[3];
            double pi23 = Wmunu[4];
           
            double bulkPi;
            if (turn_on_bulk == 1) {
                status = std::fread(&bulkPi, sizeof(double), 1, fin);
                if (status == 0) {
                    cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                         << "file format is read in wrong" << endl;
                    exit(1);
                }
            }

            double qmu[3] = {0., 0., 0.};
            if (turn_on_diff == 1) {
                status = std::fread(&qmu, sizeof(double), 3, fin);
                if (status == 0) {
                    cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                         << "file format is read in wrong" << endl;
                    exit(1);
                }
            }

            if (itau_max < idx[0])
                itau_max = idx[0];
            newCell.itau = idx[0];
            newCell.ix = idx[1];
            newCell.iy = idx[2];
            newCell.ieta = idx[3];
            newCell.temperature = ideal_variables[0];
            newCell.ux = ideal_variables[1];
            newCell.uy = ideal_variables[2];
            newCell.ueta = ideal_variables[3];
            newCell.pi11 = pi11;
            newCell.pi12 = pi12;
            newCell.pi13 = pi13;
            newCell.pi22 = pi22;
            newCell.pi23 = pi23;
            newCell.bulkPi = bulkPi;
            lattice_3D_new->push_back(newCell);
            ik++;
            if (ik%50000 == 0)
                cout << "o" << flush;
        }
        cout << endl;
        std::fclose(fin);
        itaumax = itau_max;
        hydroTauMax = hydroTau0 + hydroDtau*itaumax;
    }

    // One final step for easy automation of MARTINI:
    // hydroTauMax is reset for the case where writing to evolution.dat
    // ended early (due to all cells freezing out):
    if (whichHydro == 6) {
        hydroTauMax = (
            hydroTau0 + hydroDtau*static_cast<int>(
                        static_cast<double>(lattice_3D->size())
                        /((2.*hydroXmax/hydroDx+1.)*(2.*hydroXmax/hydroDx+1.)
                        *2.*(hydro_eta_max/hydroDeta))));
        itaumax = static_cast<int>((hydroTauMax-hydroTau0)/hydroDtau+0.001);
    }
    if (whichHydro == 8 || whichHydro == 9) {
        hydroTauMax = (
            hydroTau0 + hydroDtau*static_cast<int>(
                        static_cast<double>(lattice_2D->size())
                        /((2.*hydroXmax/hydroDx)*(2.*hydroXmax/hydroDx)) - 1));
        itaumax = static_cast<int>((hydroTauMax - hydroTau0)/hydroDtau);
    }

    cout << "hydro_tau0 = " << hydroTau0 << " fm"<< endl;
    cout << "hydro_tau_max = " << hydroTauMax << " fm" << endl;
    cout << "hydry_dtau = " << hydroDtau << " fm" << endl;
    cout << "hydro_Xmax = " << hydroXmax << " fm" << endl;
    cout << "hydro_dx = " << hydroDx << " fm" << endl;
    cout << "hydro_eta_max = " << hydro_eta_max << " fm" << endl;
    cout << "hydro_deta = " << hydroDeta << " fm" << endl;
}

void Hydroinfo_MUSIC::get_hydro_cell_info_3d(int cell_id,
                                             fluidCell_3D_new *info) {
    info->itau = (*lattice_3D_new)[cell_id].itau;
    info->ix = (*lattice_3D_new)[cell_id].ix;
    info->iy = (*lattice_3D_new)[cell_id].iy;
    info->ieta = (*lattice_3D_new)[cell_id].ieta;
    info->temperature = (*lattice_3D_new)[cell_id].temperature;
    info->ux = (*lattice_3D_new)[cell_id].ux;
    info->uy = (*lattice_3D_new)[cell_id].uy;
    info->ueta = (*lattice_3D_new)[cell_id].ueta;
    info->pi11 = (*lattice_3D_new)[cell_id].pi11;
    info->pi12 = (*lattice_3D_new)[cell_id].pi12;
    info->pi13 = (*lattice_3D_new)[cell_id].pi13;
    info->pi22 = (*lattice_3D_new)[cell_id].pi22;
    info->pi23 = (*lattice_3D_new)[cell_id].pi23;
    info->bulkPi = (*lattice_3D_new)[cell_id].bulkPi;
}

void Hydroinfo_MUSIC::getHydroValues(double x, double y,
                                     double z, double t, fluidCell* info) {
// For interpolation of evolution files in tau-eta coordinates. Only the
// reading of MUSIC's evolution_xyeta.dat file is implemented here.
// For simplicity, hydro_eta_max refers to MUSIC's eta_size, and similarly for
// hydroDeta; however, x, y, z, and t are as usual to stay compatible with
// MARTINI.
    double tau, eta;
    if (use_tau_eta_coordinate == 1) {
        if (t*t > z*z) {
            tau = sqrt(t*t-z*z);
            eta = 0.5*log((t+z)/(t-z));
        } else {
            tau = 0.;
            eta = 0.;
        }
    } else {
        // if the medium is given in cartesian coordinates
        // set tau and eta to t and z
        tau = t;
        eta = z;
    }

    int ieta = floor((hydro_eta_max+eta)/hydroDeta + 0.0001);
    if (hydroWhichHydro == 8)
        ieta = 0;

    int itau = floor((tau-hydroTau0)/hydroDtau + 0.0001);
    int ix = floor((hydroXmax+x)/hydroDx + 0.0001);
    int iy = floor((hydroXmax+y)/hydroDx + 0.0001);

    double xfrac = (x - (static_cast<double>(ix)*hydroDx - hydroXmax))/hydroDx;
    double yfrac = (y - (static_cast<double>(iy)*hydroDx - hydroXmax))/hydroDx;
    double etafrac = (eta/hydroDeta - static_cast<double>(ieta)
                      + 0.5*static_cast<double>(ietamax));
    double taufrac = (tau - hydroTau0)/hydroDtau - static_cast<double>(itau);

    if (ix < 0 || ix >= ixmax) {
        cout << "[MARTINI:Hydroinfo_MUSIC::getHydroValues]: "
             << "WARNING - x out of range x=" << x
             << ", ix=" << ix << ", ixmax=" << ixmax << endl;
        cout << "x=" << x << " y=" << y << " eta=" << eta
             << " ix=" << ix << " iy=" << iy << " ieta=" << ieta << endl;
        cout << "t=" << t << " tau=" << tau
             << " itau=" << itau << " itaumax=" << itaumax << endl;

        info->temperature = 0.0;
        info->vx = 0.0;
        info->vy = 0.0;
        info->vz = 0.0;
        return;
    }
    if (iy < 0 || iy >= ixmax) {
        cout << "[MARTINI:Hydroinfo_MUSIC::getHydroValues]: "
             << "WARNING - y out of range, y=" << y << ", iy="  << iy
             << ", iymax=" << ixmax << endl;
        cout << "x=" << x << " y=" << y << " eta=" << eta
             << " ix=" << ix << " iy=" << iy << " ieta=" << ieta << endl;
        cout << "t=" << t << " tau=" << tau
             << " itau=" << itau << " itaumax=" << itaumax << endl;

        info->temperature = 0.0;
        info->vx = 0.0;
        info->vy = 0.0;
        info->vz = 0.0;
        return;
    }
    if (itau < 0 || itau > itaumax) {
        cout << "[MARTINI:Hydroinfo_MUSIC::getHydroValues]: WARNING - "
             << "tau out of range, itau=" << itau << ", itaumax=" << itaumax
             << endl;
        cout << "[MARTINI:Hydroinfo_MUSIC::getHydroValues]: tau= " << tau
             << ", hydroTauMax = " << hydroTauMax
             << ", hydroDtau = " << hydroDtau << endl;

        info->temperature = 0.0;
        info->vx = 0.0;
        info->vy = 0.0;
        info->vz = 0.0;
        return;
    }
    if (ieta < 0 || ieta >= ietamax) {
        cout << "[MARTINI:Hydroinfo_MUSIC::getHydroValues]: WARNING - "
             << "eta out of range, ieta=" << ieta << ", ietamax=" << ietamax
             << endl;
        info->temperature = 0.0;
        info->vx = 0.0;
        info->vy = 0.0;
        info->vz = 0.0;
        return;
    }

  // The array of positions on the 4-dimensional rectangle:
  int position[2][2][2][2];
  for (int ipx = 0; ipx < 2; ipx++) {
        int px;
        if (ipx == 0 || ix == ixmax-1)
            px = ix;
        else
            px = ix + 1;
        for (int ipy = 0; ipy < 2; ipy++) {
            int py;
            if (ipy == 0 || iy == ixmax-1)
                py = iy;
            else
                py = iy + 1;
            for (int ipeta = 0; ipeta < 2; ipeta++) {
                int peta;
                if (ipeta == 0 || ieta == ietamax-1)
                    peta = ieta;
                else
                    peta = ieta + 1;
                for (int iptau = 0; iptau < 2; iptau++) {
                    int ptau;
                    if (iptau == 0 || itau == itaumax-1)
                        ptau = itau;
                    else
                        ptau = itau + 1;
                    position[ipx][ipy][ipeta][iptau] = (
                                px + ixmax*(py + ixmax*(peta + ietamax*ptau)));
                }
            }
        }
    }

    // And now, the interpolation:
    double T = 0.0;
    double vx = 0.0;
    double vy = 0.0;
    double vz = 0.0;
    double ux = 0.0;
    double uy = 0.0;
    double ueta = 0.0;
    double pi00 = 0.0;
    double pi01 = 0.0;
    double pi02 = 0.0;
    double pi03 = 0.0;
    double pi11 = 0.0;
    double pi12 = 0.0;
    double pi13 = 0.0;
    double pi22 = 0.0;
    double pi23 = 0.0;
    double pi33 = 0.0;
    double bulkPi = 0.0;

    fluidCell_2D *HydroCell_2D_ptr1, *HydroCell_2D_ptr2;
    fluidCell_3D *HydroCell_3D_ptr1, *HydroCell_3D_ptr2;
    for (int iptau = 0; iptau < 2; iptau++) {
        double taufactor;
        if (iptau == 0)
            taufactor = 1. - taufrac;
        else
            taufactor = taufrac;
        for (int ipeta = 0; ipeta < 2; ipeta++) {
            double etafactor;
            if (ipeta == 0)
                etafactor = 1. - etafrac;
            else
                etafactor = etafrac;
            for (int ipy = 0; ipy < 2; ipy++) {
                double yfactor;
                if (ipy == 0)
                    yfactor = 1. - yfrac;
                else
                    yfactor = yfrac;

                double prefrac = yfactor*etafactor*taufactor;

                if (boost_invariant) {
                    HydroCell_2D_ptr1 = (
                            &(*lattice_2D)[position[0][ipy][ipeta][iptau]]);
                    HydroCell_2D_ptr2 = (
                            &(*lattice_2D)[position[1][ipy][ipeta][iptau]]);
                    T += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->temperature
                                  + xfrac*HydroCell_2D_ptr2->temperature);
                    ux += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->ux
                                    + xfrac*HydroCell_2D_ptr2->ux);
                    uy += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->uy
                                    + xfrac*HydroCell_2D_ptr2->uy);
                    ueta += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->ueta
                                    + xfrac*HydroCell_2D_ptr2->ueta);
                    pi00 += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->pi00
                                    + xfrac*HydroCell_2D_ptr2->pi00);
                    pi01 += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->pi01
                                    + xfrac*HydroCell_2D_ptr2->pi01);
                    pi02 += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->pi02
                                    + xfrac*HydroCell_2D_ptr2->pi02);
                    pi11 += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->pi11
                                    + xfrac*HydroCell_2D_ptr2->pi11);
                    pi12 += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->pi12
                                    + xfrac*HydroCell_2D_ptr2->pi12);
                    pi22 += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->pi22
                                    + xfrac*HydroCell_2D_ptr2->pi22);
                    pi33 += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->pi33
                                    + xfrac*HydroCell_2D_ptr2->pi33);
                    bulkPi += prefrac*((1. - xfrac)*HydroCell_2D_ptr1->bulkPi
                                    + xfrac*HydroCell_2D_ptr2->bulkPi);
                } else {
                    HydroCell_3D_ptr1 = (
                            &(*lattice_3D)[position[0][ipy][ipeta][iptau]]);
                    HydroCell_3D_ptr2 = (
                            &(*lattice_3D)[position[1][ipy][ipeta][iptau]]);
                    T += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->temperature
                                  + xfrac*HydroCell_3D_ptr2->temperature);
                    vx += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->vx
                                    + xfrac*HydroCell_3D_ptr2->vx);
                    vy += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->vy
                                    + xfrac*HydroCell_3D_ptr2->vy);
                    vz += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->vz
                                    + xfrac*HydroCell_3D_ptr2->vz);
                    pi00 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi00
                                    + xfrac*HydroCell_3D_ptr2->pi00);
                    pi01 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi01
                                    + xfrac*HydroCell_3D_ptr2->pi01);
                    pi02 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi02
                                    + xfrac*HydroCell_3D_ptr2->pi02);
                    pi03 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi03
                                    + xfrac*HydroCell_3D_ptr2->pi03);
                    pi11 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi11
                                    + xfrac*HydroCell_3D_ptr2->pi11);
                    pi12 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi12
                                    + xfrac*HydroCell_3D_ptr2->pi12);
                    pi13 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi13
                                    + xfrac*HydroCell_3D_ptr2->pi13);
                    pi22 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi22
                                    + xfrac*HydroCell_3D_ptr2->pi22);
                    pi23 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi23
                                    + xfrac*HydroCell_3D_ptr2->pi23);
                    pi33 += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->pi33
                                    + xfrac*HydroCell_3D_ptr2->pi33);
                    bulkPi += prefrac*((1. - xfrac)*HydroCell_3D_ptr1->bulkPi
                                    + xfrac*HydroCell_3D_ptr2->bulkPi);
                }
            }
        }
    }

    if (boost_invariant) {      // for boost invariant medium
        double eta_local = 0.5*log((t + z)/(t - z));
        double sinh_eta, cosh_eta;
        if (fabs(eta_local) < 1e-6) {
            // use Taylor expansion for small eta_s to speed up
            // avoiding to evaluate sinh and cosh
            sinh_eta = eta_local;
            cosh_eta = 1.0 + 0.5*eta_local*eta_local;
        } else {
            sinh_eta = sinh(eta_local);
            cosh_eta = cosh(eta_local);
        }
        double utau = sqrt(1. + ux*ux + uy*uy + ueta*ueta);
        double uz = utau*sinh_eta + ueta*cosh_eta;
        double ut = utau*cosh_eta + ueta*sinh_eta;
        vx = ux/ut;
        vy = uy/ut;
        vz = uz/ut;
    }

    info->temperature = T;
    info->vx = vx;
    info->vy = vy;
    info->vz = vz;

    info->ed = 1.0;                 // pi's are already divided by e+P
    info->sd = 0.0;
    info->pressure = 0.0;

    info->pi[0][0] = pi00;
    info->pi[0][1] = pi01;
    info->pi[0][2] = pi02;
    info->pi[0][3] = pi03;
    info->pi[1][0] = pi01;
    info->pi[1][1] = pi11;
    info->pi[1][2] = pi12;
    info->pi[1][3] = pi13;
    info->pi[2][0] = pi02;
    info->pi[2][1] = pi12;
    info->pi[2][2] = pi22;
    info->pi[2][3] = pi23;
    info->pi[3][0] = pi03;
    info->pi[3][1] = pi13;
    info->pi[3][2] = pi23;
    info->pi[3][3] = pi33;

    info->bulkPi = bulkPi;
    return;
}

void Hydroinfo_MUSIC::output_temperature_evolution(string filename_base) {
    fluidCell *hydroInfo = new fluidCell;
    for (int i = 0; i < itaumax; i++) {
        double tau = hydroTau0 + i*hydroDtau;
        ostringstream filename;
        filename << filename_base << "_tau_" << tau << ".dat";
        ofstream temp_evo(filename.str().c_str());
        for (int ix = 0; ix < ixmax; ix++) {
            double x_local = -hydroXmax + ix*hydroDx;
            for (int iy = 0; iy < ixmax; iy++) {
                double y_local = -hydroXmax + iy*hydroDx;
                getHydroValues(x_local, y_local, 0.0, tau, hydroInfo);
                double temp_local = hydroInfo->temperature;
                temp_evo << scientific << setw(16) << setprecision(8)
                         << temp_local << "   ";
            }
            temp_evo << endl;
        }
        temp_evo.close();
    }
    delete hydroInfo;
}

void Hydroinfo_MUSIC::update_grid_info(
    double tau0, double tau_max, double dtau,
    double x_max, double dx, double eta_max, double deta) {
    hydroTau0 = tau0;
    hydroTauMax = tau_max;
    hydroDtau = dtau;
    hydroXmax = x_max;
    hydroDx = dx;
    hydro_eta_max = eta_max;
    hydroDeta = deta;
    if (hydroWhichHydro == 8) {
        itaumax = static_cast<int>((tau_max-tau0)/dtau+0.001);
        ixmax = static_cast<int>(2*x_max/dx+0.001);
        ietamax = static_cast<int>(2*eta_max/deta+0.001);
    }
    if (hydroWhichHydro == 6) {
        itaumax = static_cast<int>((tau_max-tau0)/dtau+0.001);
        ixmax = static_cast<int>(2*x_max/dx+0.001);
        ietamax = static_cast<int>(2*eta_max/deta+0.001);
    }
}
