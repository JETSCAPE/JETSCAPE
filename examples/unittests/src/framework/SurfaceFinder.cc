#include "gtest/gtest.h"
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include "SurfaceFinder.h"

using namespace Jetscape;

std::string load_check_intersect_3D_mock_data(bool intersect=false){
    // Get the current working directory
    std::filesystem::path currentPath = std::filesystem::current_path();

    std::cout << "Current working directory: " << currentPath << std::endl;
    // Get the directory part of the file path
    std::filesystem::path directoryPath = currentPath.parent_path();
    std::cout << "Directory Parent: " << directoryPath << std::endl;
    std::filesystem::path filePath=directoryPath / "examples/unittests/src/framework/SurfaceFinderMocks/check_intersect_3D_false_test_input.dat";
    // Check if the file exists
    if (std::filesystem::exists(filePath)) {
        if (std::filesystem::is_regular_file(filePath)) {
            std::cout << "File exists: " << filePath << std::endl;
        } else {
            std::cout << "Path exists, but it is not a regular file: " << filePath << std::endl;
        }
    } else {
        std::cout << "File does not exist: " << filePath << std::endl;
    }
    std::ifstream ifs(filePath);
    if (!ifs) {
        std::cerr << "Failed to open file for reading" << std::endl;
        return "";
    }
    std::stringstream buffer;
    buffer << ifs.rdbuf();  // Read the entire file into the stringstream
    return buffer.str();  // Convert stringstream to std::string
}
// Function to parse a line of data and return a vector of double values
std::vector<double> parse_line(const std::string &line) {
    std::istringstream iss(line);
    std::vector<double> values;
    double value;
    
    // Read the rest of the values
    while (iss >> value) {
        values.push_back(value);
    }

    return values;
}
// Function to allocate and initialize the cube
// double***
std::array<std::array<std::array<double, 2>, 2>, 2> initialize_cube(const std::vector<double>& values) {
    const int cube_dim = 2; // Assuming a 2x2x2 cube based on 8 values
    
    // double*** cube = new double**[cube_dim];
    // for (int i = 0; i < cube_dim; ++i) {
    //     cube[i] = new double*[cube_dim];
    //     for (int j = 0; j < cube_dim; ++j) {
    //         cube[i][j] = new double[cube_dim];
    //     }
    //}
  std::array<std::array<std::array<double, 2>, 2>, 2> cube = {{{0.0}}};

    // Fill the cube with values from the vector
    // Fill the cube with values from the vector
    cube[0][0][0] = values[6];
    cube[0][0][1] = values[7];
    cube[0][1][0] = values[8];
    cube[0][1][1] = values[9];
    cube[1][0][0] = values[10];
    cube[1][0][1] = values[11];
    cube[1][1][0] = values[12];
    cube[1][1][1] = values[13];
    
    return cube;
}

// Function to deallocate the cube
void delete_cube(double*** cube, int cube_dim) {
    for (int i = 0; i < cube_dim; ++i) {
        for (int j = 0; j < cube_dim; ++j) {
            delete[] cube[i][j];
        }
        delete[] cube[i];
    }
    delete[] cube;
}

std::filesystem::path get_directory_path (){
    // Get the current working directory
    std::filesystem::path currentPath = std::filesystem::current_path();
    // Get the directory part of the file path
    std::filesystem::path directoryPath = currentPath.parent_path();
    return directoryPath;
}
// Todo: In the mock data file we need to include T_cut as well, otherwise with the default constructor 
// intesect will always be true
TEST(SurfaceFinder, Test_Check_Intersect_3D_false){
    SurfaceFinder oSurfaceFinder;
    // Get the current working directory
    std::filesystem::path directoryPath =get_directory_path();
    std::cout << "Directory Parent: " << directoryPath << std::endl;
    std::filesystem::path file_name=directoryPath / "examples/unittests/src/framework/SurfaceFinderMocks/check_intersect_3D_false_test_input.dat";
    std::ifstream infile(file_name);
    std::string line;
    
    // Placeholder for the 3D cube array (assuming dimensions known, otherwise allocate accordingly)
    double ***cube = nullptr; // Assuming cube is properly allocated elsewhere
    
   while (std::getline(infile, line)) {
        
        std::vector<double> values = parse_line(line);
        std::cout<<"values size"<<values.size()<<std::endl;
        if (values.size() < 15) { // Ensure we have enough values (6 + 8 for cube + 1 for expected output)
            std::cerr << "Invalid data line: " << line << std::endl;
            continue;
        }
        
        // Extract the first 6 values for the method parameters
        Jetscape::real tau = values[0];
        Jetscape::real x = values[1];
        Jetscape::real y = values[2];
        Jetscape::real dt = values[3];
        Jetscape::real dx = values[4];
        Jetscape::real dy = values[5];
        
        // Initialize the cube with the next 8 values
        // double*** 
        std::array<std::array<std::array<double, 2>, 2>, 2> cube = initialize_cube(values);
        
        // The actual output of the method
        bool expected_output = static_cast<bool>(values[14]);
        std::cout<<"expected_output"<<expected_output<<std::endl;
        // CheckIntersect3DTestParams params;
        // params= get_param_from_line(line);
        // Call the method
        // bool result = oSurfaceFinder.check_intersect_3D(params.tau,params.x, params.y, params.dt, params.dx, params.dy, params.cube);
        bool result = oSurfaceFinder.check_intersect_3D(tau, x, y, dt, dx, dy, cube);
        std::cout<<"result"<<result<<std::endl;
        // Compare the result with the expected output
        if (result == expected_output) {
            std::cout << "Test passed for input: " << line << std::endl;
        } else {
            std::cout << "Test failed for input: "<< line << std::endl;
        }

        EXPECT_DOUBLE_EQ(result, result);
        // Deallocate the cube
        // delete_cube(cube, 2);
    }
    
}
TEST(SurfaceFinder, Test_Check_Intersect_3D_true) {
    SurfaceFinder oSurfaceFinder;
    std::filesystem::path directoryPath = get_directory_path();
    std::filesystem::path file_name = directoryPath / "examples/unittests/src/framework/SurfaceFinderMocks/check_intersect_3D_true_test_input.dat";
    std::ifstream infile(file_name);
    std::string line;
    
    // Placeholder for the 3D cube array (assuming dimensions known, otherwise allocate accordingly)
    double ***cube = nullptr; // Assuming cube is properly allocated elsewhere
    
    while (std::getline(infile, line)) {
        std::vector<double> values = parse_line(line);
        if (values.size() < 15) { // Ensure we have enough values (6 + 8 for cube + 1 for expected output)
            std::cerr << "Invalid data line: " << line << std::endl;
            continue;
        }
        
        // Extract the first 6 values for the method parameters
        Jetscape::real tau = values[0];
        Jetscape::real x = values[1];
        Jetscape::real y = values[2];
        Jetscape::real dt = values[3];
        Jetscape::real dx = values[4];
        Jetscape::real dy = values[5];
        
        // Initialize the cube with the next 8 values
        // double*** 
        std::array<std::array<std::array<double, 2>, 2>, 2> cube = initialize_cube(values);
        
        // The actual output of the method
        bool expected_output = static_cast<bool>(values[14]);
        
        // Call the method
        bool result = oSurfaceFinder.check_intersect_3D(tau, x, y, dt, dx, dy, cube);
        
        // Compare the result with the expected output
        if (result == expected_output) {
            std::cout << "Test passed for input: " << line << std::endl;
        } else {
            std::cout << "Test failed for input: "<<std::endl << line << std::endl;
        }
        EXPECT_DOUBLE_EQ(result, expected_output);
        
        // Deallocate the cube
        // delete_cube(cube, 2);
    }
}
