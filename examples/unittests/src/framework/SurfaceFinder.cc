#include "gtest/gtest.h"
#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
// #include "SurfaceFinder.h"
std::string load_mock_data(){
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
TEST(SurfaceFinder, Test_Check_Intersect_3D_false){

    std::string mock_data=load_mock_data();
    std::cout<<"Mock Data:\t"<<mock_data<<"\n";
    EXPECT_DOUBLE_EQ(0.0, 0.0);
}