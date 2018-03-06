# JETSCAPE-COMP
New code repository of the JETSCAPE collaboration.
New C++ code will be stored in this repository and used for version control.

## How to contribute to the code

1. *fork* one copy of the repository to your github account. In order to do this, first ask Abhijit to give you the privilege.

2. Follow the working flows in this [link](https://guides.github.com/introduction/flow/)

3. Important notes:
>   (1). *Pull request* does not mean "fork", it means that you request Abhijit to pull new modifications from your commits.
>
>   (2). The reposity is private, so cloning the repository from your github account to local machine is not as simple as usual. First find the reposity address by clicking on the green button "Clone or download", copy the address, for me it is: https://github.com/lgpang/JETSCAPE-COMP.git . For ordinary git reposity, one can use:  git clone https://github.com/lgpang/JETSCAPE-COMP.git
>
>    But for this private reposity, you must use: git clone https://your_user_name@github.com/lgpang/JETSCAPE-COMP.git ,
>    where your_user_name is your user name of the github account.
>
>    (3) With this local copy you can create one branch first, and do modifications on that branch.
>
>    (4) After you are satisfied with your modifications, *push* the modification to github server
>
>    (5) Go to Abhijit's reposity and click "pull request" to ask for reviewing from other programers.
>
>    (6) After being reviewed, Abhijit can merge the modifications in that branch to the master branch.

## Coding standard
1. C++11
2. [CMake](https://cmake.org/cmake-tutorial/) for compiling
3. [gtest](https://github.com/google/googletest/blob/master/googletest/docs/Primer.md) for unittest


## Use google test and CMake
1. There is one main CMakeLists.txt for the whole jetscape project
2. There is one CMakeLists.txt in each directory where you store your src code or tests code.
3. For example, in src/fluid_dynamics/test/, there is one CMakeLists.txt with content,
```cmake
    ##############
    # Unit Tests
    ##############
    add_executable(runUnitTests linear_interpolation.cc)
    
    # Standard linking to gtest stuff.
    target_link_libraries(runUnitTests gtest gtest_main)
    
    # This is so you can do 'make test' to see all your tests run, instead of
    # manually running the executable runUnitTests to see those specific tests.
    add_test(NAME test1 COMMAND runUnitTests)
    
    # You can also omit NAME and COMMAND. The second argument could be some other
    # test executable.
    add_test(test2 runUnitTests)  
```
4. You can run all the unit tests by,
```bash
   mkdir build
   cd build
   cmake .. -Dtest=ON
   make
   make test
```
The results would be,
```bash
    Running tests...
    Test project /Users/lgpang/Jetscape/JETSCAPE-COMP/build
    Start 1: test1
    1/2 Test #1: test1 ............................   Passed    0.01 sec
        Start 2: test2
        2/2 Test #2: test2 ............................   Passed    0.01 sec

        100% tests passed, 0 tests failed out of 2
```
        Total Test time (real) =   0.02 sec
