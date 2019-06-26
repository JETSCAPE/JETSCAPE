// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "../src/output.h"

#include "catch.hpp"
#include "util.h"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include "../src/event.h"
#include "../src/hdf5_utils.h"
#include "../src/nucleus.h"

using namespace trento;

TEST_CASE( "output" ) {
  auto var_map = make_var_map({
    {"normalization", 1.},
    {"reduced-thickness", 0.},
    {"grid-max", 9.},
    {"grid-step", 0.3},
    {"fluctuation", 1.},
    {"cross-section", 6.4},
    {"nucleon-width", 0.5}
  });

  // create a test event
  Event event{var_map};
  NucleonProfile profile{var_map};

  auto nucleusA = Nucleus::create("Pb");
  auto nucleusB = Nucleus::create("Pb");

  auto b = 4.*std::sqrt(random::canonical<>());
  nucleusA->sample_nucleons(+.5*b);
  nucleusB->sample_nucleons(-.5*b);

  for (auto&& A : *nucleusA)
    for (auto&& B : *nucleusB)
      profile.participate(A, B);

  event.compute(*nucleusA, *nucleusB, profile);

  SECTION( "no output" ) {
    capture_stdout capture;
    Output output{make_var_map({{"quiet", true}, {"number-events", 1}})};
    output(1, b, event);
    CHECK( capture.stream.str().empty() );  // stdout should be empty
  }

  SECTION( "stdout only" ) {
    // write event 0 to stdout for the given number of events
    auto first_line = [&b, &event](int nev) {
      capture_stdout capture;
      Output output{make_var_map({{"quiet", false}, {"number-events", nev}})};
      output(0, b, event);
      std::string line;
      std::getline(capture.stream, line);
      return line;
    };

    // verify event number padding
    CHECK( first_line(1).substr(0, 1) == "0" );
    CHECK( first_line(10).substr(0, 1) == "0" );
    CHECK( first_line(11).substr(0, 2) == " 0" );
    CHECK( first_line(100).substr(0, 2) == " 0" );
    CHECK( first_line(101).substr(0, 3) == "  0" );

    // output lines for different total event numbers should be identical except
    // for the padding
    CHECK( first_line(1).substr(1) == first_line(1000).substr(3) );

    // read an output line back into separate objects
    int num, npart;
    double impact, mult, e2, e3, e4, e5;
    char end;
    {
      capture_stdout capture;
      Output output{make_var_map({{"quiet", false}, {"number-events", 1}})};
      output(0, b, event);
      capture.stream >> num >> impact >> npart >> mult >> e2 >> e3 >> e4 >> e5 >> std::ws;
      end = capture.stream.get();
    }

    // verify output data is correct
    CHECK( num == 0 );
    CHECK( impact == Approx(b) );
    CHECK( npart == event.npart() );
    CHECK( mult == Approx(event.multiplicity()) );
    CHECK( e2 == Approx(event.eccentricity().at(2)) );
    CHECK( e3 == Approx(event.eccentricity().at(3)) );
    CHECK( e4 == Approx(event.eccentricity().at(4)) );
    CHECK( e5 == Approx(event.eccentricity().at(5)) );

    // verify the end character is really the end of file
    CHECK( end == std::char_traits<char>::eof() );
  }

  SECTION( "text and stdout" ) {
    // configure for writing text files to a random temporary path
    temporary_path temp{};
    auto output_var_map = make_var_map({
      {"quiet", false},
      {"no-header", false},
      {"number-events", 50},
      {"output", temp.path}
    });
    Output output{output_var_map};

    // output directory should have been created
    CHECK( fs::exists(temp.path) );

    {
      // output two events and verify two lines were printed to stdout
      capture_stdout capture;
      output(3, b, event);
      output(27, b, event);
      int n = 0;
      std::string line;
      while (std::getline(capture.stream, line)) { ++n; }
      CHECK( n == 2 );
    }

    // verify event files were created
    CHECK( fs::exists(temp.path/"03.dat") );
    CHECK( fs::exists(temp.path/"27.dat") );

    {
      // read event file back in
      fs::ifstream ifs{temp.path/"03.dat"};
      std::string line;

      // check header
      std::getline(ifs, line);
      CHECK( line == "# event 3" );

      std::getline(ifs, line);
      CHECK( line.substr(0, 10) == "# b     = " );
      CHECK( std::stod(line.substr(10)) == Approx(b) );

      std::getline(ifs, line);
      CHECK( line.substr(0, 10) == "# npart = " );
      CHECK( std::stoi(line.substr(10)) == event.npart() );

      std::getline(ifs, line);
      CHECK( line.substr(0, 10) == "# mult  = " );
      CHECK( std::stod(line.substr(10)) == Approx(event.multiplicity()) );

      for (const auto& ecc : event.eccentricity()) {
        std::getline(ifs, line);
        CHECK( line.substr(0, 10) == ("# e" + std::to_string(ecc.first) + "    = ") );
        CHECK( std::stod(line.substr(10)) == Approx(ecc.second) );
      }

      // read the grid back in and check each element
      const auto* iter = event.reduced_thickness_grid().origin();
      double check;
      bool all_correct = false;
      while (ifs >> check)
        all_correct = (check == Approx(*(iter++))) || all_correct;
      CHECK( all_correct );

      // verify that all grid elements were checked
      const auto* grid_end = event.reduced_thickness_grid().origin() +
                             event.reduced_thickness_grid().num_elements();
      CHECK( iter == grid_end );
    }

    {
      // check event number header in second file
      fs::ifstream ifs{temp.path/"27.dat"};
      std::string line;
      std::getline(ifs, line);
      CHECK( line == "# event 27" );
    }

    // attempting to output to the same directory again should throw an error
    CHECK_THROWS_AS( Output{output_var_map}, std::runtime_error );
  }

#ifdef TRENTO_HDF5
  SECTION( "hdf5 only" ) {
    auto nev = 10;

    // configure for writing a random hdf5 file
    temporary_path temp{".hdf5"};

    auto output_var_map = make_var_map({
      {"quiet", true},
      {"number-events", nev},
      {"output", temp.path}
    });
    Output output{output_var_map};

    for (auto n = 0; n < nev; ++n)
      output(n, b, event);

    {
      H5::H5File file{temp.path.string(), H5F_ACC_RDONLY};
      CHECK( static_cast<int>(file.getNumObjs()) == nev );

      auto name = file.getObjnameByIdx(0);
      CHECK( name == "event_0" );

      auto dataset = file.openDataSet(name);

      // read back in the event grid to another array
      Event::Grid grid_check{event.reduced_thickness_grid()};
      dataset.read(grid_check.data(), H5::PredType::NATIVE_DOUBLE);

      // verify each grid element
      auto grid_correct = std::equal(
        grid_check.origin(),
        grid_check.origin() + grid_check.num_elements(),
        event.reduced_thickness_grid().origin(),
        [](const double& value_check, const double& value) {
          return value_check == Approx(value);
        }
      );
      CHECK( grid_correct );

      // verify attributes
      double double_check;
      int int_check;

      dataset.openAttribute("b").read(H5::PredType::NATIVE_DOUBLE, &double_check);
      CHECK( double_check == Approx(b) );

      dataset.openAttribute("npart").read(H5::PredType::NATIVE_INT, &int_check);
      CHECK( int_check == event.npart() );

      dataset.openAttribute("mult").read(H5::PredType::NATIVE_DOUBLE, &double_check);
      CHECK( double_check == Approx(event.multiplicity()) );

      for (const auto& ecc : event.eccentricity()) {
        dataset.openAttribute("e" + std::to_string(ecc.first))
          .read(H5::PredType::NATIVE_DOUBLE, &double_check);
        CHECK( double_check == Approx(ecc.second) );
      }

#if H5_VERSION_GE(1, 8, 14)  // causes memory leak on earlier versions
      CHECK( dataset.getNumAttrs() == 7 );
#endif
    }

    // attempting to output to the same file again should throw an error
    CHECK_THROWS_AS( Output{output_var_map}, std::runtime_error );

    // create another empty temporary file
    temporary_path temp2{".hdf5"};
    {
      fs::ofstream{temp2.path};
    }

    // should be able to write to an existing but empty file
    Output{
      make_var_map({
        {"quiet", true},
        {"number-events", 1},
        {"output", temp2.path}
      })
    }(0, b, event);

    CHECK( fs::file_size(temp2.path) > 0);
  }
#endif  // TRENTO_HDF5
}
