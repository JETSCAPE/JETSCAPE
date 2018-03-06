// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <sstream>

#include <boost/filesystem.hpp>
#include <boost/program_options/variables_map.hpp>

#include "../src/fwd_decl.h"

// testing utilities

// Factory function; create a dummy boost::program_options::variables_map.
VarMap make_var_map(std::map<std::string, boost::any>&& args);

// redirect stdout to a stringstream and safely restore upon destruction
struct capture_stdout {
  capture_stdout() {
    // replace stdout buffer
    cout_orig = std::cout.rdbuf(stream.rdbuf());
  }

  ~capture_stdout() {
    // restore stdout to working state
    std::cout.rdbuf(cout_orig);
  }

  std::streambuf* cout_orig;
  std::stringstream stream;
};

// path that deletes itself when it goes out of scope
struct temporary_path {
  temporary_path(const fs::path& ext = fs::path{})
      : path(
          fs::temp_directory_path() /
          fs::unique_path().replace_extension(ext)
        )
    {}
  ~temporary_path() {
    fs::remove_all(path);
  }
  const fs::path path;
};

#endif  // UTIL_H
