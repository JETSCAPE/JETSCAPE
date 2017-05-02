// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "util.h"

VarMap make_var_map(std::map<std::string, boost::any>&& args) {
  VarMap var_map{};
  for (auto&& a : args)
    var_map.emplace(a.first, po::variable_value{a.second, false});
  return var_map;
}
