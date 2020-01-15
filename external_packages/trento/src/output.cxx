// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// TRENTO3D: Three-dimensional extension of TRENTO by Weiyao Ke
// MIT License

#include "output.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/program_options/variables_map.hpp>

#include "event.h"
#include "hdf5_utils.h"

namespace trento {

namespace {

// These output functions are invoked by the Output class.

void write_stream(std::ostream& os, int width,
    int num, double impact_param, const Event& event) {
  using std::fixed;
  using std::setprecision;
  using std::setw;
  using std::scientific;

  // Write a nicely-formatted line of event properties.
  os << setprecision(10)
     << setw(width)            << num
     << setw(15) << fixed      << impact_param
     << setw(5)                << event.npart();
  if (event.with_ncoll()) os << setw(8)                << event.ncoll();
  os   << setw(18) << scientific << event.multiplicity()            
     << fixed;

  for (const auto& ecc : event.eccentricity())
    os << setw(14) << ecc.second;

  //for (const auto& psi : event.participant_plane())
  //  os << setw(14) << psi.second;

  os << '\n';
}

void write_text_file(const fs::path& output_dir, int width,
    int num, double impact_param, const Event& event, bool header) {
  // Open a numbered file in the output directory.
  // Pad the filename with zeros.
  std::ostringstream padded_fname{};
  padded_fname << std::setw(width) << std::setfill('0') << num << ".dat";
  fs::ofstream ofs{output_dir / padded_fname.str()};

  if (header) {
    // Write a commented header of event properties as key = value pairs.
    ofs << std::setprecision(10)
        << "# event "   << num                  << '\n'
        << "# b     = " << impact_param         << '\n'
        << "# npart = " << event.npart()        << '\n'
        << "# mult  = " << event.multiplicity() << '\n';

    for (const auto& ecc : event.eccentricity())
      ofs << "# e" << ecc.first << "    = " << ecc.second << '\n';

    for (const auto& psi : event.participant_plane())
      ofs << "# psi" << psi.first << "    = " << psi.second << '\n';
  }

  // Write IC profile as a block grid.  Use C++ default float format (not
  // fixed-width) so that trailing zeros are omitted.  This significantly
  // increases output speed and saves disk space since many grid elements are
  // zero.
  bool is3d;
  if (event.density_grid().shape()[2] == 1) is3d = false;
  else is3d = true;


  for (const auto& slice : event.density_grid()) {
    for (const auto& row : slice) {
      for (const auto& item : row) {
		 ofs << item << " ";
	  }
	  if (is3d) ofs << std::endl;
    }
	if (!is3d) ofs << std::endl;
  }
}

#ifdef TRENTO_HDF5

/// Simple functor to write many events to an HDF5 file.
class HDF5Writer {
 public:
  /// Prepare an HDF5 file for writing.
  HDF5Writer(const fs::path& filename);

  /// Write an event.
  void operator()(int num, double impact_param, const Event& event) const;

 private:
  /// Internal storage of the file object.
  H5::H5File file_;
};

// Add a simple scalar attribute to an HDF5 dataset.
//template <typename T>
//void hdf5_add_scalar_attr(
//    const H5::DataSet& dataset, const std::string& name, const T& value) {
//  const auto& datatype = hdf5::type<T>();
//  auto attr = dataset.createAttribute(name, datatype, H5::DataSpace{});
//  attr.write(datatype, &value);
//}
template <typename T>
void hdf5_add_scalar_attr(
    const H5::Group& group, const std::string& name, const T& value) {
  const auto& datatype = hdf5::type<T>();
  auto attr = group.createAttribute(name, datatype, H5::DataSpace{});
  attr.write(datatype, &value);
}

HDF5Writer::HDF5Writer(const fs::path& filename)
    : file_(filename.string(), H5F_ACC_TRUNC)
{}

void HDF5Writer::operator()(
    int num, double impact_param, const Event& event) const {
  const auto& grid1 = event.density_grid();
  const auto& grid2 = event.TAB_grid();

  // The dataset name is a prefix plus the event number.
  const std::string gp_name{"/event_" + std::to_string(num)};
  const std::string sd_name{gp_name + "/matter_density"};
  const std::string tab_name{gp_name + "/Ncoll_density"};

  // create a group called event_i
  auto group = H5::Group(file_.createGroup(gp_name));
  // Write event attributes.
  hdf5_add_scalar_attr(group, "b", impact_param);
  hdf5_add_scalar_attr(group, "npart", event.npart());
  hdf5_add_scalar_attr(group, "ncoll", event.ncoll());
  hdf5_add_scalar_attr(group, "mult", event.multiplicity());
  hdf5_add_scalar_attr(group, "dxy", event.dxy());
  hdf5_add_scalar_attr(group, "deta", event.deta());
  hdf5_add_scalar_attr(group, "Ny", grid1.shape()[0]);
  hdf5_add_scalar_attr(group, "Nx", grid1.shape()[1]);
  hdf5_add_scalar_attr(group, "Nz", grid1.shape()[2]);
  for (const auto& ecc : event.eccentricity())
    hdf5_add_scalar_attr(group, "e" + std::to_string(ecc.first), ecc.second);
  for (const auto& psi : event.participant_plane())
    hdf5_add_scalar_attr(group, "psi" + std::to_string(psi.first), psi.second);

  ////////////////////////////////////////////////////////////////////
  // Define HDF5 datatype and dataspace to match the density (3D) grid.
  
  const auto& datatype1 = hdf5::type<Event::Grid3D::element>();
  std::array<hsize_t, Event::Grid3D::dimensionality> shape1;
  std::copy(grid1.shape(), grid1.shape() + shape1.size(), shape1.begin());
  auto dataspace1 = hdf5::make_dataspace(shape1);

  // Set dataset storage properties.
  H5::DSetCreatPropList proplist1{};
  // Set chunk size to the entire grid.  For typical grid sizes (~100x100), this
  // works out to ~80 KiB, which is pretty optimal.  Anyway, it makes logical
  // sense to chunk this way, since chunks must be read contiguously and there's
  // no reason to read a partial grid.
  proplist1.setChunk(shape1.size(), shape1.data());
  // Set gzip compression level.  4 is the default in h5py.
  proplist1.setDeflate(4);

  // Create the new dataset and write the grid for matter density
  auto dataset1 = file_.createDataSet(sd_name, datatype1, dataspace1, proplist1);
  dataset1.write(grid1.data(), datatype1);

  //////////////////////////////////////////////////////////////////
  // Define HDF5 datatype and dataspace to match the Ncoll (2D) grid.
  const auto& datatype2 = hdf5::type<Event::Grid::element>();
  std::array<hsize_t, Event::Grid::dimensionality> shape2;
  std::copy(grid2.shape(), grid2.shape() + shape2.size(), shape2.begin());
  auto dataspace2 = hdf5::make_dataspace(shape2);
  
  // Set dataset storage properties.
  H5::DSetCreatPropList proplist2{};
  // Set chunk size to the entire grid.  For typical grid sizes (~100x100), this
  // works out to ~80 KiB, which is pretty optimal.  Anyway, it makes logical
  // sense to chunk this way, since chunks must be read contiguously and there's
  // no reason to read a partial grid.
  proplist1.setChunk(shape2.size(), shape2.data());
  // Set gzip compression level.  4 is the default in h5py.
  proplist1.setDeflate(4);

  // Create the new dataset and write the grid for matter density
  auto dataset2 = file_.createDataSet(tab_name, datatype2, dataspace2, proplist2);
  dataset2.write(grid2.data(), datatype2);
}

#endif  // TRENTO_HDF5

}  // unnamed namespace

Output::Output(const VarMap& var_map){
  // Determine the required width (padding) of the event number.  For example if
  // there are 10 events, the numbers are 0-9 and so no padding is necessary.
  // However given 11 events, the numbers are 00-10 with padded 00, 01, ...
  auto nevents = var_map["number-events"].as<int>();
  auto width = static_cast<int>(std::ceil(std::log10(nevents)));

  // Write to stdout unless the quiet option was specified.
  if (!var_map["quiet"].as<bool>()) {
    writers_.emplace_back(
      [width](int num, double impact_param, const Event& event) {
        write_stream(std::cout, width, num, impact_param, event);
      }
    );
  }

  // Possibly write to text or HDF5 files.
  if (var_map.count("output")) {
    const auto& output_path = var_map["output"].as<fs::path>();
    if (hdf5::filename_is_hdf5(output_path)) {
#ifdef TRENTO_HDF5
      if (fs::exists(output_path) && !fs::is_empty(output_path))
        throw std::runtime_error{"file '" + output_path.string() +
                                 "' exists, will not overwrite"};
      writers_.emplace_back(HDF5Writer{output_path});
#else
      throw std::runtime_error{"HDF5 output was not compiled"};
#endif  // TRENTO_HDF5
    } else {
      // Text files are all written into a single directory.  Require the
      // directory to be initially empty to avoid accidental overwriting and/or
      // spewing files into an already-used location.  If the directory does not
      // exist, create it.
      if (fs::exists(output_path)) {
        if (!fs::is_empty(output_path)) {
          throw std::runtime_error{"output directory '" + output_path.string() +
                                   "' must be empty"};
        }
      } else {
        fs::create_directories(output_path);
      }
      auto header = !var_map["no-header"].as<bool>();
      writers_.emplace_back(
        [output_path, width, header](
            int num, double impact_param, const Event& event) {
          write_text_file(output_path, width, num, impact_param, event, header);
        }
      );
    }
  }
}

}  // namespace trento
