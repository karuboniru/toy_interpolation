#include "hkkm_reader.hxx"
#include "hkkm_reader_3D.hxx"
#include "interpolation.hxx"
#include "spline.hxx"

#include <TSystem.h>

#include <boost/python.hpp>
#include <boost/python/init.hpp>
#include <boost/python/numpy.hpp>
#include <boost/python/numpy/ndarray.hpp>
#include <boost/python/overloads.hpp>

#include <array>
#include <cmath>
#include <format>
#include <ostream>
#include <print>

namespace np = boost::python::numpy;

namespace {
const bool initializer = []() {
  Py_Initialize();
  np::initialize();
  gSystem->ResetSignal(kSigBus);
  gSystem->ResetSignal(kSigSegmentationViolation);
  gSystem->ResetSignal(kSigIllegalInstruction);
  return true;
}();

constexpr size_t n_logE_points = 101;
constexpr size_t n_costh_bins = 20;
constexpr size_t n_costh_points = n_costh_bins + 1;
constexpr size_t n_phi_bins = 12;
constexpr size_t n_phi_points = n_phi_bins + 1;

constexpr axis_object logE_points{
    .min = -1, .max = 4, .n_points = n_logE_points};
constexpr axis_object costh_points{
    .min = -1, .max = 1, .n_points = n_costh_points};
constexpr axis_object phi_points{
    .min = 0, .max = M_PI * 2, .n_points = n_phi_points};

size_t pdg2idx(int pdg) {
  switch (pdg) {
  case 12:
    return 0;
  case 14:
    return 1;
  case -12:
    return 2;
  case -14:
    return 3;
  default:
    throw std::invalid_argument("Invalid PDG code");
  }
}
constexpr std::array<int, 4> pdg_list{12, 14, -12, -14};

class hkkm_2d {
public:
  hkkm_2d(const char *fluxFilePath) {
    HKKM_READER_2D reader(fluxFilePath);
    for (const auto pdg : pdg_list) {
      auto &flux_hist = reader[pdg];
      auto &interop_obj = interp[pdg2idx(pdg)];
      for (size_t i = 0; i < n_logE_points; ++i) {
        double cdf_along_costh = 0;
        interop_obj[{i, 0}] = 0;
        for (size_t j = 0; j < n_costh_bins; ++j) {
          cdf_along_costh += flux_hist.GetBinContent(i + 1, j + 1);
          interop_obj[{i, j + 1}] = cdf_along_costh * 0.1;
        }
      }
    }
  }
  [[nodiscard]] double get_flux(double E, double costh, int pdg) const {
    auto logE = std::log10(E);
    return interp[pdg2idx(pdg)].do_interpolation({logE, costh}, {false, true});
  }

  [[nodiscard]] np::ndarray get_flux(const np::ndarray &E,
                                     const np::ndarray &costh, int pdg) const {
    namespace np = np;
    auto shape = E.get_shape();
    auto dimension = E.get_nd();
    auto result = np::empty(dimension, shape, np::dtype::get_builtin<double>());
    auto E_data = reinterpret_cast<double *>(E.get_data());
    auto costh_data = reinterpret_cast<double *>(costh.get_data());
    auto result_data = reinterpret_cast<double *>(result.get_data());
    auto size = shape[0];
    for (auto nd = 1; nd < dimension; ++nd) {
      size *= shape[nd];
    }
#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < size; ++i) {
      result_data[i] = get_flux(E_data[i], costh_data[i], pdg);
    }
    return result;
  }

private:
  std::array<interpolate<3, 4>, 4> interp{
      interpolate<3, 4>{{logE_points, costh_points}},
      interpolate<3, 4>{{logE_points, costh_points}},
      interpolate<3, 4>{{logE_points, costh_points}},
      interpolate<3, 4>{{logE_points, costh_points}}};
};

class hkkm_3d {
public:
  hkkm_3d(const char *fluxFilePath) {
    HKKM_READER_3D reader(fluxFilePath);
    for (const auto pdg : pdg_list) {
      auto &flux_hist = reader[pdg];
      auto &interop_obj = interp[pdg2idx(pdg)];
      for (size_t i = 0; i < n_logE_points; ++i) {
        interop_obj[{i, 0, 0}] = 0;
        for (size_t j = 0; j < n_costh_bins; ++j) {
          interop_obj[{i, j + 1, 0}] = 0;
          for (size_t k = 0; k < n_phi_bins; ++k) {
            interop_obj[{i, 0, k + 1}] = 0;
            interop_obj[{i, j + 1, k + 1}] =
                flux_hist.GetBinContent(i + 1, j + 1, k + 1) * 0.1 *
                (M_PI * 2 / n_phi_bins);
            interop_obj[{i, j + 1, k + 1}] += interop_obj[{i, j + 1, k}];
          }
          for (size_t k = 0; k < n_phi_bins; ++k) {
            interop_obj[{i, j + 1, k + 1}] += interop_obj[{i, j, k + 1}];
          }
        }
      }
    }
  }
  [[nodiscard]] double get_flux(double E, double costh, double phi,
                                int pdg) const {
    auto logE = std::log10(E);
    return interp[pdg2idx(pdg)].do_interpolation(
        {logE, costh, phi}, {false, true, true}, {false, false, true});
  }

  // takes 3 numpy array and return a numpy array
  [[nodiscard]] np::ndarray get_flux(const np::ndarray &E,
                                     const np::ndarray &costh,
                                     const np::ndarray &phi, int pdg) const {
    namespace np = np;
    auto shape = E.get_shape();
    auto dimension = E.get_nd();
    auto result = np::empty(dimension, shape, np::dtype::get_builtin<double>());
    auto E_data = reinterpret_cast<double *>(E.get_data());
    auto costh_data = reinterpret_cast<double *>(costh.get_data());
    auto phi_data = reinterpret_cast<double *>(phi.get_data());
    auto result_data = reinterpret_cast<double *>(result.get_data());
    auto size = shape[0];
    for (auto nd = 1; nd < dimension; ++nd) {
      size *= shape[nd];
    }
#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < size; ++i) {
      result_data[i] = get_flux(E_data[i], costh_data[i], phi_data[i], pdg);
    }
    return result;
  }

private:
  std::array<interpolate<3, 4, 4>, 4> interp{
      interpolate<3, 4, 4>{{logE_points, costh_points, phi_points}},
      interpolate<3, 4, 4>{{logE_points, costh_points, phi_points}},
      interpolate<3, 4, 4>{{logE_points, costh_points, phi_points}},
      interpolate<3, 4, 4>{{logE_points, costh_points, phi_points}}};
};

class spline_reader_py : public spline_reader {
public:
  using spline_reader::spline_reader;
  [[nodiscard]] np::ndarray get_cross_section(int neutrino_id, int target_id,
                                              const np::ndarray &energy,
                                              const std::string &channel_name) {
    auto shape = energy.get_shape();
    auto dimension = energy.get_nd();
    auto result = np::empty(dimension, shape, np::dtype::get_builtin<double>());
    auto energy_data = reinterpret_cast<double *>(energy.get_data());
    auto result_data = reinterpret_cast<double *>(result.get_data());
    auto size = shape[0];
    for (auto nd = 1; nd < dimension; ++nd) {
      size *= shape[nd];
    }
    // just trigger the cache to avoid multithreading issue
    result_data[0] = spline_reader::get_cross_section(
        neutrino_id, target_id, energy_data[0], channel_name);
#pragma omp parallel for schedule(dynamic)
    for (size_t i = 1; i < size; ++i) {
      result_data[i] = spline_reader::get_cross_section(
          neutrino_id, target_id, energy_data[i], channel_name);
    }
    return result;
  }
};
} // namespace

BOOST_PYTHON_MODULE(hkkm_interpolation) {
  boost::python::class_<hkkm_2d>("hkkm_2d", boost::python::init<const char *>())
      .def("get_flux",
           static_cast<double (hkkm_2d::*)(double, double, int) const>(
               &hkkm_2d::get_flux))
      .def("get_flux",
           static_cast<np::ndarray (hkkm_2d::*)(
               const np::ndarray &, const np::ndarray &, int) const>(
               &hkkm_2d::get_flux));

  boost::python::class_<hkkm_3d>("hkkm_3d", boost::python::init<const char *>())
      .def("get_flux",
           static_cast<double (hkkm_3d::*)(double, double, double, int) const>(
               &hkkm_3d::get_flux))
      .def("get_flux",
           static_cast<np::ndarray (hkkm_3d::*)(
               const np::ndarray &, const np::ndarray &, const np::ndarray &,
               int) const>(&hkkm_3d::get_flux));

  boost::python::class_<spline_reader_py>("genie_spline",
                                          boost::python::init<const char *>())
      .def("get_cross_section", &spline_reader::get_cross_section)
      .def("get_cross_section", &spline_reader_py::get_cross_section);
}