#include "hkkm_reader.hxx"
#include "hkkm_reader_3D.hxx"
#include "interpolation.hxx"
#include "spline.hxx"

#include <array>
#include <boost/python.hpp>
#include <boost/python/init.hpp>
#include <cmath>

namespace {
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

private:
  std::array<interpolate<3, 4, 4>, 4> interp{
      interpolate<3, 4, 4>{{logE_points, costh_points, phi_points}},
      interpolate<3, 4, 4>{{logE_points, costh_points, phi_points}},
      interpolate<3, 4, 4>{{logE_points, costh_points, phi_points}},
      interpolate<3, 4, 4>{{logE_points, costh_points, phi_points}}};
};
} // namespace

BOOST_PYTHON_MODULE(hkkm_interpolation) {
  boost::python::class_<hkkm_2d>("hkkm_2d", boost::python::init<const char *>())
      .def("get_flux", &hkkm_2d::get_flux);
  boost::python::class_<hkkm_3d>("hkkm_3d", boost::python::init<const char *>())
      .def("get_flux", &hkkm_3d::get_flux);
  boost::python::class_<spline_reader>("genie_spline",
                                       boost::python::init<const char *>())
      .def("get_cross_section", &spline_reader::get_cross_section);
}