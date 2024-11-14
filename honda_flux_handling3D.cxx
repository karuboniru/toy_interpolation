#include "hkkm_reader_3D.hxx"
#include "interpolation.hxx"
#include <algorithm>
#include <boost/program_options.hpp>
#include <cmath>
#include <print>

#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TF2.h>
#include <TStyle.h>

int main(int argc, char **argv) {
  namespace po = boost::program_options;
  po::options_description desc("Allowed options");
  desc.add_options()("help", "produce help message")(
      "input-file", po::value<std::string>()->required(), "input file");
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  auto input_file = vm["input-file"].as<std::string>();

  HKKM_READER_3D reader(input_file);
  auto &numu = reader[14];
  constexpr size_t n_logE_points = 101;
  constexpr size_t n_costh_bins = 20;
  constexpr size_t n_costh_points = n_costh_bins + 1;
  constexpr size_t n_phi_bins = 12;
  constexpr size_t n_phi_points = n_phi_bins + 1;

  // here is where the interpolation method initializes
  axis_object logE_points{.min = -1, .max = 4, .n_points = n_logE_points};
  axis_object costh_points{.min = -1, .max = 1, .n_points = n_costh_points};
  axis_object phi_points{
      .min = 0, .max = 2 * TMath::Pi(), .n_points = n_phi_points};

  interpolate<2, 4, 4> interp{{logE_points, costh_points, phi_points}};
  for (size_t i = 0; i < n_logE_points; ++i) {
    interp[{i, 0, 0}] = 0;
    for (size_t j = 0; j < n_costh_bins; ++j) {
      interp[{i, j + 1, 0}] = 0;
      for (size_t k = 0; k < n_phi_bins; ++k) {
        interp[{i, 0, k + 1}] = 0;
        interp[{i, j + 1, k + 1}] = numu.GetBinContent(i + 1, j + 1, k + 1) *
                                    0.1 * (M_PI * 2 / n_phi_bins);
        interp[{i, j + 1, k + 1}] += interp[{i, j + 1, k}];
      }
      for (size_t k = 0; k < n_phi_bins; ++k) {
        interp[{i, j + 1, k + 1}] += interp[{i, j, k + 1}];
      }
    }
  }

  // then plot
  auto get_flux = [&](double E, double costh, double phi) {
    auto logE = std::log10(E);
    return interp.do_interpolation({logE, costh, phi}, {false, true, true});
  };
  auto get_flux_no_interop = [&](double E, double costh, double phi) {
    auto logE = std::log10(E);
    // return numu.Interpolate(logE, costh, phi);
    auto binx = numu.GetXaxis()->FindBin(logE);
    auto biny = numu.GetYaxis()->FindBin(costh);
    auto binz = numu.GetZaxis()->FindBin(phi);
    return numu.GetBinContent(binx, biny, binz);
  };

  TApplication app("app", nullptr, nullptr);
  TCanvas c1("c1", "c1", 800 * 3, 600);
  c1.Divide(3, 1);
  c1.cd(1);

  constexpr double fixed_costh = 0.85;
  constexpr double fixed_energy = 1.2;
  constexpr double fixed_phi = 1.2 * M_PI;

  TF1 flux_given_costh{
      "flux",
      [&](double *x, double *) {
        return std::log(get_flux(std::pow(10, x[0]), fixed_costh, fixed_phi));
      },
      -1, 4, 0};
  flux_given_costh.SetTitle(
      std::format(
          "Flux along cos#theta = {:.1f}, #phi = {:.1f}; log_{{10}}(E/GeV); ",
          fixed_costh, fixed_phi)
          .c_str());
  flux_given_costh.SetNpx(1000);
  const double costh_thres = 1.0;
  TF1 flux_given_E{
      "flux", [&](double *x, double *) { return get_flux(1, x[0], fixed_phi); },
      -costh_thres, costh_thres, 0};
  flux_given_E.SetNpx(1000);
  flux_given_E.SetTitle(
      std::format("Flux along E = {} GeV, #phi = {:.1f}; cos#theta",
                  fixed_energy, fixed_phi)
          .c_str());

  TF1 flux_given_costh_no_interop{
      "flux",
      [&](double *x, double *) {
        return std::log(
            get_flux_no_interop(std::pow(10, x[0]), fixed_costh, fixed_phi));
      },
      -1, 4, 0};
  flux_given_costh_no_interop.SetNpx(1000);
  TF1 flux_given_E_no_interop{"flux",
                              [&](double *x, double *) {
                                return get_flux_no_interop(1, x[0], fixed_phi);
                              },
                              -costh_thres, costh_thres, 0};
  flux_given_E_no_interop.SetNpx(1000);

  TF1 flux_given_E_costh{"flux",
                         [&](double *x, double *) {
                           return std::log(get_flux(1, fixed_costh, x[0]));
                         },
                         0, 2 * M_PI, 0};
  flux_given_E_costh.SetNpx(1000);
  flux_given_E_costh.SetTitle(
      std::format("Flux along E = 1 GeV, cos#theta = {:.1f}; phi", fixed_costh)
          .c_str());
  TF1 flux_given_E_costh_no_interop{
      "flux",
      [&](double *x, double *) {
        return std::log(get_flux_no_interop(1, fixed_costh, x[0]));
      },
      0, 2 * M_PI, 0};
  flux_given_E_costh_no_interop.SetNpx(1000);

  flux_given_costh.SetMaximum(std::max(
      flux_given_costh.GetMaximum(), flux_given_costh_no_interop.GetMaximum()));
  flux_given_costh.Draw();
  flux_given_costh_no_interop.SetLineColor(kBlue);
  flux_given_costh_no_interop.SetLineStyle(kDashed);
  flux_given_costh_no_interop.SetLineWidth(4);
  flux_given_costh_no_interop.Draw("same");

  c1.cd(2);
  // flux_given_E.SetMinimum(0);
  flux_given_E.SetMaximum(1.1 * std::max(flux_given_E.GetMaximum(),
                                         flux_given_E_no_interop.GetMaximum()));
  flux_given_E.Draw();
  flux_given_E_no_interop.SetLineColor(kBlue);
  flux_given_E_no_interop.SetLineStyle(kDashed);
  flux_given_E_no_interop.SetLineWidth(4);
  flux_given_E_no_interop.Draw("same");

  c1.cd(3);
  // flux_given_E_costh.SetMinimum(0);
  // flux_given_E_costh.SetMaximum(1.1 * std::max(
  //     flux_given_E_costh.GetMaximum(),
  //     flux_given_E_costh_no_interop.GetMaximum()));
  flux_given_E_costh.Draw();
  flux_given_E_costh_no_interop.SetLineColor(kBlue);
  flux_given_E_costh_no_interop.SetLineStyle(kDashed);
  flux_given_E_costh_no_interop.SetLineWidth(4);
  flux_given_E_costh_no_interop.Draw("same");

  c1.Draw();

  app.Run();
}