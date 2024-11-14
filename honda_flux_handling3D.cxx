#include "hkkm_reader_3D.hxx"
#include "interpolation.hxx"
#include <TAttLine.h>
#include <algorithm>
#include <boost/program_options.hpp>
#include <cmath>
#include <print>

#include <TCanvas.h>
#include <TF1.h>
#include <TF2.h>
#include <TStyle.h>

int main(int argc, char **argv) {
  namespace po = boost::program_options;
  po::options_description desc("Allowed options");
  desc.add_options()("help", "produce help message")                       //
      ("input-file,i", po::value<std::string>()->required(), "input file") //
      ("fix-E,e", po::value<double>()->default_value(1.6), "fixed energy") //
      ("fix-costh,c", po::value<double>()->default_value(0.85),
       "fixed costh") //
      ("fix-phi,p", po::value<double>()->default_value(1.0),
       "fixed phi") //
      ("output,o", po::value<std::string>()->default_value(""),
       "Prefix For output plots");
  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    std::cerr << desc << '\n';
    return 1;
  }
  if (vm.count("help")) {
    std::cout << desc << '\n';
    return 1;
  }
  auto input_file = vm["input-file"].as<std::string>();
  const double fixed_costh = vm["fix-costh"].as<double>();
  const double fixed_energy = vm["fix-E"].as<double>();
  const double fixed_phi = vm["fix-phi"].as<double>();

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

  interpolate<3, 4, 4> interp{{logE_points, costh_points, phi_points}};
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
    auto binx = numu.GetXaxis()->FindBin(logE);
    auto biny = numu.GetYaxis()->FindBin(costh);
    auto binz = numu.GetZaxis()->FindBin(phi);
    return numu.GetBinContent(binx, biny, binz);
  };

  // TApplication app("app", nullptr, nullptr);
  constexpr int res_factor = 1;
  TCanvas c1("c1", "c1", 800 * 3 * res_factor, 600 * res_factor);
  c1.Divide(3, 1);
  c1.cd(1);

  TF1 flux_given_costh{
      "flux",
      [&](double *x, double *) {
        return std::log(get_flux(std::pow(10, x[0]), fixed_costh, fixed_phi));
      },
      -1, 4, 0};
  flux_given_costh.SetTitle(
      std::format("Flux along cos #theta = {:.1f},  #phi = {:.1f}; "
                  "log_{{10}}(E/GeV); ln (#Phi/(m^{{2}} sec sr GeV) ",
                  fixed_costh, fixed_phi)
          .c_str());
  // flux_given_costh.SetNpx(1000);
  const double costh_thres = 1.0;
  TF1 flux_given_E{"flux",
                   [&](double *x, double *) {
                     return get_flux(fixed_energy, x[0], fixed_phi);
                   },
                   -costh_thres, costh_thres, 0};
  // flux_given_E.SetNpx(1000);
  flux_given_E.SetTitle(std::format("Flux along E = {} GeV,  #phi = {:.1f}; "
                                    "cos #theta; #Phi (m^{{2}} sec sr GeV)",
                                    fixed_energy, fixed_phi)
                            .c_str());

  TF1 flux_given_costh_no_interop{
      "flux",
      [&](double *x, double *) {
        return std::log(
            get_flux_no_interop(std::pow(10, x[0]), fixed_costh, fixed_phi));
      },
      -1, 4, 0};
  // flux_given_costh_no_interop.SetNpx(1000);
  TF1 flux_given_E_no_interop{"flux",
                              [&](double *x, double *) {
                                return get_flux_no_interop(fixed_energy, x[0],
                                                           fixed_phi);
                              },
                              -costh_thres, costh_thres, 0};
  // flux_given_E_no_interop.SetNpx(1000);

  TF1 flux_given_E_costh{"flux",
                         [&](double *x, double *) {
                           //  return std::log(
                           return get_flux(fixed_energy, fixed_costh, x[0]);
                         },
                         0, 2 * M_PI, 0};
  // flux_given_E_costh.SetNpx(1000);
  flux_given_E_costh.SetTitle(
      std::format("Flux along E = {} GeV, cos  #theta = "
                  "{:.1f}; #phi; #Phi (m^{{2}} sec sr GeV)",
                  fixed_energy, fixed_costh)
          .c_str());
  TF1 flux_given_E_costh_no_interop{"flux",
                                    [&](double *x, double *) {
                                      return get_flux_no_interop(
                                          fixed_energy, fixed_costh, x[0]);
                                    },
                                    0, 2 * M_PI, 0};
  // flux_given_E_costh_no_interop.SetNpx(1000);

  auto do_plot = [&](int id, TF1 &f, TF1 &f_no_interop) {
    f.SetNpx(500);
    f_no_interop.SetNpx(1000);
    auto xmin = f.GetXmin();
    auto xmax = f.GetXmax();
    xmin = xmin + 0.01 * (xmax - xmin);
    xmax = xmax - 0.01 * (xmax - xmin);

    auto min =
        std::min(f.GetMinimum(xmin, xmax), f_no_interop.GetMinimum(xmin, xmax));
    auto max =
        std::max(f.GetMaximum(xmin, xmax), f_no_interop.GetMaximum(xmin, xmax));
    min = min - 0.1 * std::abs(max - min);
    max = max + 0.1 * std::abs(max - min);

    f.SetMinimum(min);
    f.SetMaximum(max);
    c1.cd(id);
    f.Draw();
    f_no_interop.SetLineColor(kBlue);
    f_no_interop.SetLineStyle(kDotted);
    f_no_interop.SetLineWidth(3);
    f_no_interop.Draw("same");
  };

  do_plot(1, flux_given_costh, flux_given_costh_no_interop);
  do_plot(2, flux_given_E, flux_given_E_no_interop);
  do_plot(3, flux_given_E_costh, flux_given_E_costh_no_interop);

  c1.Draw();
  auto output_str = vm["output"].as<std::string>();
  if (output_str.size()) {
    c1.SaveAs(std::format("{}_flux.svg", output_str).c_str());
    c1.SaveAs(std::format("{}_flux.pdf", output_str).c_str());
    c1.SaveAs(std::format("{}_flux.eps", output_str).c_str());
  }
  // app.Run();
}