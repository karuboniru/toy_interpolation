#include "hkkm_reader.hxx"
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

  HKKM_READER_2D reader(input_file);
  auto &numu = reader[14];
  constexpr size_t n_logE_points = 101;
  constexpr size_t n_costh_bins = 20;
  constexpr size_t n_costh_points = n_costh_bins + 1;

  // here is where the interpolation method initializes
  axis_object logE_points{.min = -1, .max = 4, .n_points = n_logE_points};
  axis_object costh_points{.min = -1, .max = 1, .n_points = n_costh_points};
  interpolate<2, 4> interp{{logE_points, costh_points}};
  for (size_t i = 0; i < n_logE_points; ++i) {
    double cdf_along_costh = 0;
    interp[{i, 0}] = 0;
    for (size_t j = 0; j < n_costh_bins; ++j) {
      cdf_along_costh += numu.GetBinContent(i + 1, j + 1);
      interp[{i, j + 1}] = cdf_along_costh * 0.1;
    }
  }
  // then plot

  auto get_flux = [&](double E, double costh) {
    auto logE = std::log10(E);
    return interp.do_interpolation({logE, costh}, {false, true});
  };
  auto get_flux_no_interop = [&](double E, double costh) {
    auto logE = std::log10(E);
    // return numu.Interpolate(logE, costh);
    auto binx = numu.GetXaxis()->FindBin(logE);
    auto biny = numu.GetYaxis()->FindBin(costh);
    return numu.GetBinContent(binx, biny);
  };
  TApplication app("app", nullptr, nullptr);
  TCanvas c1("c1", "c1", 1600, 600);
  c1.Divide(2, 1);
  c1.cd(1);
  const double fixed_costh = 0.9;
  TF1 flux_given_costh{"flux",
                       [&](double *x, double *) {
                         return std::log(
                             get_flux(std::pow(10, x[0]), fixed_costh));
                       },
                       -1, 4, 0};
  flux_given_costh.SetTitle(
      std::format("Flux along cos#theta = {:.1f}; log_{{10}}(E/GeV); ln (#Phi/(m^{{2}} sec sr GeV)) ",
                  fixed_costh)
          .c_str());
  flux_given_costh.SetNpx(1000);
  TF1 flux_given_E{
      "flux", [&](double *x, double *) { return get_flux(1, x[0]); }, -1, 1, 0};
  flux_given_E.SetNpx(1000);
  flux_given_E.SetTitle("Flux along E = 1 GeV; cos#theta; #Phi (m^{{2}} sec sr GeV)");
  TF1 flux_given_costh_no_interop{
      "flux",
      [&](double *x, double *) {
        return std::log(get_flux_no_interop(std::pow(10, x[0]), fixed_costh));
      },
      -1, 4, 0};
  flux_given_costh_no_interop.SetNpx(1000);
  TF1 flux_given_E_no_interop{
      "flux", [&](double *x, double *) { return get_flux_no_interop(1, x[0]); },
      -1, 1, 0};
  flux_given_E_no_interop.SetNpx(1000);

  flux_given_costh.Draw();
  flux_given_costh_no_interop.SetLineColor(kBlue);
  flux_given_costh_no_interop.SetLineStyle(kDashed);
  flux_given_costh_no_interop.SetLineWidth(2);
  flux_given_costh_no_interop.Draw("same");
  c1.cd(2);
  // flux_given_E.SetMinimum(0);
  flux_given_E.Draw();
  flux_given_E_no_interop.SetLineColor(kBlue);
  flux_given_E_no_interop.SetLineStyle(kDashed);
  flux_given_E_no_interop.SetLineWidth(2);
  flux_given_E_no_interop.Draw("same");
  c1.Draw();

  app.Run();
}