#include "interpolation.hxx"

#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>

int main() {
  axis_object x_points{.min = 0.0, .max = 1.0, .n_points = 11};
  axis_object y_points{.min = 0.0, .max = 1.0, .n_points = 11};
  auto func = [](double x, double y) { return (std::sin(x) * std::cos(y)); };
  auto func_derivative_x = [](double x, double y) {
    return -std::cos(x) * std::sin(y);
  };

  try {
    interpolate<4, 4> interp({x_points, y_points});
    for (size_t i = 0; i < x_points.n_points; ++i) {
      for (size_t j = 0; j < y_points.n_points; ++j) {
        interp[{i, j}] = func(
            x_points.min +
                (i * (x_points.max - x_points.min) / (x_points.n_points - 1)),
            y_points.min +
                (j * (y_points.max - y_points.min) / (y_points.n_points - 1)));
      }
    }
    double x = 0.9322;
    double y = 0.52321;
    // std::cout << interp.do_interpolation({x, y}) << '\n';
    // std::cout << "expected: " << func(x, y) << '\n';
    double interpolated = interp.do_interpolation({x, y});
    double expected = func(x, y);
    std::println("Interpolated: {}, Expected: {}, RelErr: {}", interpolated,
                 expected, std::abs(interpolated - expected) / expected);
    interpolated = interp.do_interpolation({x, y}, {true, true});
    expected = func_derivative_x(x, y);
    std::println("Interpolated: {}, Expected: {}, RelErr: {}", interpolated,
                 expected, std::abs(interpolated - expected) / expected);

    TApplication app("app", nullptr, nullptr);
    TCanvas c1("c1", "c1", 1600, 600);
    c1.Divide(2, 1);
    c1.cd(1);
    TF1 orig_func{"orig_func",
                  [&](double *x, double *) { return func(x[0], y); }, 0, 1, 0};
    TF1 func(
        "func",
        [&](double *x, double *) { return interp.do_interpolation({x[0], y}); },
        0, 1, 0);
    orig_func.SetLineColor(kRed);
    orig_func.Draw();
    func.SetLineColor(kBlue);
    func.SetLineStyle(kDashed);
    func.SetLineWidth(4);
    func.Draw("same");
    c1.cd(2);
    TF1 orig_func_derivative_x{
        "orig_func_derivative_x",
        [&](double *x, double *) { return func_derivative_x(x[0], y); }, 0, 1,
        0};
    TF1 func_derivative_x{"func_derivative_x",
                          [&](double *x, double *) {
                            return interp.do_interpolation({x[0], y},
                                                           {true, false});
                          },
                          0, 1, 0};
    orig_func_derivative_x.SetLineColor(kRed);
    orig_func_derivative_x.Draw();
    func_derivative_x.SetLineColor(kBlue);
    func_derivative_x.SetLineStyle(kDashed);
    func_derivative_x.SetLineWidth(4);
    func_derivative_x.Draw("same");
    c1.cd(0);
    c1.Update();
    c1.Draw();
    app.Run();
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
  }
  return 0;
}