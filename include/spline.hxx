#pragma once

#include <TF1.h>
#include <TFile.h>
#include <map>
#include <memory>
#include <string>
#include <tuple>

#include "interpolation.hxx"

using xsec_interpolater = interpolate<3>;

class spline_reader {
public:
  spline_reader(const char *filename);
  spline_reader(spline_reader &&) = default;
  spline_reader(const spline_reader &) = default;
  spline_reader &operator=(const spline_reader &) = default;
  spline_reader &operator=(spline_reader &&) = default;
  ~spline_reader();

  double
  get_cross_section(int neutrino_id, int target_id, double energy,
                    const std::string &channel_name = "total_comparable");

private:
  std::shared_ptr<TFile> file;
  std::map<
      std::tuple<int, int, std::string> /* neutrino_id, target_id, channel */,
      xsec_interpolater>
      cross_section_objects;
};