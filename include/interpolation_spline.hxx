#include <array>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <print>
#include <vector>

#include <TSpline.h>

#include <Eigen/Dense>

struct axis_object {
  double min{};
  double max{};
  size_t n_points{};
};

// template parameter specifies the order of the polynomial
// to use for interpolation for each axis
template <size_t dimension> class interpolate {
public:
  // static constexpr size_t dimension = sizeof...(poly_orders);
  // static constexpr std::array<size_t, dimension> orders_array = {
  //     (poly_orders + 1)...};
  constexpr static size_t npoints = 5;
  // just for copy-pasting
  constexpr static size_t order = npoints;

  interpolate(const std::array<axis_object, dimension> &axes_)
      : axes(axes_), values([](const auto &axes) {
          size_t n_values = 1;
          for (const auto &axis : axes) {
            n_values *= axis.n_points;
          }
          return std::vector<double>(n_values);
        }(axes_)) {}

  interpolate(const interpolate &) = default;
  interpolate(interpolate &&) = default;
  interpolate &operator=(const interpolate &) = default;
  interpolate &operator=(interpolate &&) = default;
  ~interpolate() = default;

  auto &operator[](const std::array<size_t, dimension> &index_array) {
    size_t flat_index = 0;
    for (size_t i = 0; i < dimension; ++i) {
      flat_index *= axes[i].n_points;
      flat_index += index_array[i];
    }
    return values[flat_index];
  }

  const auto &
  operator[](const std::array<size_t, dimension> &index_array) const {
    size_t flat_index = 0;
    for (size_t i = 0; i < dimension; ++i) {
      flat_index *= axes[i].n_points;
      flat_index += index_array[i];
    }
    return values[flat_index];
  }

  double
  do_interpolation(const std::array<double, dimension> &coordinates,
                   const std::array<bool, dimension> &do_derivative = {},
                   const std::array<bool, dimension> &do_periodic = {}) const {
    auto grid = build_grid(coordinates, do_derivative, do_periodic);
    auto value = do_interpolation_from_grid(grid);
    for (size_t i = 0; i < dimension; ++i) {
      if (do_derivative[i]) {
        value /= (axes[i].max - axes[i].min) / (axes[i].n_points - 1) *
                 (order - 1) / 2;
      }
    }
    return value;
  }

private:
  std::array<axis_object, dimension> axes;
  std::vector<double> values;
  template <typename T, size_t order>
  void static add_grid(
      std::tuple<std::array<T, order>, bool, double> &grid1,
      const std::tuple<std::array<T, order>, bool, double> &grid2) {
    if constexpr (std::is_same_v<T, double>) {
      for (size_t i = 0; i < order; ++i) {
        std::get<0>(grid1)[i] += std::get<0>(grid2)[i];
      }
    } else {
      for (size_t i = 0; i < order; ++i) {
        add_grid(std::get<0>(grid1)[i], std::get<0>(grid2)[i]);
      }
    }
  }

  template <size_t dimension_index = 0>
  auto build_grid(const std::array<double, dimension> &coordinates,
                  const std::array<bool, dimension> &do_derivative,
                  const std::array<bool, dimension> &do_periodic,
                  const std::array<size_t, dimension_index> &index =
                      std::array<size_t, 0>{}) const {
    // constexpr size_t order = orders_array[dimension_index];
    double normalized_coordinates =
        (coordinates[dimension_index] - axes[dimension_index].min) /
        (axes[dimension_index].max - axes[dimension_index].min) *
        (axes[dimension_index].n_points - 1);
    auto this_do_derivative = do_derivative[dimension_index];
    auto this_do_periodic = do_periodic[dimension_index];

    auto closest_indices =
        this_do_periodic
            ? get_closest_indices_periodic(normalized_coordinates,
                                           axes[dimension_index].n_points)
            : get_closest_indices(normalized_coordinates,
                                  axes[dimension_index].n_points);

    // notice this value is defined in [-1, 1]
    double normalized_value{};
    if (!this_do_periodic)
      normalized_value =
          -1 + (2 * (normalized_coordinates - closest_indices[0]) /
                (closest_indices[order - 1] - closest_indices[0]));
    else {
      auto center_in_big_coord = order % 2 == 0 ? 0.5 : 1;
      auto shift_in_big_coord = normalized_coordinates - center_in_big_coord;
      shift_in_big_coord -= std::round(shift_in_big_coord);
      normalized_value = shift_in_big_coord * 2 / (order - 1);
    }

    std::array<size_t, dimension_index + 1> new_coord{};
    for (size_t i = 0; i < dimension_index; ++i) {
      new_coord[i] = index[i];
    }
    if constexpr (dimension_index == dimension - 1) {
      std::array<double, order> sub_grids;
      double CDF_shift = 0;
      for (size_t i = 0; i < order; ++i) {
        new_coord[dimension_index] = closest_indices[i];
        sub_grids[i] = (*this)[new_coord];
        if (this_do_derivative && this_do_periodic) {
          // in case from order - 2 to 0
          if (i > 0 && closest_indices[i] == 0 &&
              closest_indices[i - 1] == axes[dimension_index].n_points - 2) {
            new_coord[dimension_index] = axes[dimension_index].n_points - 1;
            CDF_shift = (*this)[new_coord];
          }
          sub_grids[i] += CDF_shift;
        }
      }
      return std::make_tuple(sub_grids, this_do_derivative, normalized_value);
    } else {
      using sub_grid_type =
          decltype(build_grid(coordinates, do_derivative, do_periodic,
                              std::array<size_t, dimension_index + 1>{}));
      std::array<sub_grid_type, order> ret{};
      sub_grid_type CDF_shift;
      for (size_t i = 0; i < order; ++i) {
        new_coord[dimension_index] = closest_indices[i];
        ret[i] = build_grid(coordinates, do_derivative, do_periodic, new_coord);
        if (this_do_derivative && this_do_periodic) {
          // in case from order - 2 to 0
          if (i > 0 && closest_indices[i] == 0 &&
              closest_indices[i - 1] == axes[dimension_index].n_points - 2) {
            new_coord[dimension_index] = axes[dimension_index].n_points - 1;
            CDF_shift =
                build_grid(coordinates, do_derivative, do_periodic, new_coord);
          }
          add_grid(ret[i], CDF_shift);
        }
      }
      return std::make_tuple(ret, this_do_derivative, normalized_value);
    }
  }

  template <class T>
  static double do_interpolation_from_grid(
      const std::tuple<std::array<T, npoints>, bool, double> &grid_values) {
    const auto &[sub_grids, do_derivative, value] = grid_values;
    std::array<double, npoints> t_value, yvalue;
    for (size_t i = 0; i < npoints; ++i) {
      if constexpr (std::is_same_v<T, double>) {
        yvalue[i] = sub_grids[i];
      } else {
        yvalue[i] = do_interpolation_from_grid(sub_grids[i]);
      }
      t_value[i] = -1 + 2. * i / (npoints - 1);
    }
    TSpline3 spline{"spline", t_value.data(), yvalue.data(), npoints};
    if (do_derivative) {
      return spline.Derivative(value);
    }
    return spline.Eval(value);
  }

  static std::array<size_t, npoints> get_closest_indices(double this_index_d,
                                                         size_t n_points) {
    std::array<size_t, npoints> ret{};
    std::iota(ret.begin(), ret.end(), 0);
    auto mid = (npoints - 1) / 2.0;
    auto shift = std::round(this_index_d - mid);
    if (shift < 0) {
      shift = 0;
    } else if (shift > n_points - npoints) {
      shift = n_points - npoints;
    }

    for (auto &i : ret) {
      i += shift;
    }
    return ret;
  }
  static std::array<size_t, npoints>
  get_closest_indices_periodic(double this_index_d, size_t n_points) {
    std::array<size_t, npoints> ret{};
    std::iota(ret.begin(), ret.end(), 0);
    auto mid = (npoints - 1) / 2.0;
    int shift = std::round(this_index_d - mid);
    auto n_points_used = n_points - 1;
    for (auto &i : ret) {
      if (shift < 0) {
        i = (i + shift + n_points_used) % n_points_used;
      } else {
        i = (i + shift) % n_points_used;
      }
    }
    return ret;
  }
};
