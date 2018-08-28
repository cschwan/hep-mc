#ifndef HEP_MC_DISTRIBUTION_PARAMETERS_HPP
#define HEP_MC_DISTRIBUTION_PARAMETERS_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2015-2018  Christopher Schwan
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstddef>
#include <string>

namespace hep
{

/// \addtogroup distributions
/// @{

/// Defines the parameters of a one- or two-dimensional distribution.
template <typename T>
class distribution_parameters
{
public:
    /// Constructor. Constructs a one-dimensional distribution.
    distribution_parameters(std::size_t bins, T x_min, T x_max, std::string const& name)
        : bins_x_{bins}
        , bins_y_{1}
        , x_min_{x_min}
        , y_min_{T()}
        , bin_size_x_{(x_max - x_min) / T(bins)}
        , bin_size_y_{T(1.0)}
        , name_{name}
    {
    }

    /// Constructor. Constructs a two dimensional distribution.
    distribution_parameters(
        std::size_t bins_x,
        std::size_t bins_y,
        T x_min,
        T x_max,
        T y_min,
        T y_max,
        std::string const& name
    )
        : bins_x_{bins_x}
        , bins_y_{bins_y}
        , x_min_{x_min}
        , y_min_{y_min}
        , bin_size_x_{(x_max - x_min) / T(bins_x)}
        , bin_size_y_{(y_max - y_min) / T(bins_y)}
        , name_{name}
    {
    }

    /// Returns the number of bins in x-direction.
    std::size_t bins_x() const
    {
        return bins_x_;
    }

    /// Returns the number of bins in y-direction.
    std::size_t bins_y() const
    {
        return bins_y_;
    }

    /// Name of the distribution.
    std::string const& name() const
    {
        return name_;
    }

    /// Smallest x-value of a point that would still be accumulated into the distribution.
    T x_min() const
    {
        return x_min_;
    }

    /// Smallest y-value of a point that would still be accumulated into the distribution.
    T y_min() const
    {
        return x_min_;
    }

    /// Size of the bins in x-direction.
    T bin_size_x() const
    {
        return bin_size_x_;
    }

    /// Size of the bins in y-direction.
    T bin_size_y() const
    {
        return bin_size_y_;
    }

private:
    std::size_t bins_x_;
    std::size_t bins_y_;
    T x_min_;
    T y_min_;
    T bin_size_x_;
    T bin_size_y_;
    std::string name_;
};

/// Shortcut for calling the constructor that automatically determines the numeric type.
template <typename T>
distribution_parameters<T> make_dist_params(
    std::size_t bins,
    T x_min,
    T x_max,
    std::string const& name = ""
) {
    return distribution_parameters<T>(bins, x_min, x_max, name);
}

/// @}

}

#endif
