#ifndef HEP_MC_DISTRIBUTION_RESULT_HPP
#define HEP_MC_DISTRIBUTION_RESULT_HPP

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

#include "hep/mc/distribution_parameters.hpp"
#include "hep/mc/mc_result.hpp"

#include <istream>
#include <ostream>
#include <vector>

namespace hep
{

/// \addtogroup distributions
/// @{

/// Captures the result of the integration of a distribution.
template <typename T>
class distribution_result
{
public:
    /// Constructor.
    distribution_result(
        distribution_parameters<T> const& parameters,
        std::vector<mc_result<T>> const& results
    )
        : parameters_(parameters)
        , results_(results)
    {
    }

    /// Deserialization constructor.
    explicit distribution_result(std::istream& in)
        : parameters_{in}
    {
        std::size_t const size = parameters_.bins_x() * parameters_.bins_y();
        results_.reserve(size);

        for (std::size_t i = 0; i != size; ++i)
        {
            results_.emplace_back(in);
        }
    }

    /// Returns the parameters associated with this distribution.
    distribution_parameters<T> const& parameters() const
    {
        return parameters_;
    }

    /// Returns the result for each bin, corresponding to the bin positions returned by \ref
    /// mid_points_x and \ref mid_points_y.
    std::vector<mc_result<T>> const& results() const
    {
        return results_;
    }

    /// Serializes this object.
    void serialize(std::ostream& out) const
    {
        parameters_.serialize(out);

        for (std::size_t i = 0; i != results_.size(); ++i)
        {
            out << '\n';
            results_.at(i).serialize(out);
        }
    }

private:
    distribution_parameters<T> parameters_;
    std::vector<mc_result<T>> results_;
};

/// Returns the middle point of each bin of this distribution in the x-direction.
template <typename T>
inline std::vector<T> mid_points_x(distribution_result<T> const& result)
{
    auto const& parameters = result.parameters();

    std::vector<T> mid_points;
    mid_points.reserve(parameters.bins_x() * parameters.bins_y());

    for (std::size_t bin_y = 0; bin_y != parameters.bins_y(); ++bin_y)
    {
        T x = parameters.x_min() + T(0.5) * parameters.bin_size_x();

        for (std::size_t bin = 0; bin != parameters.bins_x(); ++bin)
        {
            mid_points.push_back(x);
            x += parameters.bin_size_x();
        }
    }

    return mid_points;
}

/// Returns the middle point of each bin of this distribution in the y-direction.
template <typename T>
inline std::vector<T> mid_points_y(distribution_result<T> const& result)
{
    auto const& parameters = result.parameters();

    std::vector<T> mid_points;
    mid_points.reserve(parameters.bins_x() * parameters.bins_y());

    T y = parameters.y_min() + T(0.5) * parameters.bin_size_y();

    for (std::size_t bin_y = 0; bin_y != parameters.bins_y(); ++bin_y)
    {
        for (std::size_t bin = 0; bin != parameters.bins_x(); ++bin)
        {
            mid_points.push_back(y);
        }

        y += parameters.bin_size_y();
    }

    return mid_points;
}

/// @}

}

#endif
