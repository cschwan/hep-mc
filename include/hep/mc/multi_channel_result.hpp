#ifndef HEP_MC_MULTI_CHANNEL_RESULT_HPP
#define HEP_MC_MULTI_CHANNEL_RESULT_HPP

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

#include "hep/mc/plain_result.hpp"

#include <cstddef>
#include <iomanip>
#include <ios>
#include <istream>
#include <limits>
#include <ostream>
#include <vector>

namespace hep
{

/// \addtogroup results
/// @{

/// Result of a multi-channel integration.
template <typename T>
class multi_channel_result : public plain_result<T>
{
public:
    /// Constructor. The parameter `adjustment_data` must contain two additional values, being the
    /// sum of the sum of squares of the integration.
    multi_channel_result(
        plain_result<T> const& result,
        std::vector<T> const& adjustment_data,
        std::vector<T> const& channel_weights
    )
        : plain_result<T>(result)
        , adjustment_data_(adjustment_data)
        , channel_weights_(channel_weights)
    {
    }

    /// Deserialization constructor.
    multi_channel_result(std::istream& in)
        : plain_result<T>(in)
    {
        std::size_t channels;
        in >> channels;
        adjustment_data_.resize(channels);
        channel_weights_.resize(channels);

        for (std::size_t i = 0; i != channels; ++i)
        {
            in >> adjustment_data_.at(i) >> channel_weights_.at(i);
        }
    }

    /// This is the data used by \ref multi_channel_refine_weights to refine the \ref
    /// channel_weights used in the same iteration. The refined weights are then used in a
    /// subsequent iteration.
    std::vector<T> const& adjustment_data() const
    {
        return adjustment_data_;
    }

    /// The weight for each channel.
    std::vector<T> const& channel_weights() const
    {
        return channel_weights_;
    }

    /// Serializes this object.
    virtual void serialize(std::ostream& out) const override
    {
        plain_result<T>::serialize(out);
        out << '\n' << channel_weights_.size();

        for (std::size_t i = 0; i != channel_weights_.size(); ++i)
        {
            out << '\n' << std::scientific
                << std::setprecision(std::numeric_limits<T>::max_digits10 - 1)
                << adjustment_data_.at(i) << ' ' << channel_weights_.at(i);
        }
    }

private:
    std::vector<T> adjustment_data_;
    std::vector<T> channel_weights_;
};

/// @}

}

#endif
