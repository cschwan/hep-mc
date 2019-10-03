#ifndef HEP_MC_MULTI_CHANNEL_POINT_HPP
#define HEP_MC_MULTI_CHANNEL_POINT_HPP

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

#include "hep/mc/mc_point.hpp"
#include "hep/mc/multi_channel_map.hpp"

#include <cstddef>
#include <vector>

namespace hep
{

/// \addtogroup integrands
/// @{

/// Point in the unit-hypercube for Multi Channel Monte Carlo integration.
template <typename T>
class multi_channel_point : public mc_point<T>
{
public:
    /// Constructor.
    multi_channel_point(
        std::vector<T> const& point,
        T weight,
        std::vector<T>& coordinates,
        std::size_t channel
    )
        : mc_point<T>(point, weight)
        , channel_(channel)
        , coordinates_(coordinates)
    {
    }

    /// There is no copy constructor.
    multi_channel_point(multi_channel_point<T> const&) = delete;

    /// There is no move constructor.
    multi_channel_point(multi_channel_point<T>&&) = delete;

    /// There is no copy assignment operator.
    multi_channel_point& operator=(multi_channel_point<T> const&) = delete;

    /// There is no move assignment operator.
    multi_channel_point& operator=(multi_channel_point<T>&&) = delete;

    /// Destructor.
    ~multi_channel_point() override = default;

    /// The selected channel for this point.
    std::size_t channel() const
    {
        return channel_;
    }

    /// The point in the hypercube transformed by the current \ref channel.
    std::vector<T> const& coordinates() const
    {
        return coordinates_;
    }

private:
    std::size_t channel_;
    std::vector<T>& coordinates_;
};

/// Point in the unit-hypercube for multi-channel Monte Carlo integration. This type also captures
/// the map that is used to generate `coordinates`.
template <typename T, typename M>
class multi_channel_point2 : public multi_channel_point<T>
{
public:
    /// Constructor.
    multi_channel_point2(
        std::vector<T> const& point,
        std::vector<T>& coordinates,
        std::size_t channel,
        std::vector<T>& densities,
        std::vector<T> const& channel_weights,
        std::vector<std::size_t> const& enabled_channels,
        M& map
    )
        : multi_channel_point<T>(point, T(), coordinates, channel)
        , densities_(densities)
        , channel_weights_(channel_weights)
        , enabled_channels_(enabled_channels)
        , map_(map)
    {
    }

    /// The map function that constructed this point. See \ref multi_channel_iteration for reference
    /// on the signature of this function.
    M const& map() const
    {
        return map_;
    }

    /// Returns the weight for this Monte Carlo point.
    T weight() const override
    {
        if (this->weight_ == T())
        {
            // lazy evaluation of the jacobian of `map` and `densities`
            this->weight_ = map_(
                this->channel(),
                this->point(),
                // should be OK, since `coordinates` is actually non-const as forced by the c'tor
                const_cast <std::vector<T>&> (this->coordinates()),
                enabled_channels_,
                densities_,
                multi_channel_map::calculate_densities
            );

            T total_density = T();

            for (std::size_t j = 0; j != channel_weights_.size(); ++j)
            {
                total_density += channel_weights_[j] * densities_[j];
            }

            this->weight_ /= total_density;
        }

        return this->weight_;
    }

private:
    std::vector<T>& densities_;
    std::vector<T> const& channel_weights_;
    std::vector<std::size_t> const& enabled_channels_;
    M& map_;
};

/// @}

}

#endif
