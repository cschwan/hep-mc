#ifndef HEP_MC_MULTI_CHANNEL_CHKPT_HPP
#define HEP_MC_MULTI_CHANNEL_CHKPT_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2018-2019  Christopher Schwan
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

#include "hep/mc/chkpt.hpp"
#include "hep/mc/multi_channel_refine_weights.hpp"
#include "hep/mc/multi_channel_result.hpp"

#include <cassert>
#include <cstddef>
#include <iomanip>
#include <ios>
#include <limits>
#include <random>
#include <vector>

namespace hep
{

/// \addtogroup checkpoints
/// @{

template <typename T>
class multi_channel_chkpt : public chkpt<multi_channel_result<T>>
{
public:
    multi_channel_chkpt(T min_weight, T beta)
        : min_weight_{min_weight}
        , beta_{beta}
    {
    }

    multi_channel_chkpt(std::vector<T> const& channel_weights, T min_weight, T beta)
        : first_channel_weights_(channel_weights)
        , min_weight_{min_weight}
        , beta_{beta}
    {
    }

    /// Deserialization constructor. This creates a checkpoint by reading from the stream `in`.
    multi_channel_chkpt(std::istream& in)
        : chkpt<multi_channel_result<T>>{in}
    {
        in >> beta_ >> min_weight_;

        if (this->results().empty())
        {
            std::size_t channels = 0;
            in >> channels;
            first_channel_weights_.reserve(channels);

            for (std::size_t i = 0; i != channels; ++i)
            {
                T weight;
                in >> weight;
                first_channel_weights_.push_back(weight);
            }
        }
    }

    /// Returns the channel weights for the next iteration.
    std::vector<T> const channel_weights() const
    {
        auto const& results = this->results();

        if (results.empty())
        {
            return first_channel_weights_;
        }

        return multi_channel_refine_weights(results.back().channel_weights(),
            results.back().adjustment_data(), min_weight_, beta_);
    }

    /// Sets the number of channels.
    void channels(std::size_t channels)
    {
        if (first_channel_weights_.empty())
        {
            first_channel_weights_.assign(channels, T(1.0) / channels);
        }

        assert( this->results().empty() ||
            (this->results().back().channel_weights().size() == channels) );
    }

    T beta() const
    {
        return beta_;
    }

    T min_weight() const
    {
        return min_weight_;
    }

    virtual void serialize(std::ostream& out) const override
    {
        chkpt<multi_channel_result<T>>::serialize(out);

        out << '\n' << std::scientific
            << std::setprecision(std::numeric_limits<T>::max_digits10 - 1) << beta_ << ' '
            << min_weight_;

        if (this->results().empty())
        {
            out << '\n' << first_channel_weights_.size();

            for (auto const weight : first_channel_weights_)
            {
                out << ' ' << std::scientific
                    << std::setprecision(std::numeric_limits<T>::max_digits10 - 1) << weight;
            }
        }
    }

private:
    std::vector<T> first_channel_weights_;
    T min_weight_;
    T beta_;
};

///
template <typename RandomNumberEngine, typename T>
using multi_channel_chkpt_with_rng = chkpt_with_rng<RandomNumberEngine, multi_channel_chkpt<T>>;

///
template <typename T, typename RandomNumberEngine = std::mt19937>
multi_channel_chkpt_with_rng<RandomNumberEngine, T> make_multi_channel_chkpt(
    T beta = T(0.25),
    T min_fraction = T(),
    RandomNumberEngine const& rng = std::mt19937()
) {
    return multi_channel_chkpt_with_rng<RandomNumberEngine, T>{rng, beta, min_fraction};
}

///
template <typename T, typename RandomNumberEngine = std::mt19937>
multi_channel_chkpt_with_rng<RandomNumberEngine, T> make_multi_channel_chkpt(
    std::vector<T> const& channel_weights,
    T beta = T(0.25),
    T min_fraction = T(),
    RandomNumberEngine const& rng = std::mt19937()
) {
    return multi_channel_chkpt_with_rng<RandomNumberEngine, T>{rng, channel_weights, beta,
        min_fraction};
}

///
template <typename T>
using default_multi_channel_chkpt = decltype (make_multi_channel_chkpt<T>());

/// @}

}

#endif
