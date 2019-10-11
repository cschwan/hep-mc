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

/// Class capturing the complete internal state of the \ref multi_channel_group.
template <typename T>
class multi_channel_chkpt : public chkpt<multi_channel_result<T>>
{
public:
    /// Constructor. Do not use directly, but instead use \ref make_multi_channel_chkpt.
    multi_channel_chkpt(T min_weight, T beta)
        : beta_{beta}
        , min_weight_{min_weight}
    {
    }

    /// Constructor. Do not use directly, but instead use \ref make_multi_channel_chkpt.
    multi_channel_chkpt(std::vector<T> const& channel_weights, T min_weight, T beta)
        : beta_{beta}
        , min_weight_{min_weight}
        , first_channel_weights_(multi_channel_refine_weights(channel_weights,
            std::vector<T>(channel_weights.size(), T(1.0)), min_weight_, beta_))
    {
    }

    /// Deserialization constructor. Do not use directly, but instead use \ref
    /// make_multi_channel_chkpt.
    explicit multi_channel_chkpt(std::istream& in)
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
    std::vector<T> channel_weights() const
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

    /// Returns the parameter \f$ \beta \f$ used to refine the weights with \ref
    /// multi_channel_refine_weights.
    T beta() const
    {
        return beta_;
    }

    /// The smallest weight allowed for each channel. This does not concern weights which are
    /// exactly zero, which disable channels.
    T min_weight() const
    {
        return min_weight_;
    }

    void serialize(std::ostream& out) const override
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
    T beta_;
    T min_weight_;
    std::vector<T> first_channel_weights_;
};

/// Multi channel checkpoint with random number generators.
template <typename RandomNumberEngine, typename T>
using multi_channel_chkpt_with_rng = chkpt_with_rng<RandomNumberEngine, multi_channel_chkpt<T>>;

/// Creates a checkpoint that can be used to start a multi channel integration.
template <typename T, typename RandomNumberEngine = std::mt19937>
multi_channel_chkpt_with_rng<RandomNumberEngine, T> make_multi_channel_chkpt(
    T min_weight = T(),
    T beta = T(0.25),
    RandomNumberEngine const& rng = RandomNumberEngine()
) {
    return multi_channel_chkpt_with_rng<RandomNumberEngine, T>{rng, min_weight, beta};
}

/// Creates a checkpoint that can be used to start a multi channel integration.
template <typename T, typename RandomNumberEngine = std::mt19937>
multi_channel_chkpt_with_rng<RandomNumberEngine, T> make_multi_channel_chkpt(
    std::vector<T> const& channel_weights,
    T min_weight = T(),
    T beta = T(0.25),
    RandomNumberEngine const& rng = RandomNumberEngine()
) {
    return multi_channel_chkpt_with_rng<RandomNumberEngine, T>{rng, channel_weights, min_weight,
         beta};
}

/// Helper function create a checkpoint reading from the stream `in`. Note the the numeric type `T`
/// as well as the type of the random number generator, `RandomNumberEngine` have to explicitly
/// stated.
template <typename T, typename RandomNumberEngine>
multi_channel_chkpt_with_rng<RandomNumberEngine, T> make_multi_channel_chkpt(std::istream& in)
{
    if (in.peek() == std::istream::traits_type::eof())
    {
        return make_multi_channel_chkpt<T, RandomNumberEngine>();
    }

    return multi_channel_chkpt_with_rng<RandomNumberEngine, T>{in};
}

/// Return type of \ref make_multi_channel_chkpt with default parameters.
template <typename T>
using default_multi_channel_chkpt = decltype (make_multi_channel_chkpt<T>());

/// @}

}

#endif
