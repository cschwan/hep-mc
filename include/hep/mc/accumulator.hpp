#ifndef HEP_MC_ACCUMULATOR_HPP
#define HEP_MC_ACCUMULATOR_HPP

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

#include "hep/mc/accumulator_fwd.hpp"
#include "hep/mc/distribution_parameters.hpp"
#include "hep/mc/distribution_result.hpp"
#include "hep/mc/plain_result.hpp"
#include "hep/mc/projector.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace hep
{

/// \cond INTERNAL

template <typename T>
inline void accumulate(T& sum, T& sum_of_squares, T& compensation, T value)
{
    T const y = value - compensation;
    T const t = sum + y;
    compensation = (t - sum) - y;
    sum = t;

    sum_of_squares += value * value;
}

template <typename T>
class accumulator<T, true>
{
public:
    explicit accumulator(std::vector<hep::distribution_parameters<T>> const& parameters)
        : parameters_(parameters)
        , sums_()
        , compensations_()
    {
        std::size_t index = 2;
        indices_.reserve(parameters.size());

        for (auto const& params : parameters)
        {
            indices_.push_back(index);
            index += 2 * params.bins_x() * params.bins_y();
        }

        sums_.resize(2 + index);
        compensations_.resize(index / 2);
        non_zero_calls_.resize(index / 2);
        finite_calls_.resize(index / 2);
    }

    template <typename I, typename P>
    T invoke(I& integrand, P const& point)
    {
        using std::isfinite;

        hep::projector<T> projector(this, point);

        // call the integrand function with the supplied point. Distributions
        // are generated here
        T value = integrand.function()(point, projector);

        if (value != T())
        {
            value *= point.weight();

            if (isfinite(value))
            {
                accumulate(sums_[0], sums_[1], compensations_[0], value);
                ++finite_calls_[0];
            }
            else
            {
                value = T();
            }

            ++non_zero_calls_[0];
        }

        return value;
    }

    void add_to_1d_distribution(std::size_t index, T x, T value)
    {
        using std::isfinite;

        if (!isfinite(value))
        {
            return;
        }

        // TODO: index might be larger than the than allowed; throw?
        auto const parameters = parameters_.at(index);

        T const shifted_x = x - parameters.x_min();

        if (shifted_x < T())
        {
            // point is outside the binning range
            return;
        }

        std::size_t const bin_x = shifted_x / parameters.bin_size_x();

        if (bin_x >= parameters.bins_x())
        {
            // point is right of the range that we are binning
            return;
        }

        std::size_t const new_index = indices_.at(index) + 2 * bin_x;

        accumulate(
            sums_.at(new_index),
            sums_.at(new_index + 1),
            compensations_.at(new_index / 2),
            value
        );

        // FIXME: if this function is called more than once, the values are
        // incorrect
        ++non_zero_calls_.at(new_index / 2);
        ++finite_calls_.at(new_index / 2);
    }

    void add_to_2d_distribution(std::size_t index, T x, T y, T value)
    {
        using std::isfinite;

        if (!isfinite(value))
        {
            return;
        }

        // TODO: index might be larger than the than allowed; throw?
        auto const parameters = parameters_.at(index);

        T const shifted_x = x - parameters.x_min();

        if (shifted_x < T())
        {
            // point is outside the binning range
            return;
        }

        T const shifted_y = y - parameters.y_min();

        if (shifted_y < T())
        {
            // point is outside the binning range
            return;
        }

        std::size_t const bin_x = shifted_x / parameters.bin_size_x();

        if (bin_x >= parameters.bins_x())
        {
            // point is right of the range that we are binning
            return;
        }

        std::size_t const bin_y = shifted_y / parameters.bin_size_y();

        if (bin_y >= parameters.bins_y())
        {
            return;
        }

        std::size_t const new_index = indices_.at(index) + 2 * (bin_y *
            parameters.bins_x() + bin_x);

        accumulate(
            sums_.at(new_index),
            sums_.at(new_index + 1),
            compensations_.at(new_index / 2),
            value
        );

        // FIXME: if this function is called more than once, the values are
        // incorrect
        ++non_zero_calls_.at(new_index / 2);
        ++finite_calls_.at(new_index / 2);
    }

    hep::plain_result<T> result(std::size_t calls) const
    {
        std::vector<hep::distribution_result<T>> result;
        result.reserve(parameters_.size());

        std::size_t index = 2;

        // loop over all distributions
        for (auto const& params : parameters_)
        {
            std::vector<hep::mc_result<T>> bin_results;
            bin_results.reserve(params.bins_x());

            T const inv_bin_size = T(1.0) / params.bin_size_x() / params.bin_size_y();
            std::size_t const bins = params.bins_x() * params.bins_y();

            // loop over the bins of the current distribution
            for (std::size_t bin = 0; bin != bins; ++bin)
            {
                bin_results.emplace_back(
                    calls,
                    non_zero_calls_[index / 2],
                    finite_calls_[index / 2],
                    inv_bin_size                * sums_[index],
                    inv_bin_size * inv_bin_size * sums_[index + 1]
                );

                index += 2;
            }

            result.emplace_back(params, bin_results);
        }

        return hep::plain_result<T>(
            result,
            calls,
            non_zero_calls_[0],
            finite_calls_[0],
            sums_[0],
            sums_[1]
        );
    }

private:
    std::vector<hep::distribution_parameters<T>> parameters_;
    std::vector<std::size_t> indices_;
    std::vector<T> sums_;
    std::vector<T> compensations_;
    std::vector<std::size_t> non_zero_calls_;
    std::vector<std::size_t> finite_calls_;
};

template <typename T>
class accumulator<T, false>
{
public:
    explicit accumulator(std::vector<hep::distribution_parameters<T>> const& /*parameters*/)
        : sums_()
        , non_zero_calls_{}
        , finite_calls_{}
    {
    }

    template <typename I, typename P>
    T invoke(I& integrand, P const& point)
    {
        using std::isfinite;

        // call the integrand function with the supplied point. No distributions
        // are generated here
        T value = integrand.function()(point);

        if (value != T())
        {
            value *= point.weight();

            if (isfinite(value))
            {
                accumulate(sums_[0], sums_[1], sums_[2], value);
                ++finite_calls_;
            }
            else
            {
                value = T();
            }

            ++non_zero_calls_;
        }

        return value;
    }

    hep::plain_result<T> result(std::size_t calls) const
    {
        return hep::plain_result<T>(
            std::vector<hep::distribution_result<T>>{},
            calls,
            non_zero_calls_,
            finite_calls_,
            sums_[0],
            sums_[1]
        );
    }

private:
    std::array<T, 3> sums_;
    std::size_t non_zero_calls_;
    std::size_t finite_calls_;
};

template <typename I>
inline accumulator<typename I::numeric_type, I::has_distributions> make_accumulator(
    I const& integrand
) {
    using T = typename I::numeric_type;
    constexpr bool has_distributions = I::has_distributions;

    return accumulator<T, has_distributions>(integrand.parameters());
}

template <typename T>
inline void projector<T>::add(std::size_t index, T x, T value)
{
    accumulator_->add_to_1d_distribution(index, x, value * point_.weight());
}

template <typename T>
inline void projector<T>::add(std::size_t index, T x, T y, T value)
{
    accumulator_->add_to_2d_distribution(index, x, y, value * point_.weight());
}

/// \endcond

}

#endif
