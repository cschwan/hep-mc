#ifndef HEP_MC_VEGAS_HPP
#define HEP_MC_VEGAS_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2012-2019  Christopher Schwan
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

#include "hep/mc/accumulator.hpp"
#include "hep/mc/callback.hpp"
#include "hep/mc/integrand.hpp"
#include "hep/mc/vegas_chkpt.hpp"
#include "hep/mc/vegas_pdf.hpp"
#include "hep/mc/vegas_point.hpp"
#include "hep/mc/vegas_result.hpp"

#include <cstddef>
#include <limits>
#include <random>
#include <vector>

namespace hep
{

/// \addtogroup vegas_group
/// @{

/// Performs one VEGAS iteration. This integrates `function` over the unit-hypercube using `calls`
/// function evaluations with random numbers generated by `generator`. The generator is not seeded.
/// The `pdf` is used to implement importance sampling; stratified sampling is not used. The
/// dimension of the function is determined by `pdf.size()`.
///
/// The parameter `total_calls` determines the sample size \f$ N \f$ of the iteration \ref
/// vegas_iteration performs and therefore usually has usually the same value as `calls`. If VEGAS
/// is run in parallel, then this function will be called multiple times with a differently seeded
/// generator and with `calls` parameters each smaller than `total_calls` but their sum being equal
/// to `total_calls`.
template <typename I, typename R>
inline vegas_result<numeric_type_of<I>> vegas_iteration(
    I&& integrand,
    std::size_t calls,
    vegas_pdf<numeric_type_of<I>> const& pdf,
    R& generator
) {
    using T = numeric_type_of<I>;

    auto accumulator = make_accumulator(integrand);

    std::size_t const dimensions = pdf.dimensions();
    std::size_t const bins       = pdf.bins();

    std::vector<T> adjustment_data(dimensions * bins);
    std::vector<T> random_numbers(dimensions);
    std::vector<std::size_t> bin(dimensions);

    for (std::size_t i = 0; i != calls; ++i)
    {
        for (std::size_t j = 0; j != dimensions; ++j)
        {
            random_numbers[j] = std::generate_canonical<T,
                std::numeric_limits<T>::digits>(generator);
        }

        vegas_point<T> const point(random_numbers, bin, pdf);

        T const value = accumulator.invoke(integrand, point);
        T const square = value * value;

        // save square for each bin in order to refine the pdf later
        for (std::size_t j = 0; j != dimensions; ++j)
        {
            adjustment_data[j * bins + point.bin()[j]] += square;
        }
    }

    return vegas_result<T>(accumulator.result(calls), pdf, adjustment_data);
}

/// Integrates `function` by performing `iteration_calls.size()` iterations of the VEGAS algorithm,
/// with as many function calls for each iteration specified by the corresponding value in
/// `iteration_calls`. The pdf refinement is done using \ref vegas_refine_pdf which is called with
/// the \f$ \alpha \f$-parameter given by `alpha`.
///
/// This function can be used to start from an already adapted pdf, e.g. one by \ref
/// vegas_result.pdf obtained by a previous \ref vegas call.
template <typename I, typename Checkpoint = default_vegas_chkpt<numeric_type_of<I>>,
    typename Callback = callback<Checkpoint>>
inline Checkpoint vegas(
    I&& integrand,
    std::vector<std::size_t> const& iteration_calls,
    Checkpoint chkpt = make_vegas_chkpt<numeric_type_of<I>>(),
    Callback callback = hep::callback<Checkpoint>()
) {
    chkpt.dimensions(integrand.dimensions());

    auto generator = chkpt.generator();

    // perform iterations
    for (auto const calls : iteration_calls)
    {
        auto const& pdf = chkpt.pdf();
        auto const& result = vegas_iteration(integrand, calls, pdf, generator);

        chkpt.add(result, generator);

        if (!callback(chkpt))
        {
            break;
        }
    }

    return chkpt;
}

/// @}

}

#endif
