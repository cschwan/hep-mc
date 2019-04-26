#ifndef HEP_MC_PLAIN_CALLBACK_HPP
#define HEP_MC_PLAIN_CALLBACK_HPP

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

#include "hep/mc/mc_helper.hpp"
#include "hep/mc/plain_chkpt.hpp"

#include <cmath>
#include <iostream>
#include <vector>

namespace hep
{

/// \addtogroup callbacks
/// @{

/// The default callback function. This function does nothing and always returns `true`.
template <typename T>
inline bool plain_silent_callback(plain_chkpt<T> const&)
{
    return true;
}

/// Callback function that prints a detailed summary about every iteration performed so far. This
/// function always returns `true`.
template <typename T>
inline bool plain_verbose_callback(plain_chkpt<T> const& chkpt)
{
    using std::fabs;

    auto const& results = chkpt.results();

    std::cout << "iteration " << (results.size() - 1) << " finished.\n";

    T const relative_error_percent = (T(100.0) * results.back().error() /
        fabs(results.back().value()));

    T const efficiency = T(100.0) * T(results.back().non_zero_calls()) / T(results.back().calls());

    std::size_t const number_of_non_finite_calls = results.back().non_zero_calls() -
        results.back().finite_calls();

    // print result for this iteration
    std::cout << "this iteration: N=" << results.back().calls() << " E=" << results.back().value()
        << " +- " << results.back().error() << " (" << relative_error_percent << "%) eff="
        << efficiency << "% nnf=" << number_of_non_finite_calls << "\n";

    // compute cumulative results
    auto const result = accumulate<weighted_with_variance>(results.begin(), results.end());
    T const chi = chi_square_dof<weighted_with_variance>(results.begin(), results.end());

    T const relative_error_percent_all = (T(100.0) * result.error() / fabs(result.value()));

    // print the combined result
    std::cout << "all iterations: N=" << result.calls() << " E=" << result.value() << " +- "
        << result.error() << " (" << relative_error_percent_all << "%) chi^2/dof=" << chi << "\n\n";

    std::cout.flush();

    return true;
}

/// @}

}

#endif
