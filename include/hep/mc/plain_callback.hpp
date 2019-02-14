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
#include <functional>
#include <iostream>
#include <vector>

namespace hep
{

/// \addtogroup callbacks
/// @{

/// The default callback function. This function does nothing and always returns `true`.
///
/// \see \ref plain_callback
template <typename T>
inline bool plain_default_callback(plain_chkpt<T> const&)
{
    return true;
}

/// Callback function that prints a detailed summary about every iteration performed so far. This
/// function always returns `true`.
///
/// \see \ref plain_callback
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

/// The type of callback function that can be set by the user with \ref plain_callback.
template <typename T>
using plain_callback_type = std::function<bool(plain_chkpt<T> const&)>;

/// Sets the plain `callback` function and returns it. This function is called after each iteration
/// performed by \ref plain. The default callback is \ref plain_default_callback. The function can
/// e.g. be set to \ref plain_verbose_callback which prints after each iteration. If the callback
/// function returns `false` the integration is stopped.
///
/// If this function is called without any argument, the previous function is retained.
template <typename T>
inline plain_callback_type<T> plain_callback(plain_callback_type<T> callback = nullptr)
{
    static plain_callback_type<T> object = plain_default_callback<T>;

    if (callback != nullptr)
    {
        object = callback;
    }

    return object;
}

/// @}

}

#endif
