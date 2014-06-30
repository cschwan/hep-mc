#ifndef HEP_MC_VEGAS_CALLBACK_HPP
#define HEP_MC_VEGAS_CALLBACK_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2014  Christopher Schwan
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

#include "hep/mc/vegas_iteration_result.hpp"
#include "hep/mc/mc_helper.hpp"

#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

namespace hep
{

/// \addtogroup vegas_callback_group
/// @{

/**
 * The default callback function. This function does nothing and always returns `true`.
 *
 * \see vegas_callback
 */
template <typename T>
inline bool vegas_default_callback(std::vector<vegas_iteration_result<T>> const&)
{
	return true;
}

/**
 * Callback function that prints a detailed summary about every iteration performed so far. This
 * function always returns `true`.
 *
 * \see vegas_callback
 */
template <typename T>
inline bool vegas_verbose_callback(std::vector<vegas_iteration_result<T>> const& results)
{
	std::cout << "iteration " << (results.size()-1) << " finished.\n";

	// print result for this iteration
	std::cout << "this iteration: N=" << results.back().calls << " E=" << results.back().value;
	std::cout << " +- " << results.back().error << " (";
	std::cout << (T(100.0) * results.back().error / std::abs(results.back().value)) << "%)\n";

	// compute cumulative results
	auto const result = cumulative_result(results.begin(), results.end());
	T const chi = chi_square_dof(results.begin(), results.end());

	// print the combined result
	std::cout << "all iterations: N=" << result.calls << " E=" << result.value << " +- ";
	std::cout << result.error << " (" << (T(100.0) * result.error / std::abs(result.value));
	std::cout << "%) chi^2/dof=" << chi << "\n\n";

	std::cout.flush();

	return true;
}

/// The type of callback function that can be set by the user with \ref vegas_callback.
template <typename T>
using vegas_callback_type = std::function<bool(std::vector<vegas_iteration_result<T>>)>;

/**
 * Sets the vegas `callback` function and returns it. This function is called after each iteration
 * performed by \ref vegas. The default callback is \ref vegas_default_callback. The function can
 * e.g. be set to \ref vegas_verbose_callback which prints after each iteration. If the callback
 * function returns `false` the integration is stopped.
 *
 * If this function is called without any argument, the previous function is retained.
 */
template <typename T>
inline vegas_callback_type<T> vegas_callback(vegas_callback_type<T> callback = nullptr)
{
	static vegas_callback_type<T> object = vegas_default_callback<T>;

	if (callback != nullptr)
	{
		object = callback;
	}

	return object;
}

/// @}

}

#endif
