#ifndef HEP_MC_MC_HELPER_HPP
#define HEP_MC_MC_HELPER_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2013  Christopher Schwan
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

#include <hep/mc/mc_result.hpp>

namespace hep
{

/**
 * \defgroup mc_helper Monte Carlo Helper Functions
 * @{
 */

/**
 * Computes a cumulative result using a range of results pointed to by `begin`
 * and `end`. If \f$ E_i, S_i, N_i \f$ are the estimate, error and number of
 * calls of each iteration in the range defined by [begin, end) then the
 * cumulative result is computed as:
 * \f{align}{
 * E &= S^2 \sum_i \frac{E_i}{S_i^2} \\
 * S &= \left( \sum_i \frac{1}{S_i^2} \right)^{-\frac{1}{2}} \\
 * N &= \sum_i N_i
 * \f}
 */
template <typename T, typename MCResultIterator>
mc_result<T> cumulative_result(MCResultIterator begin, MCResultIterator end)
{
	std::size_t calls = 0;
	T estimate = T();
	T variance = T();

	for (MCResultIterator i = begin; i != end; ++i)
	{
		T const tmp = T(1.0) / (i->error * i->error);
		calls += i->calls;
		variance += tmp;
		estimate += tmp * i->value;
	}

	variance  = T(1.0) / variance;
	estimate *= variance;

	return mc_result<T>(
		calls,
		T(calls) * estimate,
		T(calls) * (estimate * estimate + T(calls) * variance)
	);
}

/**
 * Returns an approximation for the \f$ \chi^2 \f$ per degree of freedom using
 * the results pointed to by the range [begin, end). The cumulative value
 * \f$ E \f$ is obtained by calling \ref cumulative_result():
 * \f[
 * \chi^2 / \mathrm{dof} \approx \frac{1}{n-1} \sum_{i=1}^n \frac{\left( E_i -
 * E \right)^2}{S_i^2}
 * \f]
 */
template <typename T, typename MCResultIterator>
T chi_square_dof(MCResultIterator begin, MCResultIterator end)
{
	mc_result<T> const result = cumulative_result<T>(begin, end);
	T sum = T();
	std::size_t n = 0;

	for (MCResultIterator i = begin; i != end; ++i)
	{
		T const tmp = i->value - result.value;
		sum += tmp * tmp / (i->error * i->error);
		++n;
	}

	return sum / T(n - 1);
}

/**
 * @}
 */

}

#endif
