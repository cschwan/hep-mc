#ifndef HEP_MC_MC_HELPER_HPP
#define HEP_MC_MC_HELPER_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2013-2015  Christopher Schwan
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

#include "hep/mc/mc_result.hpp"

#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>

namespace
{

template <typename Iterator>
using get_T = typename std::iterator_traits<Iterator>::value_type::numeric_type;

}

namespace hep
{

/// \addtogroup results
/// @{

/**
 * Computes a cumulative result using a range of results pointed to by `begin`
 * and `end`. If \f$ E_i, S_i, N_i \f$ are the estimate, error and number of
 * calls of each iteration in the range defined by [`begin`, `end`) and \f$ M
 * \f$ the number of results, then the cumulative result is computed as:
 * \f{align}{
 *     E &= S^2 \sum_{i=1}^M \frac{E_i}{S_i^2} \\
 *     S &= \left( \sum_{i=1}^M \frac{1}{S_i^2} \right)^{-\frac{1}{2}} \\
 *     N &= \sum_{i=1}^M N_i
 * \f}
 */
template <typename Iterator>
inline mc_result<get_T<Iterator>> cumulative_result0(
	Iterator begin, Iterator end
) {
	typedef get_T<Iterator> T;

	std::size_t calls = 0;
	T estimate = T();
	T variance = T();

	for (Iterator i = begin; i != end; ++i)
	{
		T const tmp = T(1.0) / i->variance();
		calls += i->calls();
		variance += tmp;
		estimate += tmp * i->value();
	}

	variance = T(1.0) / variance;
	estimate *= variance;

	return create_result(calls, estimate, std::sqrt(variance));
}

/**
 * Computes a cumulative result using a range of results pointed to by `begin`
 * and `end`. If \f$ E_i, S_i, N_i \f$ are the estimate, error and number of
 * calls of each iteration in the range defined by [`begin`, `end`) and \f$ M
 * \f$ the number of results, then the cumulative result is computed as:
 * \f{align}{
 *     E &= \frac{1}{M} \sum_{i=1}^M E_i \\
 *     S &= \left( \frac{1}{M} \frac{1}{M-1} \sum_{i=1}^M \left( E_i - E
 *          \right)^2 \right)^{-\frac{1}{2}} \\
 *     N &= \sum_{i=1}^M N_i
 * \f}
 * Note that this function weighs the result of every iteration equally,
 * independent from the sample size \f$ N_i \f$.
 */
template <typename Iterator>
inline mc_result<get_T<Iterator>> cumulative_result1(
	Iterator begin, Iterator end
) {
	typedef get_T<Iterator> T;

	std::size_t const m = std::distance(begin, end);

	switch (m)
	{
	case 0:
		return mc_result<T>(0, T(), T());

	case 1:
		return *begin;
	}

	std::size_t calls = 0;
	T sum = T();
	T sum_of_squares = T();

	for (Iterator i = begin; i != end; ++i)
	{
		T const tmp = i->value();
		calls += i->calls();
		sum += tmp;
		sum_of_squares += tmp * tmp;
	}

	T const value = sum / m;
	T const error = std::sqrt((sum_of_squares / m - value * value) / T(m - 1));

	return create_result(calls, value, error);
}

/// Returns an approximation for the \f$ \chi^2 \f$ per degree of freedom using
/// the results \f$ (E_i, S_i) \f$ pointed to by the range [`begin`, `end`). The
/// cumulative value \f$ E \f$ is given by the parameter `result`. The \f$
/// \chi^2 \f$ is then computed as:
/// \f[
///     \chi^2 / \mathrm{dof} \approx \frac{1}{n-1} \sum_{i=1}^n \frac{\left(
///     E_i - E \right)^2}{S_i^2}
/// \f]
/// If the range [`begin`, `end`) is empty, the result is zero. If it contains
/// one element the result is infinity.
template <typename Iterator, typename T>
inline T chi_square_dof(
	Iterator begin,
	Iterator end,
	mc_result<T> const& result
) {
	T sum = T();
	std::size_t n = 0;

	if (std::distance(begin, end) == 1)
	{
		return std::numeric_limits<T>::infinity();
	}

	for (Iterator i = begin; i != end; ++i)
	{
		T const tmp = i->value() - result.value();
		sum += tmp * tmp / i->variance();
		++n;
	}

	return sum / T(n - 1);
}

/// Wrapper function for \ref chi_square_dof with parameter `result` computed by
/// \ref cumulative_result0.
template <typename Iterator>
inline get_T<Iterator> chi_square_dof0(Iterator begin, Iterator end)
{
	return chi_square_dof(begin, end, cumulative_result0(begin, end));
}

/// Wrapper function for \ref chi_square_dof with parameter `result` computed by
/// \ref cumulative_result1.
template <typename Iterator>
inline get_T<Iterator> chi_square_dof1(Iterator begin, Iterator end)
{
	return chi_square_dof(begin, end, cumulative_result1(begin, end));
}

/// @}

}

#endif
