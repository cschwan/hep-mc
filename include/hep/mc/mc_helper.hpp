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

#include "hep/mc/distribution_result.hpp"
#include "hep/mc/mc_result.hpp"
#include "hep/mc/plain_result.hpp"

#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <type_traits>
#include <vector>

namespace
{

template <typename Iterator>
using hep_numeric_type =
	typename std::iterator_traits<Iterator>::value_type::numeric_type;

template <typename Iterator>
using hep_mc_result = hep::mc_result<hep_numeric_type<Iterator>>;

}

namespace hep
{

/// \addtogroup results
/// @{

/**
 * Accumulates results weighted with their inverse variance. If \f$ E_i, S_i,
 * N_i \f$ are the estimate, error and number of calls of each iteration in the
 * range defined by [`begin`, `end`) and \f$ M \f$ the number of results, then
 * the cumulative result is computed as:
 * \f{align}{
 *     E &= S^2 \sum_{i=1}^M \frac{E_i}{S_i^2} \\
 *     S &= \left( \sum_{i=1}^M \frac{1}{S_i^2} \right)^{-\frac{1}{2}} \\
 *     N &= \sum_{i=1}^M N_i
 * \f}
 */
template <typename IteratorOverMcResults>
struct weighted_with_variance
{
	/// Performs the accumulation of results in the interval `[begin, end)`.
	hep_mc_result<IteratorOverMcResults> operator()(
		IteratorOverMcResults begin,
		IteratorOverMcResults end
	) const {
		using std::sqrt;
		using T = hep_numeric_type<IteratorOverMcResults>;

		std::size_t calls = 0;
		T estimate = T();
		T variance = T();

		for (IteratorOverMcResults i = begin; i != end; ++i)
		{
			T const tmp = T(1.0) / i->variance();
			calls += i->calls();
			variance += tmp;
			estimate += tmp * i->value();
		}

		variance = T(1.0) / variance;
		estimate *= variance;

		return create_result(calls, estimate, sqrt(variance));
	}
};

/**
 * Accumulates results with the same weight. Computes a cumulative result using
 * a range of results pointed to by `begin` and `end`. If \f$ E_i, S_i, N_i \f$
 * are the estimate, error and number of calls of each iteration in the range
 * defined by [`begin`, `end`) and \f$ M \f$ the number of results, then the
 * cumulative result is computed as:
 * \f{align}{
 *     E &= \frac{1}{M} \sum_{i=1}^M E_i \\
 *     S &= \left( \frac{1}{M} \frac{1}{M-1} \sum_{i=1}^M \left( E_i - E
 *          \right)^2 \right)^{-\frac{1}{2}} \\
 *     N &= \sum_{i=1}^M N_i
 * \f}
 * Note that this function weighs the result of every iteration equally,
 * independent from the sample size \f$ N_i \f$.
 */
template <typename IteratorOverMcResults>
struct weighted_equally
{
	/// Performs the accumulation of results in the interval `[begin, end)`.
	hep_mc_result<IteratorOverMcResults> operator()(
		IteratorOverMcResults begin,
		IteratorOverMcResults end
	) const {
		using std::sqrt;
		using T = hep_numeric_type<IteratorOverMcResults>;

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

		for (IteratorOverMcResults i = begin; i != end; ++i)
		{
			T const tmp = i->value();
			calls += i->calls();
			sum += tmp;
			sum_of_squares += tmp * tmp;
		}

		T const value = sum / m;
		T const error = sqrt((sum_of_squares / m - value * value) / T(m - 1));

		return create_result(calls, value, error);
	}
};

/// Accumulates the results in the interval [`begin`, `end`) using an instance
/// of the type `Accumulator`. This type can be \ref weighted_with_variance,
/// \ref weighted_equally, or a similar type.
template <template <typename> class Accumulator,
	typename IteratorOverMcResults>
inline hep_mc_result<IteratorOverMcResults> accumulate(
	IteratorOverMcResults begin,
	IteratorOverMcResults end
) {
	return Accumulator<IteratorOverMcResults>{}(begin, end);
}

/// Returns an approximation for the \f$ \chi^2 \f$ per degree of freedom using
/// the results \f$ (E_i, S_i) \f$ pointed to by the range [`begin`, `end`). The
/// cumulative value \f$ E \f$ is calculated using an instance of `Accumulator.
/// The \f$ \chi^2 \f$ is then computed as:
/// \f[
///     \chi^2 / \mathrm{dof} \approx \frac{1}{n-1} \sum_{i=1}^n \frac{\left(
///     E_i - E \right)^2}{S_i^2}
/// \f]
/// If the range [`begin`, `end`) is empty, the result is zero. If it contains
/// one element the result is infinity.
template <template <typename> class Accumulator, typename IteratorOverMcResults>
inline hep_numeric_type<IteratorOverMcResults> chi_square_dof(
	IteratorOverMcResults begin,
	IteratorOverMcResults end
) {
	using T = hep_numeric_type<IteratorOverMcResults>;

	T sum = T();
	std::size_t n = 0;

	if (std::distance(begin, end) == 1)
	{
		return std::numeric_limits<T>::infinity();
	}

	auto const result = Accumulator<IteratorOverMcResults>{}(begin, end);

	for (IteratorOverMcResults i = begin; i != end; ++i)
	{
		T const tmp = i->value() - result.value();
		sum += tmp * tmp / i->variance();
		++n;
	}

	return sum / T(n - 1);
}

/// @}

}

#endif
