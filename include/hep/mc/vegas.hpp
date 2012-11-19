#ifndef HEP_MC_VEGAS_HPP
#define HEP_MC_VEGAS_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2012  Christopher Schwan
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

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <random>
#include <vector>

namespace
{

template <typename T>
void refine_grid(
	T alpha,
	std::vector<T>& grid,
	std::vector<T>& binned_sum_of_squares
) {
	std::size_t const bins = grid.size();

	// smooth the f^2 value stored for each bin
	T previous = binned_sum_of_squares[0];
	T current = binned_sum_of_squares[1];
	binned_sum_of_squares[0] = T(0.5) * (previous + current);
	T norm = binned_sum_of_squares[0];

	for (std::size_t bin = 1; bin < bins - 1; ++bin)
	{
		T const sum = previous + current;
		previous = current;
		current = binned_sum_of_squares[bin + 1];
		binned_sum_of_squares[bin] = (sum + current) / T(3.0);
		norm += binned_sum_of_squares[bin];
	}
	binned_sum_of_squares[bins - 1] = T(0.5) * (previous + current);
	norm += binned_sum_of_squares[bins - 1];

	// if norm is zero there is nothing to do here
	if (norm == T())
	{
		return;
	}

	norm = T(1.0) / norm;

	// compute the importance function for each bin
	T average_per_bin = T();

	std::vector<T> imp(bins);

	for (std::size_t bin = 0; bin < bins; ++bin)
	{
		if (binned_sum_of_squares[bin] > T())
		{
			T const r = binned_sum_of_squares[bin] * norm;
			T impfun = std::pow((r - T(1.0)) / std::log(r), alpha);
			average_per_bin += impfun;
			imp[bin] = impfun;
		}
	}
	average_per_bin /= bins;

	// redefine the size of each bin
	current = T();
	T this_bin = T();

	std::vector<T> newgrid(bins);
	int bin = -1;

	for (std::size_t new_bin = 0; new_bin != bins - 1; ++new_bin)
	{
		while (this_bin < average_per_bin)
		{
			this_bin += imp[++bin];
			previous = current;
			current = grid[bin];
		}

		this_bin -= average_per_bin;
		T const delta = (current - previous) * this_bin;
		newgrid[new_bin] = current - delta / imp[bin];
	}

	newgrid[bins - 1] = T(1.0);
	grid = newgrid;
}

}

namespace hep
{

template <typename T>
struct vegas_sample
{
	vegas_sample(
		std::size_t dimensions,
		std::size_t steps,
		T const& random,
		std::vector<std::vector<T>> const& grid
	)
		: point(dimensions)
		, bin(dimensions)
		, weight(1.0 / steps)
	{
		std::size_t const bins = grid[0].size();

		for (std::size_t i = 0; i != dimensions; ++i)
		{
			// compute position of 'random' in bins, as a floating point number
			T const bin_position = random * bins;

			// in which bin is 'random' (integer) ?
			std::size_t const position = bin_position;

			// compute value of grid at the previous position
			T const grid_previous =
				(position == 0) ? T() : grid[i][position - 1];

			// compute difference of grid values at 'position'
			T const difference = grid[i][position] - grid_previous;

			// TODO: explain
			point[i] = grid_previous + (bin_position - position) * difference;

			// save the index of the bin in which point lies
			bin[i] = position;

			// multiply weight for each dimension
			weight *= difference * bins;
		}
	}

	std::vector<T> point;
	std::vector<std::size_t> bin;
	T weight;
};

/**
 * The result of a VEGAS Monte Carlo integration.
 */
template <typename T>
class vegas_result
{
public:
	/**
	 *
	 */
	vegas_result(std::vector<std::size_t> const& steps)
		: values(steps.size())
		, errors(steps.size())
		, steps(steps)
		, sum_of_inv_variances()
		, sum_of_averages()
	{
		values.clear();
		errors.clear();
	}

	/**
	 *
	 */
	void add_iteration(
		std::size_t samples,
		T const& sum,
		T const& sum_of_squares
	) {
		T const tmp = std::sqrt(sum_of_squares * samples);

		// compute inverse variance
		T const inv_variance = T(samples - 1.0) /
			std::max((tmp + sum) * (tmp - sum), T(0x1p-104));

		sum_of_inv_variances += inv_variance;
		T const variance = T(1.0) / sum_of_inv_variances;
		sum_of_averages += inv_variance * sum;
		T const average = variance * sum_of_averages;

		values.push_back(average);
		errors.push_back(std::sqrt(variance));
	}

	/**
	 * The computed approximations for the integral for each iteration.
	 */
	std::vector<T> values;

	/**
	 * The errors for each iteration.
	 */
	std::vector<T> errors;

	/**
	 * A std::vector containing the number of steps performed by VEGAS for each
	 * iteration.
	 */
	std::vector<std::size_t> steps;

private:
	T sum_of_inv_variances;
	T sum_of_averages;
};

/**
 * VEGAS Monte Carlo integrator. \c vegas integrates the specified \c function
 * with the specified \c dimensions using <tt>steps.size()</tt> iterations
 * each  having the number integrand evaluation specified in \c steps. The grid
 * has a resolution specified with \c bins in each dimension.
 */
template <typename T, typename F, typename A, typename R = std::mt19937>
vegas_result<T> vegas(
	std::size_t dimensions,
	std::vector<std::size_t> steps,
	std::size_t batch_size,
	std::size_t bins,
	F& function,
	A const& aux_variable,
	std::size_t seed = 0,
	R&& generator = std::mt19937()
) {
	// set parameters to one if they are zero
	bins = (bins == 0) ? 1 : bins;
	batch_size = (batch_size == 0) ? 1 : batch_size;

	// seed random number generator
	generator.seed(seed);

	// distribution [0, 1] for the random number generator
	std::uniform_real_distribution<T> distribution;

	// initialize grid
	std::vector<std::vector<T>> grid(dimensions, std::vector<T>(bins));

	for (std::size_t i = 0; i != dimensions; ++i)
	{
		for (std::size_t j = 0; j != bins; ++j)
		{
			grid[i][j] = T(1.0 + j) / bins;
		}
	}

	// container holding samples
	std::vector<vegas_sample<T>> samples;
	samples.reserve(batch_size);

	// container holding weighted function values
	std::vector<T> values;
	values.reserve(batch_size);

	std::vector<std::vector<T>> binned_sum_of_squares(dimensions,
		std::vector<T>(bins));

	//
	vegas_result<T> result{steps};

	// loop over iterations
	for (std::size_t i : steps)
	{
		std::size_t samples_done = 0;
		std::size_t size = std::min(i - samples_done, batch_size);

		// initialize to zero before each iteration
		T sum = T();
		T sum_of_squares = T();
		T compensation = T();

		// reset values to zero
		for (std::size_t j = 0; j != dimensions; ++j)
		{
			for (std::size_t k = 0; k != bins; ++k)
			{
				binned_sum_of_squares[j][k] = T();
			}
		}

		do
		{
			// fill samples
			for (std::size_t j = 0; j != size; ++j)
			{
				samples.push_back(vegas_sample<T>(
					dimensions, i, distribution(generator), grid)
				);
			}

			// generate corresponding (weighted) function values
			for (vegas_sample<T> const& sample : samples)
			{
				values.push_back(
					function(sample, aux_variable) * sample.weight
				);
			}

			// compute sum, sum of squares, and sum of squares for each bin
			for (std::size_t j = 0; j != size; ++j)
			{
				// perform kahan summation 'sum += values[j]' - this improves
				// precision if T is single precision and 'steps' is large
				T y = values[j] - compensation;
				T t = sum + y;
				compensation = (t - sum) - y;
				sum = t;

				T const square = values[j] * values[j];

				// no kahan summation needed
				sum_of_squares += square;

				// save sum_of_squares for each bin to later refine the grid
				for (std::size_t k = 0; k != dimensions; ++k)
				{
					binned_sum_of_squares[k][samples[j].bin[k]] += square;
				}
			}

			samples_done += size;

			// clear samples
			samples.clear();
		}
		while (samples_done != i);

		// add result for this iteration
		result.add_iteration(i, sum, sum_of_squares);

		// refine grid for every dimension
		for (std::size_t j = 0; j != dimensions; ++j)
		{
			refine_grid(T(1.5), grid[j], binned_sum_of_squares[j]);
		}
	}

	return result;
}

}

#endif
