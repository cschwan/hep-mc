#ifndef HEP_MC_IMPL_VEGAS_HPP
#define HEP_MC_IMPL_VEGAS_HPP

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

#include <hep/mc/vegas.hpp>
#include <hep/mc/vegas_sample.hpp>

#include <algorithm>

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

// 		if (true)
// 		{
			newgrid[new_bin] = current - delta / imp[bin];
// 		}
// 		else
// 		{
// 			std::size_t const previous_bin = (bin == 0) ? 0 : bin - 1;
// 
// 			T new_current = std::max(
// 				new_current,
// 				current - T(2.0) * delta / (imp[bin] + imp[previous_bin])
// 			);
// 
// 			newgrid[new_bin] = new_current;
// 		}
	}

	newgrid[bins - 1] = T(1.0);
	grid = newgrid;
}

}

namespace hep
{

template <typename T, typename F, typename A, typename R>
vegas_result<T> vegas(
	std::size_t dimensions,
	std::vector<std::size_t> steps,
	std::size_t batch_size,
	std::size_t bins,
	F& function,
	A const& aux_variable,
	std::size_t seed,
	R&& generator
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
				samples.push_back(
					vegas_sample<T>(dimensions, i,
						distribution(generator), grid)
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

// 	for (std::size_t i = 0; i < bins - 1; ++i)
// 	{
// 		T previous = (i == 0) ? T(0.0) : grid[0][i - 1];
// 		T current = grid[0][i];
// 
// 		std::cout << (T(1.0) / (bins * (current - previous))) << std::endl;
// 	}

	return result;
}

}

#endif
