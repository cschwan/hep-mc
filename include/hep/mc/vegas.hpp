#ifndef HEP_MC_VEGAS_HPP
#define HEP_MC_VEGAS_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2012-2013  Christopher Schwan
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

#include <cstddef>
#include <random>
#include <vector>

#include <hep/mc/mc_point.hpp>
#include <hep/mc/mc_result.hpp>

namespace hep
{

/**
 * The result of a single \ref vegas_iteration.
 */
template <typename T>
struct vegas_iteration_result : public mc_result<T>
{
	/**
	 * Constructor.
	 */
	vegas_iteration_result(
		std::size_t samples,
		std::vector<std::vector<T>> const& grid,
		std::vector<T> const& grid_refinement_data
	)
		: mc_result<T>(
			samples,
			grid_refinement_data[grid_refinement_data.size() - 2],
			grid_refinement_data[grid_refinement_data.size() - 1]
		 )
		, grid(grid)
		, grid_refinement_data(grid_refinement_data)
	{
	}

	/**
	 * The grid used to obtain this result.
	 */
	std::vector<std::vector<T>> grid;

	/**
	 * The data that will used to refine `grid` for a subsequent iteration.
	 */
	std::vector<T> grid_refinement_data;
};

/**
 *
 */
template <typename T>
struct vegas_point : public mc_point<T>
{
	/**
	 * Constructor.
	 */
	vegas_point(
		std::size_t total_samples,
		std::vector<T> const& random_numbers,
		std::vector<std::vector<T>> const& grid
	)
		: mc_point<T>(total_samples, grid.size())
		, bin(grid.size())
	{
		std::size_t const dimensions = grid.size();
		std::size_t const bins = grid[0].size();

		for (std::size_t i = 0; i != dimensions; ++i)
		{
			// compute position of 'random' in bins, as a floating point number
			T const bin_position = random_numbers[i] * bins;

			// in which bin is 'random' (integer) ?
			std::size_t const position = bin_position;

			// compute value of grid at the previous position
			T const grid_previous =
				(position == 0) ? T() : grid[i][position - 1];

			// compute difference of grid values at 'position'
			T const difference = grid[i][position] - grid_previous;

			// TODO: explain
			this->point[i] = grid_previous +
				(bin_position - position) * difference;

			// save the index of the bin in which point lies
			bin[i] = position;

			// multiply weight for each dimension
			this->weight *= difference * bins;
		}
	}

	/**
	 * The indices that determine the bin of the sample in the grid.
	 */
	std::vector<std::size_t> bin;
};

/**
 * Refine the `grid` using `grid_refinement_data`. The process can be controlled
 * by adjusting the parameter `alpha`. This function's code is based on the
 * function `refine_grid` from the CUBA VEGAS implementation from Thomas Hahn.
 */
template <typename T>
std::vector<std::vector<T>> vegas_refine_grid(
	T alpha,
	std::vector<std::vector<T>> const& grid,
	std::vector<T> const& grid_refinement_data
) {
	std::size_t const bins = grid[0].size();
	std::size_t const dimensions = grid.size();

	std::vector<std::vector<T>> new_grid(dimensions, std::vector<T>(bins));

	for (std::size_t i = 0; i != dimensions; ++i)
	{
		std::vector<T> smoothed(
			grid_refinement_data.begin() + (i+0) * bins,
			grid_refinement_data.begin() + (i+1) * bins
		);

		// smooth the entries stored in grid for dimension i
		T previous = smoothed[0];
		T current = smoothed[1];
		smoothed[0] = T(0.5) * (previous + current);
		T norm = smoothed[0];

		for (std::size_t bin = 1; bin < bins - 1; ++bin)
		{
			T const sum = previous + current;
			previous = current;
			current = smoothed[bin];
			smoothed[bin] = (sum + current) / T(3.0);
			norm += smoothed[bin];
		}
		smoothed[bins - 1] = T(0.5) * (previous + current);
		norm += smoothed[bins - 1];

		// if norm is zero there is nothing to do here
		if (norm == T())
		{
			continue;
		}

		norm = T(1.0) / norm;

		// compute the importance function for each bin
		T average_per_bin = T();

		std::vector<T> imp(bins);

		for (std::size_t bin = 0; bin < bins; ++bin)
		{
			if (smoothed[bin] > T())
			{
				T const r = smoothed[bin] * norm;
				T const impfun = std::pow((r - T(1.0)) / std::log(r), alpha);
				average_per_bin += impfun;
				imp[bin] = impfun;
			}
		}
		average_per_bin /= bins;

		// redefine the size of each bin
		current = T();
		T this_bin = T();

		int bin = -1;

		for (std::size_t new_bin = 0; new_bin != bins - 1; ++new_bin)
		{
			while (this_bin < average_per_bin)
			{
				this_bin += imp[++bin];
				previous = current;
				current = grid[i][bin];
			}

			this_bin -= average_per_bin;
			T const delta = (current - previous) * this_bin;
			new_grid[i][new_bin] = current - delta / imp[bin];
		}

		new_grid[i][bins - 1] = T(1.0);
	}

	return new_grid;
}

/**
 * Integrates `function` over the unit-hypercube using `samples` function
 * evaluations with random numbers generated by `generator`. The `grid` is used
 * to implement importance sampling; stratified sampling is not used. The
 * dimension of the function is determined by `grid.size()`.
 *
 * If `total_samples` is larger than `samples` this means the iteration is split
 * up in multiple `vegas_iteration` calls, typically run in parallel. In this
 * case `total_samples` is the total number of samples used.
 */
template <typename T, typename F, typename R = std::mt19937>
vegas_iteration_result<T> vegas_iteration(
	std::size_t samples,
	std::size_t total_samples,
	std::vector<std::vector<T>> const& grid,
	F function,
	R&& generator = std::mt19937()
) {
	// generates random number in the range [0,1]
	std::uniform_real_distribution<T> distribution;

	T average = T();
	T averaged_squares = T();

	// for kahan summation
	T compensation = T();

	std::size_t const dimensions = grid.size();
	std::size_t const bins = grid[0].size();

	std::vector<T> grid_refinement_data(dimensions * bins + 2);
	std::vector<T> random_numbers(dimensions);

	for (std::size_t i = 0; i != samples; ++i)
	{
		for (std::size_t j = 0; j != dimensions; ++j)
		{
			random_numbers[j] = distribution(generator);
		}

		vegas_point<T> const point(total_samples, random_numbers, grid);

		// evaluate function at the specified point and multiply with its weight
		T const value = function(point) * point.weight;

		// perform kahan summation 'sum += value' - this improves precision if
		// T is single precision and many values are added
		T const y = value - compensation;
		T const t = average + y;
		compensation = (t - average) - y;
		average = t;

		T const square = value * value;

		// no kahan summation needed, because it only affects the result
		// indirectly via the grid recomputation
		averaged_squares += square;

		// save square for each bin in order to refine the grid later
		for (std::size_t j = 0; j != dimensions; ++j)
		{
			grid_refinement_data[j * bins + point.bin[j]] += square;
		}
	}

	// save 'sum' and 'sum_of_squares' by rescaling the variables
	grid_refinement_data[dimensions * bins + 0] =
		average * T(total_samples);
	grid_refinement_data[dimensions * bins + 1] =
		averaged_squares * T(total_samples * total_samples);

	return vegas_iteration_result<T>(samples, grid, grid_refinement_data);
}

/**
 * Creates an equally-spaced VEGAS grid with the specified number of `bins` for
 * a function with `dimensions` parameters. Using this grid together with
 * \ref vegas_iteration is equivalent to an integration performed by the \ref
 * plain algorithm.
 */
template <typename T>
std::vector<std::vector<T>> vegas_grid(std::size_t dimensions, std::size_t bins)
{
	std::vector<std::vector<T>> grid(dimensions, std::vector<T>(bins));

	for (std::size_t i = 0; i != dimensions; ++i)
	{
		for (std::size_t j = 0; j != bins; ++j)
		{
			grid[i][j] = T(1.0 + j) / bins;
		}
	}

	return grid;
}

/**
 * Implements the VEGAS algorithm. In particular, this function calls \ref
 * vegas_iteration for every element in `iteration_samples` which determines the
 * `samples` parameter for each iteration. After each iteration the grid is
 * refined using \ref vegas_refine_grid. The grid refinement itself can be
 * controlled by the parameter `alpha`. The number of bins of the grid is
 * specified by `bins`.
 */
template <typename T, typename F, typename R = std::mt19937>
std::vector<vegas_iteration_result<T>> vegas(
	std::size_t dimensions,
	std::vector<std::size_t> const& iteration_samples,
	F function,
	std::size_t bins = 30,
	T alpha = T(1.5),
	R&& generator = std::mt19937()
) {
	// create a fresh grid
	auto grid = vegas_grid<T>(dimensions, bins);

	// vector holding all iteration results
	std::vector<vegas_iteration_result<T>> results;
	results.reserve(iteration_samples.size());

	// perform iterations
	for (auto i = iteration_samples.begin(); i != iteration_samples.end(); ++i)
	{
		auto const result = vegas_iteration(*i, *i, grid, function, generator);
		results.push_back(result);

		grid = vegas_refine_grid(alpha, grid, result.grid_refinement_data);
	}

	return results;
}

}

#endif
