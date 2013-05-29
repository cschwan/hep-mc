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

#include <hep/mc/default_vegas_parallelizer.hpp>

namespace
{

template <typename T>
void refine_grid(
	T alpha,
	std::vector<T>& grid,
	typename std::vector<T>::iterator bin_iterator
) {
	std::size_t const bins = grid.size();

	// smooth the squared function values stored for each bin
	T previous = *bin_iterator;
	T current = *(bin_iterator + 1);
	*bin_iterator = T(0.5) * (previous + current);
	T norm = *bin_iterator;

	for (std::size_t bin = 1; bin < bins - 1; ++bin)
	{
		T const sum = previous + current;
		previous = current;
		current = *(bin_iterator + bin + 1);
		*(bin_iterator + bin) = (sum + current) / T(3.0);
		norm += *(bin_iterator + bin);
	}
	*(bin_iterator + bins - 1) = T(0.5) * (previous + current);
	norm += *(bin_iterator + bins - 1);

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
		if (*(bin_iterator + bin) > T())
		{
			T const r = *(bin_iterator + bin) * norm;
			T const impfun = std::pow((r - T(1.0)) / std::log(r), alpha);
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

	// overwrite the old grid with the new one
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

	/**
	 *
	 */
	std::vector<T> point;

	/**
	 * The indices that determine the bin of the sample in the grid.
	 */
	std::vector<std::size_t> bin;

	/**
	 *
	 */
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
 * An implementation of the VEGAS Monte Carlo integration algorithm using
 * importance sampling only.
 *
 * \c vegas integrates \c function over the hypercube \f$ [0,1]^d \f$ with
 * \c d given by \c dimensions. The function must take a single argument of the
 * type \ref vegas_sample and return the function value of type \c T at the
 * sampled point:
 * \code
 * double square(hep::vegas_sample const& sample)
 * {
 *     return sample.point[0] * sample.point[0];
 * }
 * \endcode
 * The number of iterations and the samples for each iteration is specified
 * within \c iteration_samples. The grid being used to implement importance
 * sampling has a resolution specified by \c bins for every dimension. Note that
 * this parameter may have a huge impact on the precision of every iteration's
 * result. Settings this parameter too small may yield unprecise results, e.g.
 * setting it to \c 1 makes this algorithm equivalent to the \ref plain
 * algorithm.
 *
 * The parameter \c batch_size does not influence the result of \c vegas, but
 * can be used to optimize the performance when many samples are drawn. Each
 * iteration is divided into subiterations of where \c batch_size samples are
 * drawn at once.
 *
 * With \c generator you may specify a different random number generator than
 * the default mersenne twister.
 *
 * The \c parallelizer defines how this algorithm is parallelized, the default
 * choice does not parallelize.
 */
template <
	typename T,
	typename F,
	typename R = std::mt19937,
	typename P = default_vegas_parallelizer<T>>
vegas_result<T> vegas(
	std::size_t dimensions,
	std::vector<std::size_t> iteration_samples,
	std::size_t batch_size,
	std::size_t bins,
	F function,
	T const& alpha = T(1.5),
	R&& generator = std::mt19937(),
	P&& parallelizer = default_vegas_parallelizer<T>()
) {
	// set parameters to one if they are zero
	bins = (bins == 0) ? 1 : bins;
	batch_size = (batch_size == 0) ? 1 : batch_size;

	// distribution [0, 1] for the random number generator
	std::uniform_real_distribution<T> distribution;

	// get the number of processes
	std::size_t const world_size = parallelizer.world_size();
	// get the index of this process
	std::size_t const rank = parallelizer.rank();

	if (world_size != 1)
	{
		std::size_t r = 5 * rank;
		std::seed_seq sequence{r, r + 1, r + 2, r + 3, r + 4};
		generator.seed(sequence);
	}

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

	// store the following variables in a single vector; this is important for
	// fast MPI code
	std::vector<T> data(dimensions * bins + 2);

	// vector containing the binned sum of squared function values
	std::vector<T>& grid_info = data;

	T& sum = data[dimensions * bins + 0];
	T& sum_of_squares = data[dimensions * bins + 1];

	vegas_result<T> result(iteration_samples);

	// loop over iterations
	for (auto i = iteration_samples.begin(); i != iteration_samples.end(); ++i)
	{
		std::size_t const remainder = *i % world_size;
		std::size_t k = (*i / world_size);
		k += (rank < remainder) ? 1 : 0;
		std::size_t samples_done = 0;

		// initialize to zero before each iteration
		std::fill(data.begin(), data.end(), T());

		T compensation = T();

		do
		{
			std::size_t size = std::min(k - samples_done, batch_size);

			// generate samples and corresponding values
			for (std::size_t j = 0; j != size; ++j)
			{
				samples.push_back(vegas_sample<T>(
					dimensions, *i, distribution(generator), grid)
				);
				T const evaluation = function(samples.back());
				values.push_back(evaluation * samples.back().weight);
			}

			// compute sum, sum of squares, and sum of squares for each bin
			for (std::size_t j = 0; j != size; ++j)
			{
				// perform kahan summation 'sum += values[j]' - this improves
				// precision if T is single precision and many samples are added
				T y = values[j] - compensation;
				T t = sum + y;
				compensation = (t - sum) - y;
				sum = t;

				T const square = values[j] * values[j];

				// no kahan summation needed, because it only affects the result
				// indirectly via the grid recomputation
				sum_of_squares += square;

				// save sum_of_squares for each bin to later refine the grid
				for (std::size_t k = 0; k != dimensions; ++k)
				{
					grid_info[k * bins + samples[j].bin[k]] += square;
				}
			}

			samples_done += size;

			// clear samples
			samples.clear();
		}
		while (samples_done != k);

		// take 'data' from all processes, sum them elementwise into 'data'
		// again and share the result between all processes. we thereby achieve
		// a 'synchronized macro-parallelization'
		parallelizer.all_reduce(data);

		// add result for this iteration
		result.add_iteration(*i, sum, sum_of_squares);

		// refine grid for every dimension
		for (std::size_t j = 0; j != dimensions; ++j)
		{
			refine_grid(alpha, grid[j], grid_info.begin() + j * bins);
		}
	}

	return result;
}

}

#endif
