#ifndef HEP_MC_VEGAS_PDF_HPP
#define HEP_MC_VEGAS_PDF_HPP

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

#include "hep/mc/global_configuration.hpp"

#include <cmath>
#include <cstddef>
#include <iosfwd>
#include <vector>

namespace hep
{

/// \addtogroup vegas
/// @{

/**
 * A class for generating random numbers according to a piecewise constant probability distribution
 * function with varying bin widths. The PDF generates random numbers for any dimension by using
 * a separable pdf, i.e. for \f$ n \f$-dimensions the PDF is
 * \f[
 *     p \left( x_1, x_2, \ldots, x_n \right) = \prod_{i=1}^n p_i \left( x_i \right)
 * \f]
 * with \f$ n \f$ different piecewise constant one-dimensional PDFs \f$ p_i \left( x_i \right) \f$.
 * The value of the PDF is determined by the inverse size of the bins.
 */
template <typename T>
class vegas_pdf
{
public:
	/**
	 * Constructor. Constructs a piecewise constant PDF with the given `dimensions`, each dimension
	 * subdivided by given number of `bins`. Each bin has the same size and therefore this PDF
	 * generates random numbers uniformly.
	 */
	vegas_pdf(std::size_t dimensions, std::size_t bins)
		: x(dimensions, std::vector<T>(bins + 1))
	{
		std::vector<T> one_dimensional_grid(bins + 1);
		for (std::size_t i = 0; i != bins + 1; ++i)
		{
			one_dimensional_grid[i] = T(i) / T(bins);
		}

		for (std::size_t i = 0; i != dimensions; ++i)
		{
			x[i] = one_dimensional_grid;
		}
	}

	/// Returns the right boundary of `bin` for `dimension`.
	T operator()(std::size_t dimension, std::size_t bin) const
	{
		return x[dimension][bin + 1];
	}

	/// Returns the right boundary of `bin` for `dimension`.
	T& operator()(std::size_t dimension, std::size_t bin)
	{
		return x[dimension][bin + 1];
	}

	/// The number of bins for each dimension.
	std::size_t bins() const
	{
		return x[0].size() - 1;
	}

	/// The number of dimension of this PDF.
	std::size_t dimensions() const
	{
		return x.size();
	}

	/**
	 * Applies the inverse cumulative distribution function to `random_numbers` and updates it with
	 * the new numbers. The bin indices are written into `bin`. The number returned by this
	 * function is the corresponding weight; a weight of one means that this pdf is a uniform one.
	 */
	T icdf(std::vector<T>& random_numbers, std::vector<std::size_t>& bin) const
	{
		T weight = T(1.0);

		for (std::size_t i = 0; i != dimensions(); ++i)
		{
			T const position = random_numbers[i] * bins();

			// randomly select a bin index (integer part) ...
			std::size_t const index = position;

			// and the position inside this bin (0 -> left boundary, 1 -> right boundary)
			T const position_inside_bin = position - index;

			// save the index for later
			bin[i] = index;

			// determine the left bin boundary and its size
			T const bin_left = x[i][index];
			T const bin_size = x[i][index + 1] - bin_left;

			// apply the one-dimensional linear inverse CDF
			random_numbers[i] = bin_left + position_inside_bin * bin_size;

			// multiply weight for each dimension
			weight *= bin_size * bins();
		}

		return weight;
	}

private:
	std::vector<std::vector<T>> x;
};

/**
 * Refines the `pdf` using `data`, which must be the binned square-function values, and returns the
 * new pdf. The process can be controlled by the parameter `alpha` which is documented in the Vegas
 * publication \cite Vegas1 \cite Vegas2. This function's code is derived from Thomas Hahn's
 * `refine_grid` from the CUBA VEGAS implementation \cite Cuba.
 */
template <typename T>
inline vegas_pdf<T> vegas_refine_pdf(T alpha, vegas_pdf<T> const& pdf, std::vector<T> const& data)
{
	std::size_t const dimensions = pdf.dimensions();
	std::size_t const bins       = pdf.bins();

	vegas_pdf<T> new_pdf(dimensions, bins);
	std::vector<T> tmp(bins);

	for (std::size_t i = 0; i != dimensions; ++i)
	{
		// load the binned sum of squares into 'tmp'
		tmp.assign(data.begin() + (i+0) * bins, data.begin() + (i+1) * bins);

		// smooth the entries by averaging over the neighbor(s)
		T previous = tmp[0];
		T current = tmp[1];
		tmp[0] = T(0.5) * (previous + current);
		T norm = tmp[0];

		for (std::size_t bin = 1; bin != bins - 1; ++bin)
		{
			T const sum = previous + current;
			previous = current;
			current = tmp[bin + 1];
			tmp[bin] = (sum + current) / T(3.0);
			norm += tmp[bin];
		}
		tmp.back() = T(0.5) * (previous + current);
		norm += tmp.back();

		// if norm is zero there is nothing to do here
		if (norm == T())
		{
			continue;
		}

		// compute the importance function for each bin and overwrite it into 'tmp'
		T average_per_bin = T();

		for (std::size_t bin = 0; bin != bins; ++bin)
		{
			if (tmp[bin] != T())
			{
				T const r = tmp[bin] / norm;
				T const impfun = std::pow((r - T(1.0)) / std::log(r), alpha);
				average_per_bin += impfun;
				tmp[bin] = impfun;
			}
		}
		average_per_bin /= bins;

		T this_bin = T();
		std::size_t bin = 0;
		T new_current = T();

		// redefine the size of each bin
		for (std::size_t new_bin = 0; new_bin != bins - 1; ++new_bin)
		{
			for (; this_bin < average_per_bin; ++bin)
			{
				this_bin += tmp[bin];
			}

			T const previous = (bin != 1) ? pdf(i, bin - 2) : T();
			T const current = pdf(i, bin - 1);
			this_bin -= average_per_bin;
			T const delta = (current - previous) * this_bin;

			if (vegas_cuba_refinement())
			{
				T const average_importance = T(0.5) * (tmp[bin - 1] + tmp[bin != 1 ? bin - 2 : 0]);
				new_current = std::fmax(new_current, current - delta / average_importance);
			}
			else
			{
				new_current = current - delta / tmp[bin - 1];
			}

			new_pdf(i, new_bin) = new_current;
		}
	}

	return new_pdf;
}

/// Output operator for \ref vegas_pdf.
template <typename CharT, typename Traits, typename T>
inline std::basic_ostream<CharT, Traits>& operator<<(
	std::basic_ostream<CharT, Traits>& out,
	vegas_pdf<T> const& pdf
) {
	for (std::size_t i = 0; i != pdf.dimensions(); ++i)
	{
		for (std::size_t j = 0; j < pdf.bins() - 1; ++j)
		{
			out << pdf(i, j) << " ";
		}
		out << pdf(i, pdf.bins() - 1);

		if (i < pdf.dimensions() - 1)
		{
			out << "\n";
		}
	}

	return out;
}

/**
 * Input operator for \ref vegas_pdf. Note that the PDF is not resized and therefore must have
 * correct `dimensions` and `bins` before calling this function.
 */
template <typename CharT, typename Traits, typename T>
inline std::basic_istream<CharT, Traits>& operator>>(
	std::basic_istream<CharT, Traits>& in,
	vegas_pdf<T>& pdf
) {
	for (std::size_t i = 0; i != pdf.dimensions(); ++i)
	{
		for (std::size_t j = 0; j != pdf.bins(); ++j)
		{
			in >> pdf(i, j);
		}
	}

	return in;
}

/// @}

}

#endif
