#ifndef HEP_MC_PIECEWISE_CONSTANT_PDF_HPP
#define HEP_MC_PIECEWISE_CONSTANT_PDF_HPP

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

#include <cstddef>
#include <iosfwd>
#include <vector>

namespace hep
{

/**
 * A class for generating random numbers according to a piecewise constant probability distribution
 * function. The PDF can generate random numbers for any dimension by using a separable pdf, i.e.
 * for \f$ n \f$-dimensions the PDF looks like
 * \f[
 *     p \left( x_1, x_2, \ldots, x_n \right) = \prod_{i=1}^n p_i \left( x_i \right)
 * \f]
 * with \f$ n \f$ different piecewise constant one-dimensional PDFs \f$ p_i \left( x_i \right) \f$.
 * The value of the PDF is determined by the inverse size of the bins.
 */
template <typename T>
class piecewise_constant_pdf
{
public:
	/**
	 * Constructor. Constructs a piecewise constant PDF with the given `dimensions`, each dimension
	 * subdivided by given number of `bins`. The newly constructed PDF generates random numbers
	 * uniformly.
	 */
	piecewise_constant_pdf(std::size_t dimensions, std::size_t bins)
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

/// Output operator for \ref piecewise_constant_pdf.
template <typename CharT, typename Traits, typename T>
inline std::basic_ostream<CharT, Traits>& operator<<(
	std::basic_ostream<CharT, Traits>& out,
	piecewise_constant_pdf<T> const& pdf
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
 * Input operator for \ref piecewise_constant_pdf. Note that the PDF ist not resized and therefore
 * must have the correct `dimensions` and `bins` before calling this function.
 */
template <typename CharT, typename Traits, typename T>
inline std::basic_istream<CharT, Traits>& operator>>(
	std::basic_istream<CharT, Traits>& in,
	piecewise_constant_pdf<T>& pdf
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

}

#endif
