#ifndef HEP_MC_LINEAR_GRID_HPP
#define HEP_MC_LINEAR_GRID_HPP

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
 * A grid approximating a function. This grid is called linear because it uses
 * a factorizing function
 * \f[
 *     p(x) = p(x_1) p(x_2) \cdots p(x_d)
 * \f]
 * that grows linearly with the number of `dimensions` \f$ d \f$. Every function
 * \f$ p(x_i) \f$ is piecewise-constant function using `bins` subdivisions.
 */
template <typename T>
class linear_grid
{
public:
	/**
	 * Constructs a linear grid with `dimensions` times `bins` entries. The
	 * entries are set to zero.
	 */
	linear_grid(std::size_t dimensions, std::size_t bins)
		: x(dimensions, std::vector<T>(bins))
	{
	}

	/**
	 * Constructs a linear grid from a two-dimensional vector.
	 */
	linear_grid(std::vector<std::vector<T>> const& grid)
		: x(grid)
	{
	}

	/**
	 * Read operator for the elements of this grid.
	 */
	T operator()(std::size_t dimension, std::size_t bin) const
	{
		return x[dimension][bin];
	}

	/**
	 * Read/write operator for the elements of this grid.
	 */
	T& operator()(std::size_t dimension, std::size_t bin)
	{
		return x[dimension][bin];
	}

	/**
	 * The subdivisions of every dimension of this grid.
	 */
	std::size_t bins() const
	{
		return x[0].size();
	}

	/**
	 * The number of dimensions of the function this grid appromates.
	 */
	std::size_t dimensions() const
	{
		return x.size();
	}

private:
	std::vector<std::vector<T>> x;
};

/**
 * Output operator for \ref linear_grid.
 */
template <typename CharT, typename Traits, typename T>
inline std::basic_ostream<CharT, Traits>& operator<<(
	std::basic_ostream<CharT, Traits>& ostream,
	linear_grid<T> const& grid
) {
	for (std::size_t i = 0; i != grid.dimensions(); ++i)
	{
		for (std::size_t j = 0; j < grid.bins() - 1; ++j)
		{
			ostream << grid(i, j) << " ";
		}
		ostream << grid(i, grid.bins() - 1);

		if (i < grid.dimensions() - 1)
		{
			ostream << "\n";
		}
	}

	return ostream;
}

/**
 * Input operator for \ref linear_grid. Note that the linear_grid ist not
 * resized and therefore must have the correct `dimensions` and `bins` before
 * calling this function.
 */
template <typename CharT, typename Traits, typename T>
inline std::basic_istream<CharT, Traits>& operator>>(
	std::basic_istream<CharT, Traits>& istream,
	linear_grid<T>& grid
) {
	for (std::size_t i = 0; i != grid.dimensions(); ++i)
	{
		for (std::size_t j = 0; j != grid.bins(); ++j)
		{
			istream >> grid(i, j);
		}
	}

	return istream;
}

}

#endif
