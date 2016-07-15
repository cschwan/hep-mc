#ifndef HEP_MC_PROJECTOR_HPP
#define HEP_MC_PROJECTOR_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2016  Christopher Schwan
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

#include "hep/mc/internal/accumulator_fwd.hpp"
#include "hep/mc/distribution_parameters.hpp"
#include "hep/mc/mc_point.hpp"

#include <cstddef>

namespace hep
{

/// \addtogroup distributions
/// @{

/// Interface for generating differential distributions.
template <typename T>
class projector
{
public:
	/// \cond DOXYGEN_IGNORE
	projector(accumulator<T, true>& accumulator, mc_point<T> const& point)
		: accumulator_(accumulator)
		, point_(point)
	{
	}
	// \endcond

	/// Projects a point for the distribution with the specified `index` to the
	/// x-axis at the value `projection` and sets the value to `value`. The
	/// weight of the point is automatically multiplied with `value`. Note that
	/// index must be smaller than the size of the vector returned by \ref
	/// parameters.
	void add(std::size_t index, T projection, T value);

	/// Returns the parameters of the distributions.
	std::vector<distribution_parameters<T>> const& parameters() const;

private:
	accumulator<T, true>& accumulator_;
	mc_point<T> const& point_;
};

/// @}

}

// contains definition of `add` because of circular dependency on `accumulator`
#include "hep/mc/internal/accumulator.hpp"

#endif
