#ifndef HEP_MC_DISTRIBUTION_PROJECTOR_HPP
#define HEP_MC_DISTRIBUTION_PROJECTOR_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2015  Christopher Schwan
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

#include "hep/mc/distribution_parameters.hpp"

#include <utility>
#include <vector>

namespace hep
{

/// \addtogroup distributions
/// @{

/// Captures a projector and the corresponding parameters.
template <typename T, typename P>
class distribution_projector
{
public:
	/// Type of the projector.
	using projector_type = P;

	/// Numeric type.
	using numeric_type = T;

	/// Constructor. The projector must be a function with the following
	/// signature:
	/// \code
	/// void projector(P const& point, std::vector<T>& x);
	/// \endcode
	/// Here the type `P` is the type of the Monte Carlo point that is also
	/// passed to the integrand function, i.e. at least \ref mc_point. The type
	/// `T` is the numerical type used to perform the calculations. The
	/// parameter `x` is as large as `parameters.size()` and must be filled with
	/// the projection of `point` to the value of the x-axes of all the
	/// differential distributions.
	distribution_projector(
		P projector,
		std::vector<distribution_parameters<T>> const& parameters
	)
		: projector_(projector)
		, parameters_(parameters)
	{
	}

	/// Returns the parameters of all differential distributions.
	std::vector<distribution_parameters<T>> const& parameters() const
	{
		return parameters_;
	}

	/// Returns the projector that projects the MC point onto the values of the
	/// x-axes of all differential distributions.
	P projector() const
	{
		return projector_;
	}

private:
	P projector_;
	std::vector<distribution_parameters<T>> parameters_;
};

/// Helper function that creates a \ref distribution_projector without having to
/// specify the type of the projector function, i.e. the type `P` of
/// `projector`. The types of `parameters` must all be
/// `distribution_parameters<T>`.
template <typename T, typename P, typename... D>
inline distribution_projector<T, P> make_distribution_projector(
	P projector,
	D&&... parameters
) {
	return distribution_projector<T, P>(projector,
		std::vector<distribution_parameters<T>>{
		std::forward<D>(parameters)...});
}

/// @}

/// \cond DOXYGEN_IGNORE

struct one_bin_projector
{
};

template <typename T>
class distribution_projector<T, one_bin_projector>
{
public:
	using projector_type = one_bin_projector;

	using numeric_type = T;

	distribution_projector() = default;
};

template <typename T>
using default_projector = distribution_projector<T, one_bin_projector>;

/// \endcond

}

#endif

