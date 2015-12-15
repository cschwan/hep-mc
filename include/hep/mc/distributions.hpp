#ifndef HEP_MC_DISTRIBUTIONS_HPP
#define HEP_MC_DISTRIBUTIONS_HPP

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

#include "hep/mc/distributions.hpp"
#include "hep/mc/mc_result.hpp"

#include <cstddef>
#include <type_traits>
#include <vector>

namespace hep
{

/// \addtogroup distributions
/// @{

template <typename T, typename P>
class distributions;

/// @}

/// \addtogroup internal
/// @{

struct one_bin_projector
{
};

template <typename T>
class distributions<T, one_bin_projector>
{
public:
	using projector = one_bin_projector;
	using numeric_type = T;

	distributions() = default;
};

template <typename T>
using default_distribution = distributions<T, one_bin_projector>;

/// @}

}

#endif
