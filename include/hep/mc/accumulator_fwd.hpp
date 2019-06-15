#ifndef HEP_MC_ACCUMULATOR_FWD_HPP
#define HEP_MC_ACCUMULATOR_FWD_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2016-2019  Christopher Schwan
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

namespace hep
{

/// \cond INTERNAL

// The accumulator class is where the actual calculation of the total integral and all
// distributions, if they exist, happens. It serves three purposes:
//
// 1. If distributions should be generated, then we have to pass the integrand an additional object
//    that (indirectly) exposes the member function `add_to_distribution` of the accumulator class
//    that the user needs to define the distributions. If no distributions should be generated we
//    have a different implementation of the class, `accumulator<T, false>`, which passes only the
//    Monte Carlo point to the integrand and completely skips the potentially time-consuming
//    generation of distributions.
// 2. The accumulation is performed in a central place for all integrators.
// 3. The intermediate results (sum, sum_of_squars) for all distributions and the total integral is
//    stored in a single vector that can be easily distributed when the integration is performed on
//    multiple cores.
template <typename T, bool distributions>
class accumulator;

/// \endcond

}

#endif
