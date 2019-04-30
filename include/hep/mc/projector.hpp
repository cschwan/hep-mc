#ifndef HEP_MC_PROJECTOR_HPP
#define HEP_MC_PROJECTOR_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2016-2018  Christopher Schwan
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

#include "hep/mc/accumulator_fwd.hpp"
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
    projector(accumulator<T, true>* accumulator, mc_point<T> const& point)
        : accumulator_(accumulator)
        , point_(point)
    {
    }
    // \endcond

    /// This class has no copy constructor.
    projector(projector<T> const&) = delete;

    /// This class has no move constructor.
    projector(projector<T>&&) = delete;

    /// This class has no assignment operator.
    projector& operator=(projector<T> const&) = delete;

    /// This class has no move assignment operator.
    projector& operator=(projector<T>&&) = delete;

    /// Destructor.
    ~projector() = default;

    /// Adds the integrand denoted by `value` to the one-dimensional distribution with the
    /// corresponding `index` to the bin which is located at the point specified by `x`.
    void add(std::size_t index, T x, T value);

    /// Adds the integrand denoted by `value` to the two-dimensional distribution with the
    /// corresponding `index` to the bin which is located at the point specified by `x` and `y`.
    void add(std::size_t index, T x, T y, T value);

private:
    accumulator<T, true>* accumulator_;
    mc_point<T> const& point_;
};

/// @}

}

// contains definition of `add` because of circular dependency on `accumulator`
#include "hep/mc/accumulator.hpp"

#endif
