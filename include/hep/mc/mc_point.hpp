#ifndef HEP_MC_MC_POINT_HPP
#define HEP_MC_MC_POINT_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2012-2018  Christopher Schwan
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
#include <vector>

namespace hep
{

/// \addtogroup integrands
/// @{

/// A random point \f$ \vec{x} \in [0,1]^d \f$ in the \f$ d \f$-dimensional hypercube with an
/// associated weight.
template <typename T>
class mc_point
{
public:
    /// Constructor. The weight is computed using the inverse of `calls`, `point` determines the
    /// location of this point in the hypercube of dimension `point.size()`.
    explicit mc_point(std::vector<T> const& point, T weight = T(1.0))
        : weight_(weight)
        , point_(point)
    {
    }

    /// There is no copy constructor.
    mc_point(mc_point<T> const&) = delete;

    /// There is no move constructor.
    mc_point(mc_point<T>&&) = delete;

    /// There is no copy assignment operator.
    mc_point& operator=(mc_point<T> const&) = delete;

    /// There is no move assignment operator.
    mc_point& operator=(mc_point<T>&&) = delete;

    /// Destructor.
    virtual ~mc_point() = default;

    /// The weight \f$ w \f$ of this point. The PLAIN integrator (\ref plain) produces points that
    /// have weight equals one, i.e. are constant over the entire unit-hypercube. This also means
    /// that weight does not include the averaging factor \f$ 1 / N \f$ that is used to produce the
    /// expected value.
    virtual T weight() const
    {
        return weight_;
    }

    /// The coordinates of this point of the hypercube. The dimension can be obtained using
    /// `point.size()`.
    std::vector<T> const& point() const
    {
        return point_;
    }

protected:
    T mutable weight_;

private:
    std::vector<T> const& point_;
};

/// @}

}

#endif
