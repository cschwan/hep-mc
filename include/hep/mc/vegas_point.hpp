#ifndef HEP_MC_VEGAS_POINT_HPP
#define HEP_MC_VEGAS_POINT_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2014-2018  Christopher Schwan
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

#include "hep/mc/mc_point.hpp"
#include "hep/mc/vegas_pdf.hpp"

#include <cstddef>
#include <vector>

namespace hep
{

/// \addtogroup integrands
/// @{

/// A point \f$ \vec{x} \in U \f$ within the unit-hypercube \f$ U = [0,1]^d \f$ with the additional
/// information in which bin(s) of a \ref vegas_pdf the point lies.
template <typename T>
class vegas_point : public mc_point<T>
{
public:
    /// Creates a new point in the unit-hypercube \f$ U \f$ using the probability distribution
    /// function `pdf` and the random numbers in `random_numbers` for a Monte Carlo iteration with
    /// sample size specified by `total_calls`. For each dimension the point falls into a bin whose
    /// index is written `bin`.
    vegas_point(
        std::vector<T>& random_numbers,
        std::vector<std::size_t>& bin,
        vegas_pdf<T> const& pdf
    )
        : mc_point<T>(random_numbers, vegas_icdf(pdf, random_numbers, bin))
        , bin_(bin)
    {
    }

    /// There is no copy constructor.
    vegas_point(vegas_point<T> const&) = delete;

    /// There is no move constructor.
    vegas_point(vegas_point<T>&&) = delete;

    /// There is no copy assignment operator.
    vegas_point operator=(vegas_point<T> const&) = delete;

    /// There is no move assignment operator.
    vegas_point operator=(vegas_point<T>&&) = delete;

    /// Destructor.
    ~vegas_point() override = default;

    /// The indices that determine the bins of the point in the binned pdf.
    std::vector<std::size_t> const& bin() const
    {
        return bin_;
    }

private:
    std::vector<std::size_t> const& bin_;
};

/// @}

}

#endif
