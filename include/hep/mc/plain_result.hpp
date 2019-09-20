#ifndef HEP_MC_PLAIN_RESULT_HPP
#define HEP_MC_PLAIN_RESULT_HPP

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

#include "hep/mc/distribution_result.hpp"
#include "hep/mc/mc_result.hpp"

#include <cstddef>
#include <iosfwd>
#include <vector>

namespace hep
{

/// \addtogroup results
/// @{

/// Return type of the \ref plain MC integrator.
template <typename T>
class plain_result : public mc_result<T>
{
public:
    /// Constructor.
    plain_result(
        std::vector<distribution_result<T>> const& distributions,
        std::size_t calls,
        std::size_t non_zero_calls,
        std::size_t finite_calls,
        T sum,
        T sum_of_squares
    )
        : mc_result<T>(calls, non_zero_calls, finite_calls, sum, sum_of_squares)
        , distributions_(distributions)
    {
    }

    /// Deserialization constructor.
    explicit plain_result(std::istream& in)
        : mc_result<T>(in)
    {
        std::size_t size;
        in >> size;

        distributions_.reserve(size);

        for (std::size_t i = 0; i != size; ++i)
        {
            distributions_.emplace_back(in);
        }
    }

    /// Copy constructor.
    plain_result(plain_result<T> const&) = default;

    /// Move constructor.
    plain_result(plain_result<T>&&) noexcept = default;

    /// Assignment operator.
    plain_result& operator=(plain_result<T> const&) = default;

    /// Move assignment operator.
    plain_result& operator=(plain_result<T>&&) noexcept = default;

    /// Destructor.
    ~plain_result() override = default;

    /// Returns the differential distributions accumulated during the integration.
    std::vector<distribution_result<T>> const& distributions() const
    {
        return distributions_;
    }

    /// Serializes this object.
    void serialize(std::ostream& out) const override
    {
        mc_result<T>::serialize(out);
        out << '\n' << distributions_.size();

        for (auto const& distribution : distributions_)
        {
            out << '\n';
            distribution.serialize(out);
        }
    }

    static char const* result_name()
    {
        return "plain_result";
    }

private:
    std::vector<distribution_result<T>> distributions_;
};

/// @}

}

#endif

