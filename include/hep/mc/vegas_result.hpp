#ifndef HEP_MC_VEGAS_RESULT_HPP
#define HEP_MC_VEGAS_RESULT_HPP

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

#include "hep/mc/distribution_result.hpp"
#include "hep/mc/plain_result.hpp"
#include "hep/mc/vegas_pdf.hpp"

#include <cstddef>
#include <iomanip>
#include <ios>
#include <iosfwd>
#include <limits>
#include <ostream>
#include <vector>

namespace hep
{

/// \addtogroup results
/// @{

/// The result of a single \ref vegas_iteration.
template <typename T>
class vegas_result : public plain_result<T>
{
public:
    /// Constructor.
    vegas_result(
        plain_result<T> const& result,
        vegas_pdf<T> const& pdf,
        std::vector<T> const& adjustment_data
    )
        : plain_result<T>(result)
        , pdf_(pdf)
        , adjustment_data_(adjustment_data)
    {
    }

    /// Deserialization constructor.
    explicit vegas_result(std::istream& in)
        : plain_result<T>(in)
        , pdf_(in)
    {
        adjustment_data_.resize(pdf_.bins() * pdf_.dimensions());

        for (std::size_t i = 0; i != adjustment_data_.size(); ++i)
        {
            in >> adjustment_data_.at(i);
        }
    }

    /// Copy constructor.
    vegas_result(vegas_result<T> const&) = default;

    /// Move constructor.
    vegas_result(vegas_result<T>&&) noexcept = default;

    /// Assignment operator.
    vegas_result& operator=(vegas_result<T> const&) = default;

    /// Move assignment operator.
    vegas_result& operator=(vegas_result<T>&&) noexcept = default;

    /// Destructor.
    ~vegas_result() override = default;

    /// The pdf used to obtain this result.
    vegas_pdf<T> const& pdf() const
    {
        return pdf_;
    }

    /// The data used to adjust the \ref pdf for a subsequent iteration.
    std::vector<T> const& adjustment_data() const
    {
        return adjustment_data_;
    }

    /// Serializes this object.
    void serialize(std::ostream& out) const override
    {
        plain_result<T>::serialize(out);
        out << '\n';
        pdf_.serialize(out);
        out << '\n';

        for (std::size_t i = 0; i != adjustment_data_.size(); ++i)
        {
            out << std::scientific << std::setprecision(std::numeric_limits<T>::max_digits10 - 1)
                << adjustment_data_.at(i) << ' ';
        }
    }

    static char const* result_name()
    {
        return "vegas_result";
    }

private:
    vegas_pdf<T> pdf_;
    std::vector<T> adjustment_data_;
};

/// @}

}

#endif
