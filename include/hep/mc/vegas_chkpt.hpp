#ifndef HEP_MC_VEGAS_CHKPT_HPP
#define HEP_MC_VEGAS_CHKPT_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2018-2019  Christopher Schwan
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

#include "hep/mc/chkpt.hpp"
#include "hep/mc/vegas_pdf.hpp"
#include "hep/mc/vegas_result.hpp"

#include <cassert>
#include <cstddef>
#include <iomanip>
#include <ios>
#include <limits>
#include <random>
#include <vector>

namespace hep
{

/// \addtogroup checkpoints
/// @{

///
template <typename T>
class vegas_chkpt : public chkpt<vegas_result<T>>
{
public:
    ///
    vegas_chkpt(std::size_t bins, T alpha)
        : alpha_{alpha}
        , bins_{bins}
    {
    }

    ///
    vegas_chkpt(vegas_pdf<T> const& pdf, T alpha)
        : alpha_{alpha}
        , pdf_{pdf}
    {
    }

    ///
    vegas_chkpt(std::istream& in)
        : chkpt<vegas_result<T>>{in}
    {
        in >> alpha_;

        if (this->results().empty())
        {
            pdf_.emplace_back(in);
        }
    }

    ///
    T alpha() const
    {
        return alpha_;
    }

    ///
    void dimensions(std::size_t dimensions)
    {
        if (pdf_.empty())
        {
            pdf_.emplace_back(dimensions, bins_);
        }

        assert( this->results().empty() ||
            (this->results().back().pdf().dimensions() == dimensions) );
    }

    ///
    vegas_pdf<T> const pdf() const
    {
        if (this->results().empty())
        {
            return pdf_.front();
        }

        return this->results().back().pdf();
    }

    ///
    void serialize(std::ostream& out) const override
    {
        chkpt<vegas_result<T>>::serialize(out);

        out << std::scientific << std::setprecision(std::numeric_limits<T>::max_digits10 - 1)
            << '\n' << alpha_;

        if (this->results().empty())
        {
            out << '\n';
            pdf_.front().serialize(out);
        }
    }

private:
    T alpha_;
    std::size_t bins_;
    std::vector<vegas_pdf<T>> pdf_;
};

///
template <typename RandomNumberEngine, typename T>
using vegas_chkpt_with_rng = chkpt_with_rng<RandomNumberEngine, vegas_chkpt<T>>;

///
template <typename T, typename RandomNumberEngine = std::mt19937>
vegas_chkpt_with_rng<RandomNumberEngine, T> make_vegas_chkpt(
    std::size_t bins = 128,
    T alpha = T(1.5),
    RandomNumberEngine const& rng = std::mt19937()
) {
    return vegas_chkpt_with_rng<RandomNumberEngine, T>{rng, bins, alpha};
}

///
template <typename T, typename RandomNumberEngine = std::mt19937>
vegas_chkpt_with_rng<RandomNumberEngine, T> make_vegas_chkpt(
    vegas_pdf<T> const& pdf,
    T alpha = T(1.5),
    RandomNumberEngine const& rng = std::mt19937()
) {
    return vegas_chkpt_with_rng<RandomNumberEngine, T>{rng, pdf, alpha};
}

///
template <typename T>
using default_vegas_chkpt = decltype (make_vegas_chkpt<T>());

/// @}

}

#endif
