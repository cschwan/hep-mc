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

/// Checkpoints created by the \ref vegas_group.
template <typename T>
class vegas_chkpt : public chkpt<vegas_result<T>>
{
public:
    /// Constructor. Creates an empty checkpoint with a uniform \ref vegas_pdf with the specified
    /// number of `bins`. The argument `alpha` will be used to adapt the grids after each iteration.
    vegas_chkpt(std::size_t bins, T alpha)
        : alpha_{alpha}
        , bins_{bins}
    {
    }

    /// Constructor. Creates an empty checkpoint with the user-defined \ref vegas_pdf which can be
    /// used to integrate using a PDF different from a uniform one.
    vegas_chkpt(vegas_pdf<T> const& pdf, T alpha)
        : alpha_{alpha}
        , pdf_{pdf}
    {
    }

    /// Deserialization constructor. This creates a checkpoint by reading from the stream `in`.
    explicit vegas_chkpt(std::istream& in)
        : chkpt<vegas_result<T>>{in}
    {
        in >> alpha_;

        if (this->results().empty())
        {
            pdf_.emplace_back(in);
        }
    }

    /// Returns the parameter `alpha`, which is used to refine the PDF of VEGAS after each
    /// iteration.
    T alpha() const
    {
        return alpha_;
    }

    /// Sets the number of dimensions which is a parameter needed for the VEGAS integration
    void dimensions(std::size_t dimensions)
    {
        if (this->results().empty() && pdf_.empty())
        {
            pdf_.emplace_back(dimensions, bins_);
        }

        assert( this->results().empty() ||
            (this->results().back().pdf().dimensions() == dimensions) );
    }

    /// Returns the PDF which is used for the next iteration.
    vegas_pdf<T> pdf() const
    {
        auto const& results = this->results();

        if (results.empty())
        {
            return pdf_.front();
        }

        return vegas_refine_pdf(results.back().pdf(), alpha_, results.back().adjustment_data());
    }

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

/// Checkpoint with random number generators created by using the \ref plain_group.
template <typename RandomNumberEngine, typename T>
using vegas_chkpt_with_rng = chkpt_with_rng<RandomNumberEngine, vegas_chkpt<T>>;

/// Helper function to create an initial checkpoint to start the \ref vegas_group.
template <typename T, typename RandomNumberEngine = std::mt19937>
vegas_chkpt_with_rng<RandomNumberEngine, T> make_vegas_chkpt(
    std::size_t bins = 128,
    T alpha = T(1.5),
    RandomNumberEngine const& rng = RandomNumberEngine()
) {
    return vegas_chkpt_with_rng<RandomNumberEngine, T>{rng, bins, alpha};
}

/// Helper function to create an initial checkpoint to start the \ref vegas_group.
template <typename T, typename RandomNumberEngine = std::mt19937>
vegas_chkpt_with_rng<RandomNumberEngine, T> make_vegas_chkpt(
    vegas_pdf<T> const& pdf,
    T alpha = T(1.5),
    RandomNumberEngine const& rng = RandomNumberEngine()
) {
    return vegas_chkpt_with_rng<RandomNumberEngine, T>{rng, pdf, alpha};
}

/// Helper function create a checkpoint reading from the stream `in`. Note the the numeric type `T`
/// as well as the type of the random number generator, `RandomNumberEngine` have to explicitly
/// stated.
template <typename T, typename RandomNumberEngine>
vegas_chkpt_with_rng<RandomNumberEngine, T> make_vegas_chkpt(std::istream& in)
{
    if (in.peek() == std::istream::traits_type::eof())
    {
        return make_vegas_chkpt<T, RandomNumberEngine>();
    }

    return vegas_chkpt_with_rng<RandomNumberEngine, T>{in};
}

/// Return type of \ref make_vegas_chkpt with default arguments.
template <typename T>
using default_vegas_chkpt = decltype (make_vegas_chkpt<T>());

/// @}

}

#endif
