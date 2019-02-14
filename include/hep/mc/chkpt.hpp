#ifndef HEP_MC_CHKPT_HPP
#define HEP_MC_CHKPT_HPP

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

#include <cassert>
#include <istream>
#include <ostream>
#include <vector>

namespace hep
{

/// \addtogroup checkpoints
/// @{

/// Class representing checkpoints, which capture the state of each iteration performed by an
/// integrator.
template <typename Result>
class chkpt
{
public:
    using result_type = Result;

    /// Default constructor.
    chkpt() = default;

    /// Destructor.
    virtual ~chkpt() = default;

    /// Deserialization constructor.
    chkpt(std::istream& in)
    {
        std::size_t size = 0;
        in >> size;

        results_.reserve(size);

        for (std::size_t i = 0; i != size; ++i)
        {
            results_.emplace_back(in);
        }
    }

    /// Returns all results.
    std::vector<Result> const& results() const
    {
        return results_;
    }

    /// Serializes this object.
    virtual void serialize(std::ostream& out) const
    {
        std::size_t const size = results_.size();

        out << size;

        for (std::size_t i = 0; i != results_.size(); ++i)
        {
            out << '\n';
            results_.at(i).serialize(out);
        }

        serialize_generators(out);
    }

protected:
    /// Writes a textual representation of all the generators to `out`.
    virtual void serialize_generators(std::ostream& out) const = 0;

    std::vector<Result> results_;
};

template <typename RandomNumberEngine, typename Checkpoint>
class chkpt_with_rng : public Checkpoint
{
public:
    chkpt_with_rng(RandomNumberEngine generator = RandomNumberEngine())
        : Checkpoint()
        , generators_{generator}
    {
    }

    chkpt_with_rng(std::istream& in)
        : Checkpoint(in)
    {
        std::size_t const size = this->results().size() + 1;

        for (std::size_t i = 0; i != size; ++i)
        {
            RandomNumberEngine rne;
            in >> rne;
            generators_.push_back(rne);
        }
    }

    void add(typename Checkpoint::result_type const& result, RandomNumberEngine const& generator)
    {
        this->results_.push_back(result);
        generators_.push_back(generator);
    }

    RandomNumberEngine generator() const
    {
        return generators_.back();
    }

protected:
    void serialize_generators(std::ostream& out) const override
    {
        assert( generators_.size() == (this->results().size() + 1) );

        for (auto const& generator : generators_)
        {
            out << '\n' << generator;
        }
    }

private:
    std::vector<RandomNumberEngine> generators_;
};

/// @}

}

#endif
