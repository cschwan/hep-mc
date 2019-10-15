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
#include <limits>
#include <stdexcept>
#include <utility>
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
    /// The type of results this class stores. See \ref results for the different result types.
    using result_type = Result;

    /// Default constructor.
    chkpt() = default;

    /// Deserialization constructor. This creates a checkpoint by reading from the stream `in`.
    explicit chkpt(std::istream& in)
    {
        // if the first line starts with a header, ignore it
        if (in.peek() == '#')
        {
            in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }

        std::size_t size = 0;
        in >> size;

        results_.reserve(size);

        for (std::size_t i = 0; i != size; ++i)
        {
            results_.emplace_back(in);
        }
    }

    /// Copy constructor.
    chkpt(chkpt<Result> const&) = default;

    /// Move constructor.
    chkpt(chkpt<Result>&&) noexcept = default;

    /// Assignment operator.
    chkpt& operator=(chkpt<Result> const&) = default;

    /// Move assignment operator.
    chkpt& operator=(chkpt<Result>&&) noexcept = default;

    /// Returns all results.
    std::vector<Result> const& results() const
    {
        return results_;
    }

    /// Destructor.
    virtual ~chkpt() = default;

    /// Resets the state of this checkpoint to the one described by `iteration`. Calling this
    /// function with `iteration` set to zero will remove all results, calling it with `iteration`
    /// set to the number of results will leave the checkpoint unchanged.
    virtual void rollback(std::size_t iteration)
    {
        if (iteration > results_.size())
        {
            throw std::out_of_range("parameter `iteration` too large");
        }

        results_.erase(results_.begin() + iteration, results_.end());
    }

    /// Serializes this object. This writes a textual representation of this class to the stream
    /// `out`.
    virtual void serialize(std::ostream& out) const
    {
        // write header containing version and numeric type info
        out << "# " << result_type::result_name() << " 1 "
            << std::numeric_limits<typename Result::numeric_type>::max_digits10 << '\n';

        std::size_t const size = results_.size();

        out << size;

        for (std::size_t i = 0; i != results_.size(); ++i)
        {
            out << '\n';
            results_.at(i).serialize(out);
        }
    }

protected:
    std::vector<Result> results_;
};

/// Class representing a checkpoint together with a random number generators which were used to
/// generate the results.
template <typename RandomNumberEngine, typename Checkpoint>
class chkpt_with_rng : public Checkpoint
{
public:
    /// Constructor.
    template <typename... Args>
    explicit chkpt_with_rng(RandomNumberEngine const& generator, Args&&... args)
        : Checkpoint(std::forward<Args>(args)...)
        , generators_{generator}
    {
    }

    /// Deserialization constructor. This creates a checkpoint by reading from the stream `in`.
    explicit chkpt_with_rng(std::istream& in)
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

    /// Adds a result and a random number generator to this checkpoint. The argument `result` must
    /// correspond to the result created with the random number generator returned previously with
    /// \ref generator. The argument `generator` will be random number generator for the next
    /// iteration.
    void add(typename Checkpoint::result_type const& result, RandomNumberEngine const& generator)
    {
        this->results_.push_back(result);
        generators_.push_back(generator);
    }

    /// Returns the random number generator which will be used in the next iteration.
    RandomNumberEngine generator() const
    {
        return generators_.back();
    }

    void rollback(std::size_t iteration) override
    {
        Checkpoint::rollback(iteration);
        generators_.erase(generators_.begin() + iteration, generators_.end());
    }

    void serialize(std::ostream& out) const override
    {
        Checkpoint::serialize(out);

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
