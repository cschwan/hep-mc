#ifndef HEP_MC_MC_RESULT_HPP
#define HEP_MC_MC_RESULT_HPP

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

#include <cmath>
#include <cstddef>
#include <iomanip>
#include <ios>
#include <istream>
#include <limits>
#include <ostream>
#include <vector>

namespace hep
{

/// \addtogroup results
/// @{

/**
 * The estimation of a Monte Carlo integration. Every Monte Carlo integrator returns one or more
 * instances of this class. The PLAIN Monte Carlo integrator, for example, calculates the parameters
 * as follows:
 * \f{align}{
 *     E &= \frac{1}{N} \sum_{i=1}^N f ( \vec{x}_i ) \\
 *     S^2 &= \frac{1}{N-1} \left[ \frac{1}{N} \sum_{i=1}^N f^2 ( \vec{x}_i ) - E^2  \right]
 * \f}
 */
template <typename T>
class mc_result
{
public:
    /// The numeric type used for member variables.
    using numeric_type = T;

    /// Constructor.
    mc_result(
        std::size_t calls,
        std::size_t non_zero_calls,
        std::size_t finite_calls,
        T sum,
        T sum_of_squares
    )
        : calls_(calls)
        , non_zero_calls_(non_zero_calls)
        , finite_calls_(finite_calls)
        , sum_(sum)
        , sum_of_squares_(sum_of_squares)
    {
    }

    /// Deserialization constructor.
    explicit mc_result(std::istream& in)
    {
        in >> calls_ >> non_zero_calls_ >> finite_calls_ >> sum_ >> sum_of_squares_;
    }

    /// Copy constructor.
    mc_result(mc_result<T> const&) = default;

    /// Move constructor.
    mc_result(mc_result<T>&&) noexcept = default;

    /// Assignment operator.
    mc_result& operator=(mc_result<T> const&) = default;

    /// Move assignment operator.
    mc_result& operator=(mc_result<T>&&) noexcept = default;

    /// Destructor.
    virtual ~mc_result() = default;

    /// The number of function evaluations \f$ N \f$ performed to obtain this
    /// result.
    std::size_t calls() const
    {
        return calls_;
    }

    /// Expectation value \f$ E \f$ of this result.
    T value() const
    {
        return sum_ / T(calls_);
    }

    /// Variance \f$ S^2 \f$ of the expectation value.
    T variance() const
    {
        return (sum_of_squares_ - sum_ * sum_ / T(calls_)) / T(calls_) / T(calls_ - 1);
    }

    /// Standard deviation \f$ S \f$ of the expectation value.
    T error() const
    {
        using std::sqrt;

        return sqrt(variance());
    }

    /// Returns the number integrand evaluations that were finite.
    std::size_t finite_calls() const
    {
        return finite_calls_;
    }

    /// Returns the number of integrand evaluations that were not zero. This includes the both
    /// finite and non-finite numbers.
    std::size_t non_zero_calls() const
    {
        return non_zero_calls_;
    }

    /// Returns the sum, i.e. \f$ \sum_{i=1}^N f ( \vec{x}_i ) \f$.
    T sum() const
    {
        return sum_;
    }

    /// Returns the sum of squares, i.e. \f$ \sum_{i=1}^N f^2 ( \vec{x}_i ) \f$.
    T sum_of_squares() const
    {
        return sum_of_squares_;
    }

    /// Serializes this object.
    virtual void serialize(std::ostream& out) const
    {
        out << calls_ << ' ' << non_zero_calls_ << ' ' << finite_calls_ << ' ' << std::scientific
            << std::setprecision(std::numeric_limits<T>::max_digits10 - 1) << sum_ << ' '
            << sum_of_squares_;
    }

    static char const* result_name()
    {
        return "mc_result";
    }

private:
    std::size_t calls_;
    std::size_t non_zero_calls_;
    std::size_t finite_calls_;
    T sum_;
    T sum_of_squares_;
};

/// Creates a \ref mc_result using the parameters `calls`, `non_zero_calls`, `value` and `error`.
template <typename T>
inline mc_result<T> create_result(
    std::size_t calls,
    std::size_t non_zero_calls,
    std::size_t finite_calls,
    T value,
    T error
) {
    T sum = T(calls) * value;
    T sum_of_squares = T(calls) * (value * value + T(calls - 1) * error * error);

    return mc_result<T>(calls, non_zero_calls, finite_calls, sum, sum_of_squares);
}

/// @}

}

#endif
