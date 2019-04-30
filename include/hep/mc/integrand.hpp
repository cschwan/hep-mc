#ifndef HEP_MC_INTEGRAND_HPP
#define HEP_MC_INTEGRAND_HPP

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

#include "hep/mc/distribution_parameters.hpp"

#include <cstddef>
#include <type_traits>
#include <utility>
#include <vector>

namespace hep
{

/// \addtogroup integrands
/// @{

/// Class representing a function that can be integrated using the PLAIN-like algorithms, which
/// currently are PLAIN and VEGAS.
template <typename T, typename F, bool distributions>
class integrand
{
public:
    /// Numeric type of the integrand.
    using numeric_type = T;

    /// The type of the integrand function that is integrated.
    using function_type = F;

    /// Signals whether this integrand wants to generate distributions or not.
    static constexpr bool has_distributions = distributions;

    /// Constructor. Instead of using the constructor directly you should consider using one of the
    /// helper functions \ref make_integrand.
    template <typename G>
    integrand(
        G&& function,
        std::size_t dimensions,
        std::vector<distribution_parameters<T>> const& parameters
    )
        : function_(std::forward<G>(function))
        , parameters_(parameters)
        , dimensions_(dimensions)
    {
    }

    /// Returns the dimension of the integrand.
    std::size_t dimensions() const
    {
        return dimensions_;
    }

    /// Returns the integrand's function.
    F& function()
    {
        return function_;
    }

    /// Returns the parameters of the distribution(s).
    std::vector<distribution_parameters<T>> const& parameters() const
    {
        return parameters_;
    }

private:
    F function_;
    std::vector<distribution_parameters<T>> parameters_;
    std::size_t dimensions_;
};

/// Template alias for an \ref integrand with its type `F` decayed with `std::decay`.
template <typename T, typename F, bool distributions>
using integrand_type = integrand<T, typename std::decay<F>::type, distributions>;

/// PLAIN/VEGAS integrand constructor. This function constructs an \ref integrand using the given
/// `function` that must accept points from the \f$ d \f$-dimensional hypercube, where \f$ d \f$ is
/// given by the parameter `dimensions`. The type of the point is determined by the integration
/// algorithm later used on the integrand, e.g. for \ref vegas it is \ref vegas_point. For this case
/// the integrand would look like:
/// \code
/// T function(hep::vegas_point<T> const& x)
/// {
///     // x.point.size() is as large as specified with `dimensions`
///     return /* calculate the function value from x.point */;
/// }
/// \endcode
template <typename T, typename F>
inline integrand_type<T, F, false> make_integrand(F&& function, std::size_t dimensions)
{
    return integrand_type<T, F, false>(
        std::forward<F>(function),
        dimensions,
        std::vector<distribution_parameters<T>>()
    );
}

/// PLAIN/VEGAS distributions constructor. This function constructs an \ref integrand using the
/// given `function` that must accept points from the \f$ d \f$-dimensional hypercube and a
/// reference to a \ref projector that generates the distributions. The dimension \f$ d \f$ is given
/// by the parameter `dimension`, and `parameters` define the number and parameters of the
/// distribution(s). For the VEGAS algorithm the function would look like:
/// \code
/// T function(hep::vegas_point<T> const& x, hep::projector<T>& projector)
/// {
///     T const x0 = /* calculate the position of `x` for the distribution 0 */;
///     T const f = /* calculate the function value from x.point */;
///
///     // add the function value `f` to the zeroeth distribution at `x0`
///     projector.add(0, x0, f);
///
///     // return the value of the function for the integrator
///     return f;
/// }
/// \endcode
template <typename T, typename F, typename... D>
inline integrand_type<T, F, true> make_integrand(
    F&& function,
    std::size_t dimensions,
    D&&... parameters
) {
    return integrand_type<T, F, true>(
        std::forward<F>(function),
        dimensions,
        std::vector<distribution_parameters<T>>{std::forward<D>(parameters)...}
    );
}

/// Shortcut for accessing the numeric type of an integrand that is possibly a
/// reference.
template <typename I>
using numeric_type_of = typename std::remove_reference<I>::type::numeric_type;

/// @}

}

#endif
