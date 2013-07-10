#ifndef HEP_MC_HPP
#define HEP_MC_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2012-2013  Christopher Schwan
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

#include <hep/mc/linear_grid.hpp>
#include <hep/mc/mc_helper.hpp>
#include <hep/mc/mc_point.hpp>
#include <hep/mc/mc_result.hpp>
#include <hep/mc/plain.hpp>
#include <hep/mc/vegas.hpp>

namespace hep
{

/**
 * \mainpage Template Library for MC Integration
 *
 * Description
 * ===========
 *
 * `hep-mc` is a C++11 template library providing Monte Carlo integration
 * algorithms, currently only PLAIN and VEGAS.
 *
 * The design of library loosely follows the style that is used in the C++
 * standard library.
 *
 * See the README file on how to install this library.
 *
 * Example
 * =======
 *
 * The following example illustrates how to integrate the square-function using
 * VEGAS:
 * \include vegas_example.cpp
 */

/**
 * \defgroup integrands Integrand Functions
 *
 * Functions that are integrated by the Monte Carlo algorithms must accept a
 * single parameter of the type \ref mc_point and return a value of type `T`.
 * For example, an integrand of the form
 * \f[
 * f(x) := x^2
 * \f]
 * for double precision numbers would look like:
 * \code
 * double square(hep::mc_point<double> const& x)
 * {
 *     return x.point[0] * x.point[0];
 * }
 * \endcode
 * Some integration algorithms, e.g. \ref vegas(), supply additional information
 * that can be accessed by capturing the argument with different type, e.g. for
 * VEGAS:
 * \code
 * double square(hep::vegas_point<double> const& x)
 * {
 *     return x.point[0] * x.point[0];
 * }
 * \endcode
 *
 * If additional variables in the integrand are necessary, they can be supplied
 * by using functors:
 * \code
 * struct integrand
 * {
 *     double exponent;
 *
 *     integrand(double exponent)
 *         : exponent(exponent)
 *     {
 *     }
 *
 *     double operator()(hep::mc_point<double> const& x)
 *     {
 *         return std::pow(x.point[0], exponent);
 *     }
 * };
 * \endcode
 * which would be called using \ref plain() by
 * \code
 * hep::mc_result<double> result = hep::plain<double>(1, 1000, integrand(2.0));
 * \endcode
 */

/**
 * \example helper_example.cpp
 * \example read_linear_grid.cpp
 * \example vegas_example.cpp
 * \example vegas_grid.cpp
 */

}

#endif
