#ifndef HEP_MC_FUNCTION_VALUE_HPP
#define HEP_MC_FUNCTION_VALUE_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2016  Christopher Schwan
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

namespace hep
{

/// Captures the result of the evaluation of a function.
template <typename T>
class function_value
{
public:
	/// Constructor.
	function_value(T value)
		: value_(value)
	{
	}

	/// Returns the result of the function evaluation. Note that the function is
	/// multiplied with the corresponding weight.
	T value() const
	{
		return value_;
	}

private:
	T value_;
};

/// Captures the result of the evaluation of a function and a reference to
/// function itself.
template <typename T, typename F>
class function_value2 : public function_value<T>
{
public:
	/// Constructor.
	function_value2(F const& function, T value)
		: function_value<T>(value)
		, function_(function)
	{
	}

	/// Returns a reference to the function.
	F const& function() const
	{
		return function_;
	}

private:
	F const& function_;
};

}

#endif

