#ifndef HEP_MC_INTERNAL_CALLBACK_HELPER_HPP
#define HEP_MC_INTERNAL_CALLBACK_HELPER_HPP

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

#include "hep/mc/plain_result.hpp"

#include <iosfwd>

namespace
{

template <typename T>
void print_distributions(std::ostream& out, hep::plain_result<T> const& result)
{
	auto const& distributions = result.distributions();

	for (std::size_t dist = 0; dist != distributions.size(); ++dist)
	{
		out << "results for differential distribution " << dist << ":\n";

		auto const& distribution = distributions[dist];

		for (std::size_t bin = 0; bin != distribution.results().size(); ++bin)
		{
			out << distribution.mid_points()[bin] << "\t"
				<< distribution.results()[bin].value() << "\t"
				<< distribution.results()[bin].error() << "\n";
		}

		out << "\n";
	}
}

}

#endif
