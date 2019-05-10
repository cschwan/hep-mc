#ifndef HEP_MC_MULTI_CHANNEL_SUMMARY_HPP
#define HEP_MC_MULTI_CHANNEL_SUMMARY_HPP

/*
 * hep-mc - A Template Library for Monte Carlo Integration
 * Copyright (C) 2019  Christopher Schwan
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

#include "hep/mc/multi_channel_chkpt.hpp"
#include "hep/mc/multi_channel_max_difference.hpp"
#include "hep/mc/multi_channel_weight_info.hpp"

#include <cstddef>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

namespace hep
{

/// \cond INTERNAL

inline std::string make_list_of_ranges(std::vector<std::size_t> const& indices)
{
    std::ostringstream ranges;
    std::size_t a = 0;
    std::size_t b = 0;

    for (auto i = indices.begin(); i != indices.end(); ++i)
    {
        a = *i;
        b = a;

        if (i != indices.begin())
        {
            ranges << ',';
        }

        for (;;)
        {
            auto next = std::next(i);

            if ((next == indices.end()) || (*next != (*i + 1)))
            {
                break;
            }

            b = *next;
            ++i;
        }

        if (a == b)
        {
            ranges << a;
        }
        else
        {
            ranges << a << '-' << b;
        }
    }

    return ranges.str();
}

template <typename T>
inline void multi_channel_summary(multi_channel_chkpt<T> const& chkpt, std::ostream& out)
{
    auto const& results = chkpt.results();
    T const max_difference = multi_channel_max_difference(results.back());
    multi_channel_weight_info<T> info(results.back());
    std::size_t const channels = info.channels().size();
    std::size_t const min_channels = info.minimal_weight_count();

    out << "summary of a-priori weights: D=" << max_difference << " for " << channels << " channel"
        << (channels > 1 ? "s\n" : "\n");

    out << "wmin=" << info.weights().front() << " (N=" << info.calls().front() << ") in "
        << min_channels << " channel" << (min_channels > 1 ? "s" : "") << ": #"
        << make_list_of_ranges(minimal_weight_channels(info)) << '\n';

    auto weight_printer = [&](std::string const& prefix, std::size_t index) {
        out << prefix << info.weights().at(index) << " (N=" << info.calls().at(index)
            << ") in channel #" << info.channels().at(index) << '\n';
    };

    // number of non-minimal weights
    std::size_t printable_channels = channels - min_channels;

    if (printable_channels > 0)
    {
        // print maximum weight channel separately
        --printable_channels;
    }

    if (printable_channels > 0)
    {
        std::size_t const number = 5;

        if (printable_channels <= 2 * number + 1)
        {
            for (std::size_t i = 0; i != printable_channels; ++i)
            {
                weight_printer("   w=", min_channels + i);
            }
        }
        else
        {
            std::size_t offset = min_channels;

            for (std::size_t i = 0; i != number; ++i)
            {
                weight_printer("   w=", offset + i);
            }

            out << "     ...\n";

            offset = channels - number - 1;

            for (std::size_t i = 0; i != number; ++i)
            {
                weight_printer("   w=", offset + i);
            }
        }
    }

    if (min_channels != channels)
    {
        weight_printer("wmax=", channels - 1);
    }
}

/// \endcond

}

#endif
