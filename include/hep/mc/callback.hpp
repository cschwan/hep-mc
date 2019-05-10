#ifndef HEP_MC_CALLBACK_HPP
#define HEP_MC_CALLBACK_HPP

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

#include "hep/mc/chkpt.hpp"
#include "hep/mc/mc_helper.hpp"
#include "hep/mc/multi_channel_chkpt.hpp"
#include "hep/mc/multi_channel_summary.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <type_traits>

namespace hep
{

/// \addtogroup callbacks
/// @{

///
enum class callback_mode
{
    /// Do not print any message after an iteration.
    silent,

    /// Print a detailed message after each iteration.
    verbose,

    /// Same as \ref callback_mode::verbose, but also writes a checkpoint to disk.
    verbose_and_write_chkpt
};

template <typename Checkpoint>
class callback
{
public:
    /// Constructor.
    callback(callback_mode mode = callback_mode::verbose, std::string const& filename = "")
        : mode_{mode}
        , filename_{filename}
    {
    }

    /// Callback function that prints a detailed summary about every iteration performed so far.
    /// This function always returns `true`.
    bool operator()(Checkpoint const& chkpt)
    {
        if (mode_ == callback_mode::silent)
        {
            return true;
        }

        using std::fabs;
        using T = typename Checkpoint::result_type::numeric_type;

        auto const& results = chkpt.results();

        std::cout << "iteration " << (results.size() - 1) << " finished.\n";

//        if constexpr (std::is_base_of_v<multi_channel_chkpt<T>, Checkpoint>)
        if (std::is_base_of<multi_channel_chkpt<T>, Checkpoint>::value)
        {
//            multi_channel_summary(chkpt, std::cout);
            multi_channel_summary(dynamic_cast <multi_channel_chkpt<T> const&> (chkpt), std::cout);
        }

        // print result for this iteration

        std::size_t const nnf = results.back().non_zero_calls() - results.back().finite_calls();
        std::size_t const num = results.back().calls();
        T const val = results.back().value();
        T const err = results.back().error();
        T const eff = T(100.0) * T(results.back().non_zero_calls()) / T(num);
        T const rel_err = T(100.0) * err / fabs(val);

        std::cout << "this iteration: N=" << num << " E=" << val << " +- " << err << " (" << rel_err
            << "%) eff=" << eff << "% nnf=" << nnf << '\n';

        // print result for all iterations

        auto const result = accumulate<weighted_with_variance>(results.begin(), results.end());

        std::size_t const num_all = result.calls();
        T const val_all = result.value();
        T const err_all = result.error();
        T const rel_err_all = (T(100.0) * result.error() / fabs(result.value()));
        T const chi = chi_square_dof<weighted_with_variance>(results.begin(), results.end());

        std::cout << "all iterations: N=" << num_all << " E=" << val_all << " +- " << err_all
            << " (" << rel_err_all << "%) chi^2/dof=" << chi << '\n' << std::endl;

        if (mode_ == callback_mode::verbose_and_write_chkpt)
        {
            std::ofstream out(filename_);
            chkpt.serialize(out);
        }

        return true;
    }

private:
    callback_mode mode_;
    std::string filename_;
};

/// @}

}

#endif
