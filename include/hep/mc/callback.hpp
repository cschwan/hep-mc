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

/// Enumeration determining the behaviour of \ref callback::operator()().
enum class callback_mode
{
    /// Do not print any messages after each iteration and do not save checkpoints to disk.
    silent,

    /// Do not print any messages after each iteration, but save checkpoints to disk.
    silent_and_write_chkpt,

    /// Print a detailed message after each iteration, but do not write checkpoints to disk.
    verbose,

    /// Same as \ref callback_mode::verbose, but also writes checkpoints to disk.
    verbose_and_write_chkpt
};

/// Default callback type used by all integration algorithms.
template <typename Checkpoint>
class callback
{
public:
    /// Numeric type used by the checkpoint.
    using numeric_type = typename Checkpoint::result_type::numeric_type;

    /// Constructor. The parameter `mode` determines the behaviour of \ref operator()(). If `mode`
    /// is \ref callback_mode::verbose_and_write_chkpt, then `filename` is the file the checkpoint
    /// is written to. If `target_rel_err` is strictly larger than zero, the integration is stopped
    /// if the accumulated result has a relative precision which is better than `target_rel_err`.
    callback(
        callback_mode mode = callback_mode::verbose,
        std::string const& filename = "",
        numeric_type target_rel_err = numeric_type()
    )
        : mode_{mode}
        , filename_{filename}
        , target_rel_err_{target_rel_err}
    {
    }

    /// Callback function whose behaviour is determined by `mode` given in the constructor. This
    /// function always returns `true`.
    bool operator()(Checkpoint const& chkpt)
    {
        using std::fabs;
        using T = numeric_type;

        auto const& results = chkpt.results();
        auto const result = accumulate<weighted_with_variance>(results.begin(), results.end());

        std::size_t const num_all = result.calls();
        T const val_all = result.value();
        T const err_all = result.error();
        T const rel_err_all = err_all / fabs(val_all);

        bool const perform_more_iterations = rel_err_all > target_rel_err_;

        if ((mode_ == callback_mode::verbose) || (mode_ == callback_mode::verbose_and_write_chkpt))
        {
            std::cout << "iteration " << (results.size() - 1) << " finished.\n";

//            if constexpr (std::is_base_of_v<multi_channel_chkpt<T>, Checkpoint>)
            if (std::is_base_of<multi_channel_chkpt<T>, Checkpoint>::value)
            {
//                multi_channel_summary(chkpt, std::cout);
                multi_channel_summary(dynamic_cast <multi_channel_chkpt<T> const&> (chkpt),
                    std::cout);
            }

            // print result for this iteration

            std::size_t const nnf = results.back().non_zero_calls() - results.back().finite_calls();
            std::size_t const num = results.back().calls();
            T const val = results.back().value();
            T const err = results.back().error();
            T const eff = T(100.0) * T(results.back().non_zero_calls()) / T(num);
            T const rel_err = err / fabs(val);

            std::cout << "this iteration: N=" << num << " E=" << val << " +- " << err << " ("
                << (T(100.0) * rel_err) << "%) eff=" << eff << "% nnf=" << nnf << '\n';

            // print result for all iterations

            T const chi = chi_square_dof<weighted_with_variance>(results.begin(), results.end());

            std::cout << "all iterations: N=" << num_all << " E=" << val_all << " +- " << err_all
                << " (" << (T(100.0) * rel_err_all) << "%) chi^2/dof=" << chi << '\n' << std::endl;
        }

        if ((mode_ == callback_mode::silent_and_write_chkpt) ||
            (mode_ == callback_mode::verbose_and_write_chkpt))
        {
            std::ofstream out(filename_);
            chkpt.serialize(out);
        }

        return perform_more_iterations;
    }

    /// Sets the mode of this callback function.
    void mode(callback_mode mode)
    {
        mode_ = mode;
    }

    /// Returns the mode of this callback function.
    callback_mode mode() const
    {
        return mode_;
    }

private:
    callback_mode mode_;
    std::string filename_;
    numeric_type target_rel_err_;
};

/// @}

}

#endif
