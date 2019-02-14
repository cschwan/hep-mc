#ifndef HEP_PS_PLAIN_CHKPT_HPP
#define HEP_PS_PLAIN_CHKPT_HPP

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
#include "hep/mc/plain_result.hpp"

#include <random>

namespace hep
{

///
template <typename T>
using plain_chkpt = chkpt<plain_result<T>>;

///
template <typename RandomNumberEngine, typename T>
using plain_chkpt_with_rng = chkpt_with_rng<RandomNumberEngine, plain_chkpt<T>>;

///
template <typename T>
using plain_chkpt_with_mt19937 = chkpt_with_rng<std::mt19937, plain_chkpt<T>>;

}

#endif
