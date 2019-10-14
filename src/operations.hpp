#ifndef OPERATIONS_HPP
#define OPERATIONS_HPP

#include <functional>
#include <iosfwd>
#include <string>
#include <typeindex>
#include <unordered_map>
#include <vector>

namespace hep
{

using chkpt_operation =
    std::function<int(std::string const&, std::vector<std::string> const&, std::istream& in)>;

std::unordered_map<std::type_index, chkpt_operation> const& operations();

std::string operations_help_string();

}

#endif
