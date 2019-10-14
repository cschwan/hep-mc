#ifndef READ_TYPE_HPP
#define READ_TYPE_HPP

#include <cstddef>
#include <iosfwd>
#include <typeindex>

namespace hep
{

std::type_index read_type(std::istream& in);

}

#endif
