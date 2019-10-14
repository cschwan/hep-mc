#include "read_type.hpp"

#include "hep/mc/mc_result.hpp"
#include "hep/mc/multi_channel_result.hpp"
#include "hep/mc/plain_result.hpp"
#include "hep/mc/vegas_result.hpp"

#include <cassert>
#include <iostream>
#include <limits>
#include <string>
#include <sstream>
#include <stdexcept>
#include <typeinfo>
#include <vector>

namespace hep
{

std::type_index read_type(std::istream& in)
{
    auto const pos = in.tellg();
    std::string line;
    std::getline(in, line);
    in.seekg(pos);
    std::istringstream stream{line};
    std::string token;
    std::vector<std::string> tokens;
    tokens.reserve(4);

    while (std::getline(stream, token, ' '))
    {
        tokens.push_back(token);
    }

    if ((tokens.size() != 4) || (tokens.at(0) != "#"))
    {
        throw std::runtime_error("no header found/invalid header");
    }

    if ((tokens.at(1) != "plain_result") && (tokens.at(1) != "multi_channel_result") &&
        (tokens.at(1) != "vegas_result"))
    {
        throw std::runtime_error("invalid result type");
    }

    unsigned long version;
    unsigned long max_digits10;

    try
    {
        version = std::stoul(tokens.at(2));
    }
    catch (...)
    {
        throw std::runtime_error("could not parse `version` field");
    }

    if (version > 1)
    {
        throw std::runtime_error("version not supported");
    }

    try
    {
        max_digits10 = std::stoul(tokens.at(3));
    }
    catch (...)
    {
        throw std::runtime_error("could parse `max_digits10` field");
    }

    constexpr auto max_digits10_float = std::numeric_limits<float>::max_digits10;
    constexpr auto max_digits10_double = std::numeric_limits<double>::max_digits10;
    constexpr auto max_digits10_long_double = std::numeric_limits<long double>::max_digits10;

    if (max_digits10 <= max_digits10_float)
    {
        if (max_digits10 != max_digits10_float)
        {
            std::cerr << "Warning: `max_digits10` is " << max_digits10 << " but `float` has "
                << max_digits10_float;
        }

        if (tokens.at(1) == "multi_channel_result")
        {
            return std::type_index(typeid(hep::multi_channel_result<float>));
        }
        else if (tokens.at(1) == "plain_result")
        {
            return std::type_index(typeid(hep::plain_result<float>));
        }
        else if (tokens.at(1) == "vegas_result")
        {
            return std::type_index(typeid(hep::vegas_result<float>));
        }
    }
    else if (max_digits10 <= max_digits10_double)
    {
        if (max_digits10 != max_digits10_double)
        {
            std::cerr << "Warning: `max_digits10` is " << max_digits10 << " but `double` has "
                << max_digits10_double;
        }

        if (tokens.at(1) == "multi_channel_result")
        {
            return std::type_index(typeid(hep::multi_channel_result<double>));
        }
        else if (tokens.at(1) == "plain_result")
        {
            return std::type_index(typeid(hep::plain_result<double>));
        }
        else if (tokens.at(1) == "vegas_result")
        {
            return std::type_index(typeid(hep::vegas_result<double>));
        }
    }
    else
    {
        if (max_digits10 != max_digits10_long_double)
        {
            std::cerr << "Warning: `max_digits10` is " << max_digits10
                << " but `long double` has " << max_digits10_long_double;
        }

        if (tokens.at(1) == "multi_channel_result")
        {
            return std::type_index(typeid(hep::multi_channel_result<long double>));
        }
        else if (tokens.at(1) == "plain_result")
        {
            return std::type_index(typeid(hep::plain_result<long double>));
        }
        else if (tokens.at(1) == "vegas_result")
        {
            return std::type_index(typeid(hep::vegas_result<long double>));
        }
    }

    // if this happens, we didn't cover all the cases
    assert( false );
}

}
