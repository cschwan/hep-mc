#include "operations.hpp"
#include "stream_rng.hpp"
#include "make_chkpt.hpp"

#include "hep/mc/chkpt.hpp"
#include "hep/mc/multi_channel_chkpt.hpp"
#include "hep/mc/multi_channel_result.hpp"
#include "hep/mc/plain_chkpt.hpp"
#include "hep/mc/plain_result.hpp"
#include "hep/mc/vegas_chkpt.hpp"
#include "hep/mc/vegas_result.hpp"

#include <iostream>
#include <stdexcept>
#include <typeinfo>

namespace
{

template <typename Chkpt>
int perform_iter(
    std::vector<std::string> const& arguments,
    hep::chkpt_with_rng<hep::stream_rng, Chkpt> const& chkpt
) {
    if (!arguments.empty())
    {
        std::cerr << "Warning: additional arguments ignored\n";
    }

    std::cout << chkpt.results().size() << '\n';

    return 0;
}

template <typename Chkpt>
int dispatch_operations(
    std::string const& operation,
    std::vector<std::string> const& arguments,
    std::istream& in
) {
    auto const& chkpt = hep::make_chkpt<Chkpt>(in);

    if (operation == "iter")
    {
        return perform_iter(arguments, chkpt);
    }
    else
    {
        throw std::runtime_error("operation `" + operation + "`not recognized");
    }
}

}

namespace hep
{

std::unordered_map<std::type_index, chkpt_operation> const& operations()
{
    static std::unordered_map<std::type_index, chkpt_operation> operations_ = {
        { std::type_index(typeid(multi_channel_result<float>)),
             dispatch_operations<multi_channel_chkpt<float>> },
        { std::type_index(typeid(multi_channel_result<double>)),
             dispatch_operations<multi_channel_chkpt<double>> },
        { std::type_index(typeid(multi_channel_result<long double>)),
             dispatch_operations<multi_channel_chkpt<long double>> },
        { std::type_index(typeid(plain_result<float>)),
             dispatch_operations<plain_chkpt<float>> },
        { std::type_index(typeid(plain_result<double>)),
             dispatch_operations<plain_chkpt<double>> },
        { std::type_index(typeid(plain_result<long double>)),
             dispatch_operations<plain_chkpt<long double>> },
        { std::type_index(typeid(vegas_result<float>)),
             dispatch_operations<vegas_chkpt<float>> },
        { std::type_index(typeid(vegas_result<double>)),
             dispatch_operations<vegas_chkpt<double>> },
        { std::type_index(typeid(vegas_result<long double>)),
             dispatch_operations<vegas_chkpt<long double>> },
    };

    return operations_;
}

std::string operations_help_string()
{
    return /* 72 character limit /////////////////////////////////////////////// */
        "Usage: chkpt <command> [opt_args...] file\n"
        "  Options:\n"
        "  - iter: returns the number of iterations stored in this checkpoint file\n";
}

}
