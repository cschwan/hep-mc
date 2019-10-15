#include "operations.hpp"
#include "stream_rng.hpp"
#include "make_chkpt.hpp"

#include "hep/mc/callback.hpp"
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
int perform_print(
    std::vector<std::string> const& arguments,
    hep::chkpt_with_rng<hep::stream_rng, Chkpt> const& chkpt
) {
    auto const io_flags = std::cout.flags();
    auto const io_precision = std::cout.precision();

    for (std::size_t i = 0; i != arguments.size(); ++i)
    {
        if (arguments.at(i) == "-p")
        {
            if (i + 1 == arguments.size())
            {
                throw std::runtime_error("argument " + arguments.at(i) + " is missing the "
                    "`precision` parameter");
            }

            unsigned long precision;

            try
            {
                precision = std::stoul(arguments.at(++i));
            }
            catch (std::invalid_argument const& exception)
            {
                throw std::runtime_error("argument " + arguments.at(i) + " could not be converted "
                    " to a number");
            }

            std::cout.precision(precision);
        }
        else if (arguments.at(i) == "-s")
        {
            std::cout.setf(std::ios_base::scientific, std::ios_base::floatfield);
        }
        else
        {
            std::cerr << "Warning: additional argument `" + arguments.at(i) + "` ignored\n";
        }
    }

    using T = typename Chkpt::result_type::numeric_type;

    hep::callback<Chkpt> callback{hep::callback_mode::verbose, "", T()};

    for (std::size_t i = 0; i != chkpt.results().size(); ++i)
    {
        auto copy = chkpt;
        copy.rollback(i + 1);
        callback(copy);
    }

    std::cout.precision(io_precision);
    std::cout.flags(io_flags);

    return 0;
}

template <typename Chkpt>
int perform_dists(
    std::vector<std::string> const& arguments,
    hep::chkpt_with_rng<hep::stream_rng, Chkpt> const& chkpt
) {
    auto const io_flags = std::cout.flags();
    auto const io_precision = std::cout.precision();

    for (std::size_t i = 0; i != arguments.size(); ++i)
    {
        if (arguments.at(i) == "-p")
        {
            if (i + 1 == arguments.size())
            {
                throw std::runtime_error("argument " + arguments.at(i) + " is missing the "
                    "`precision` parameter");
            }

            unsigned long precision;

            try
            {
                precision = std::stoul(arguments.at(++i));
            }
            catch (std::invalid_argument const& exception)
            {
                throw std::runtime_error("argument " + arguments.at(i) + " could not be converted "
                    " to a number");
            }

            std::cout.precision(precision);
        }
        else if (arguments.at(i) == "-s")
        {
            std::cout.setf(std::ios_base::scientific, std::ios_base::floatfield);
        }
        else
        {
            std::cerr << "Warning: additional argument `" + arguments.at(i) + "` ignored\n";
        }
    }

    if (chkpt.results().empty() || chkpt.results().front().distributions().empty())
    {
        return 0;
    }

    std::size_t const distributions = chkpt.results().front().distributions().size();

    for (std::size_t i = 0; i != distributions; ++i)
    {
        auto const& parameters = chkpt.results().front().distributions().at(i).parameters();

        std::cout << "# " << parameters.name() << '\n';

        auto const& mid_points_x = hep::mid_points_x(chkpt.results().front().distributions().at(i));
        auto const& mid_points_y = hep::mid_points_y(chkpt.results().front().distributions().at(i));

        for (std::size_t j = 0; j != mid_points_x.size(); ++j)
        {
            std::cout << mid_points_x.at(j);

            if (parameters.bins_y() > 1)
            {
                std::cout << ' ' << mid_points_y.at(j);
            }

            for (auto const& result : chkpt.results())
            {
                auto const& mc_result = result.distributions().at(i).results().at(j);

                std::cout << ' ' << mc_result.value() << ' ' << mc_result.error() << ' '
                    << mc_result.calls();
            }

            std::cout << '\n';
        }

        if (i != (distributions - 1))
        {
            std::cout << '\n';
        }
    }

    std::cout.precision(io_precision);
    std::cout.flags(io_flags);

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
    else if (operation == "print")
    {
        return perform_print(arguments, chkpt);
    }
    else if (operation == "dists")
    {
        return perform_dists(arguments, chkpt);
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
        "  Commands:\n"
        "  - iter:\n"
        "      returns the number of iterations stored in this checkpoint file\n"
        "  - print:\n"
        "      Uses the standard callback function to print all results of this\n"
        "      checkpoint. This command accepts the optional parameters `-p` and\n"
        "      `-s`\n"
        "  - dists:\n"
        "      Prints all distributions collected in the checkpoint. This command\n"
        "      accepts the optional parameters `-p` and `-s`\n"
        "  Optional arguments:\n"
        "  - `-s`:\n"
        "      switches the output to the scientific format\n"
        "  - `-p <precision>`:\n"
        "      parameter to set the precision\n";
}

}
