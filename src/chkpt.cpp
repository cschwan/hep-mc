#include "read_type.hpp"
#include "operations.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char* argv[])
{
    std::ios_base::sync_with_stdio(false);

    if (argc < 3)
    {
        std::cerr << hep::operations_help_string();

        return 1;
    }

    // capture all arguments as strings, neglecting the first, second, and last argument
    std::vector<std::string> arguments(&argv[2], &argv[argc - 1]);

    try
    {
        std::string operation_name{argv[1]};
        std::ifstream file{argv[argc - 1]};
        auto const& type = hep::read_type(file);

        hep::operations().at(type)(operation_name, arguments, file);
    }
    catch (std::runtime_error const& exception)
    {
        std::cerr << "Error: " << exception.what() << '\n';

        return 1;
    }
}
