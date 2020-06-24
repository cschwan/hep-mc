#include "hep/mc.hpp"

#include <cstddef>
#include <iostream>
#include <vector>

double square(hep::mc_point<double> const& x)
{
    return 3.0 * x.point()[0] * x.point()[0];
}

int main()
{
    // print only 3 digits
    std::cout.precision(3);

    // perform 100 iterations with 1000 calls each _at maximum_
    auto results = hep::vegas(
        hep::make_integrand<double, double>(square, 1),
        std::vector<std::size_t>(100, 1000),
        hep::make_vegas_chkpt<double, double>(),
        // stop if error is better than 1% (=0.01)
        hep::callback<hep::vegas_chkpt<double, double>>(hep::callback_mode::verbose, "", 0.01)
    );

    return 0;
}
