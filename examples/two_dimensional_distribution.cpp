#include "hep/mc.hpp"

#include <cstddef>
#include <iostream>
#include <vector>

double const max =  3.0;
double const min = -3.0;

double gauss(
    hep::mc_point<double> const& point,
    hep::projector<double>& projector
) {
    double const range = max - min;

    double const x = range * point.point()[0] + min;
    double const y = range * point.point()[1] + min;

    double const value = std::exp(-(x * x + y * y));

    // add to the first (zeroth) distribution
    projector.add(0, x, y, value);

    return value;
}

int main()
{
    // create the integrand: We want to integrate the function `gauss` which
    // has two dimensions, and furthermore we want to generate a differential
    // distribution of this function using 100x100 bins in the interval
    // [-3,+3]x[-3,+3]
    auto integrand = hep::make_integrand<double>(
        gauss,
        2,
        hep::distribution_parameters<double>{100, 100, min, max, min, max, "gauss"}
    );

    // now integrate and record the differential distributions
    auto const chkpt = hep::plain(
        integrand,
        std::vector<std::size_t>(1, 10000000),
        hep::make_plain_chkpt<double, std::mt19937>(),
        hep::callback<hep::plain_chkpt<double>>(hep::callback_mode::verbose_and_write_chkpt, "dist_chkpt")
    );

    auto const result = chkpt.results().back();

    // integral is zero
    std::cout << "integral is I = " << result.value() << " +- " << result.error() << "\n\n"
        "to view the distribution use the checkpoint viewer: `chkpt dists dist_chkpt`\n";
}
