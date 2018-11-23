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
        hep::distribution_parameters<double>{100, 100, min, max, min, max, ""}
    );

    // now integrate and record the differential distributions
    auto const results = hep::plain(integrand, std::vector<std::size_t>(1, 10000000));
    auto const result = results.back();

    // integral is zero
    std::cout << "integral is I = " << result.value() << " +- "
        << result.error() << "\n\n";

    auto const& distribution = result.distributions()[0];
    auto const& mid_points_x = hep::mid_points_x(distribution);
    auto const& mid_points_y = hep::mid_points_y(distribution);

    std::cout.setf(std::ios_base::scientific);

    // print the distribution
    for (std::size_t i = 0; i != mid_points_x.size(); ++i)
    {
        std::cout << mid_points_x[i] << '\t' << mid_points_y[i] << '\t'
            << distribution.results()[i].value() << '\t'
            << distribution.results()[i].error() << '\n';
    }
}
