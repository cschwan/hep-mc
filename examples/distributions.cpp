#include "hep/mc.hpp"

#include <cstddef>
#include <iostream>
#include <vector>

double x_max =  5.0;
double x_min = -5.0;

double arctan(
    hep::mc_point<double> const& point,
    hep::projector<double>& projector
) {
    double const x_range = x_max - x_min;

    double const projection = x_range * point.point()[0] + x_min;
    double const value = std::atan(projection);

    // add to the first (zeroth) distribution
    projector.add(0, projection, value);

    return value;
}

int main()
{
    // create the integrand: We want to integrate the function `arctan` which
    // has a single dimension, and furthermore we want to generate a
    // differential distribution of this function using 100 bins in the interval
    // [-5,+5]
    auto integrand = hep::make_integrand<double, double>(
        arctan,
        1,
        hep::make_dist_params(100, x_min, x_max)
    );

    // now integrate and record the differential distributions
    auto const chkpt = hep::plain(integrand, std::vector<std::size_t>(1, 1000000));
    auto const result = chkpt.results().back();

    // integral is zero
    std::cout << "integral is I = " << result.value() << " +- " << result.error() << "\n\n";

    auto const& distribution = result.distributions()[0];
    auto const& mid_points = mid_points_x(distribution);

    // print the distribution - compare this with the plot of the function
    // 'atan(x)/10.0'
    for (std::size_t i = 0; i != mid_points.size(); ++i)
    {
        std::cout << mid_points[i] << "\t"
            << distribution.results()[i].value() << "\t"
            << distribution.results()[i].error() << "\n";
    }
}
