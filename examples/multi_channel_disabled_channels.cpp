#include "hep/mc.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

int main()
{
    constexpr double s[] = { -10.0, +10.0, 0.0 };

    auto const function = [&](hep::multi_channel_point<double> const& x) {
        double t = x.coordinates()[0];

        double const tms0 = t - s[0];
        double const tms1 = t - s[1];
        double const tms2 = t - s[2];

        double const value0 = 2.0 * std::exp(-tms0 * tms0);
        double const value1 = 1.0 * std::exp(-tms1 * tms1);
        double const value2 = 0.1 * std::exp(-tms2 * tms2);

        return (value0 + value1 + value2) / std::sqrt(std::acos(-1.0)) / 3.1;
    };

    auto const densities = [&](
        std::size_t channel,
        std::vector<double> const& random_numbers,
        std::vector<double>& coordinates,
        std::vector<std::size_t> const& enabled_channels,
        std::vector<double>& densities,
        hep::multi_channel_map action
    ) {
        if (action == hep::multi_channel_map::calculate_densities)
        {
            // we calculated the densities already
            return 1.0;
        }

        double t = std::tan(std::acos(-1.0) * (random_numbers[0] - 0.5));

        t += s[channel];
        coordinates[0] = t;

        double const tms[] = { t - s[0], t - s[1], t - s[2] };

        // set the channel densities
        for (auto const i : enabled_channels)
        {
            densities[i] = 1.0 / (1.0 + tms[i] * tms[i]) / std::acos(-1.0);
        }

        // global weight
        return 1.0;
    };

    hep::multi_channel_callback<double>(
        hep::multi_channel_verbose_callback<double>);

    auto results = hep::multi_channel(
        hep::make_multi_channel_integrand<double>(function, 1, densities, 1, 3),
        std::vector<std::size_t>(10, 1000000)
    ).results();

    std::cout << ">>> Disabling channel #2\n";

    auto adjustment_data = results.back().adjustment_data();
    adjustment_data.at(2) = 0.0;

    auto const weights = hep::multi_channel_refine_weights(
        results.back().channel_weights(),
        adjustment_data,
        0.0,
        0.25
    );

    std::cout << ">>> Starting new run without channel #2\n\n";

    hep::multi_channel(
        hep::make_multi_channel_integrand<double>(function, 1, densities, 1, 3),
        std::vector<std::size_t>(10, 1000000),
        hep::make_multi_channel_chkpt<double>(weights)
    ).results();

    return 0;
}
