#ifndef USE_MPI
#include "hep/mc.hpp"
#else
#include "hep/mc-mpi.hpp"
#include <mpi.h>
#endif

#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

int main(int argc, char* argv[])
{
#ifdef USE_MPI
	MPI_Init(&argc, &argv);
#else
	// suppress warning about parameters being unused
	if (argv[argc]) {}
#endif

	constexpr double s0 = -10.0;
	constexpr double s1 = +10.0;

	auto const function = [](hep::multi_channel_point<double> const& x) {
		double s = x.coordinates()[0];

		double const sms0 = s - s0;
		double const sms1 = s - s1;

		return (2.0 * std::exp(-sms0 * sms0) + std::exp(-sms1 * sms1)) 
			/ std::sqrt(std::acos(-1.0)) / 3.0;
	};

	auto const densities = [](
		std::size_t channel,
		std::vector<double> const& random_numbers,
		std::vector<double>& coordinates,
		std::vector<double>& densities,
		hep::multi_channel_map action
	) {
		if (action == hep::multi_channel_map::calculate_densities)
		{
			// we calculated the densities already
			return 1.0;
		}

		double s = std::tan(std::acos(-1.0) * (random_numbers[0] - 0.5));

		if (channel == 0)
		{
			s += s0;
		}
		else if (channel == 1)
		{
			s += s1;
		}

		coordinates[0] = s;

		double const sms0 = s - s0;
		double const sms1 = s - s1;

		// set the channel densities
		densities[1] = 1.0 / (1.0 + sms1 * sms1) / std::acos(-1.0);
		densities[0] = 1.0 / (1.0 + sms0 * sms0) / std::acos(-1.0);

		// global weight
		return 1.0;
	};

#ifndef USE_MPI
	hep::multi_channel_callback<double>(
		hep::multi_channel_verbose_callback<double>);
#else
	hep::mpi_multi_channel_callback<double>(
		hep::mpi_multi_channel_verbose_callback<double>);
#endif

#ifndef USE_MPI
	auto const results = hep::multi_channel(
#else
	auto const results = hep::mpi_multi_channel(
		MPI_COMM_WORLD,
#endif
		hep::make_multi_channel_integrand<double>(function, 1, densities, 1, 2),
		std::vector<std::size_t>(10, 10000000)
	);

#ifdef USE_MPI
	int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0)
	{
#endif

	for (auto const result : results)
	{
		std::cout << result.value() << " +- " << result.error() << "\t"
			<< result.channel_weights()[0];

		if (result.channel_weights().size() == 2)
		{
			std::cout << " " << result.channel_weights()[1];
		}

		std::cout << "\n";
	}

#ifdef USE_MPI
	}

	MPI_Finalize();
#endif

	return 0;
}
