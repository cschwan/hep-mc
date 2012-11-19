#include <boost/test/unit_test.hpp>
#include <hep/mc/vegas.hpp>

#include <cmath>
#include <cstddef>
#include <vector>

template <typename T>
T integrand(hep::vegas_sample<T> const& sample, std::size_t type)
{
	std::vector<T> const& x = sample.point;

	const T pi = std::atan(T(1.0)) * T(4.0);

	switch (type)
	{
	case 0:
		return T(1.0);

	case 1:
		return pi / T(2.0) * std::sin(pi * x[0]);

	case 2:
		return std::exp(x[0]);

	case 3:
		return x[0]*x[0];

	case 4:
		return x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3];

	case 5: // not square-integrable
		return T(1.0) / std::sqrt(x[0]);

	case 6:
		return T(3.0 / 2.0) * (T(1.0) -
			(T(2.0) * x[0] - T(1.0)) * (T(2.0) * x[0] - T(1.0)));

	case 7:
		return 2.0 * x[0];

	default:
		return T(0.0);
	}
}

BOOST_AUTO_TEST_CASE(vegas_integration_square)
{
	std::size_t const iterations = 3;
	std::size_t const steps_per_iteration = 1000000;
	std::size_t const bins = 100;

	typedef double float_type;

	auto result = hep::vegas<float_type>(
		1,
		std::vector<std::size_t>(iterations, steps_per_iteration),
		steps_per_iteration,
		bins,
		integrand<float_type>,
		6
	);

	BOOST_CHECK_CLOSE(
		result.values.back(),
		1.0,
		300.0 * result.errors.back()
	);
}
