#include <hep/mc/vegas.hpp>

#include <boost/mpl/list.hpp>
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <cstddef>
#include <vector>

template <typename T>
class integrand
{
public:
	integrand(std::size_t type)
		: type(type)
	{
	}

	T operator()(hep::vegas_sample<T> const& sample) const
	{
		std::vector<T> const& x = sample.point;

		switch (type)
		{
		case 0:
			return T(1.0);

		case 1:
			return std::acos(T(-1.0)) / T(2.0) *
				std::sin(std::acos(T(-1.0)) * x[0]);

		case 2:
			return std::exp(x[0]);

		case 3:
			return x[0] * x[0];

		case 4:
			return x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3];

		case 5: // not square-integrable
			return T(1.0) / std::sqrt(x[0]);

		case 6:
			return T(3.0) / T(2.0) * (T(1.0) -
				(T(2.0) * x[0] - T(1.0)) * (T(2.0) * x[0] - T(1.0)));

		case 7:
			return T(2.0) * x[0];

		default:
			return T(0.0);
		}
	}

private:
	std::size_t type;
};

typedef boost::mpl::list<float, double, long double> test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE(vegas_integration_square, T, test_types)
{
	std::size_t const iterations = 3;
	std::size_t const steps_per_iteration = 1000000;
	std::size_t const bins = 100;

	auto result = hep::vegas<T>(
		1,
		std::vector<std::size_t>(iterations, steps_per_iteration),
		steps_per_iteration,
		bins,
		integrand<T>(6)
	);

	BOOST_CHECK_CLOSE(
		result.values.back(),
		T(1.0),
		T(300.0) * result.errors.back()
	);
}
