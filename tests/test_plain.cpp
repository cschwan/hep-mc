#include <hep/mc/plain.hpp>

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

	T operator()(hep::mc_point<T> const& point) const
	{
		std::vector<T> const& x = point.point;

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

		default:
			return T(0.0);
		}
	}

private:
	std::size_t type;
};

typedef boost::mpl::list<float, double, long double> test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE(constant, T, test_types)
{
	std::size_t steps = 10;
	hep::mc_result<T> result = hep::plain<T>(1, steps, integrand<T>(0));

	BOOST_CHECK_EQUAL(result.value, T(1.0));
	BOOST_CHECK_SMALL(result.error, T(1e-8));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(sine, T, test_types)
{
	std::size_t steps = 1000;
	hep::mc_result<T> result = hep::plain<T>(1, steps, integrand<T>(1));

	BOOST_CHECK_CLOSE(result.value, T(1.0), T(100.0) * result.error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(exponential, T, test_types)
{
	std::size_t steps = 1000;
	hep::mc_result<double> result =
		hep::plain<double>(1, steps, integrand<double>(2));
	double reference = std::exp(1.0) - 1.0;

	BOOST_CHECK_CLOSE(result.value, reference, T(100.0) * result.error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(square, T, test_types)
{
	std::size_t steps = 10000;
	hep::mc_result<double> result =
		hep::plain<double>(1, steps, integrand<double>(3));

	BOOST_CHECK_CLOSE(result.value, T(1.0) / T(3.0), T(300.0) * result.error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(square_4d, T, test_types)
{
	std::size_t steps = 10000;
	hep::mc_result<double> result =
		hep::plain<double>(4, steps, integrand<double>(4));

	BOOST_CHECK_CLOSE(result.value, T(4.0) / T(3.0), T(200.0) * result.error);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(one_over_sqrt, T, test_types)
{
	std::size_t steps = 1000;
	hep::mc_result<double> result =
		hep::plain<double>(1, steps, integrand<double>(5));

	BOOST_CHECK_CLOSE(result.value, 2.0, T(100.0) / std::sqrt(T(steps)));
	// function is not square-integrable, error is unreliable
}
