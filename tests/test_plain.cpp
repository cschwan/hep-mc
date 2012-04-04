#include <boost/test/unit_test.hpp>
#include <hep/mc/plain.hpp>

#include <cmath>
#include <cstddef>
#include <vector>

double integrand(std::vector<double> const& x, std::size_t type)
{
	switch (type)
	{
	case 0:
		return 1.0;

	case 1:
		return M_PI/2.0 * std::sin(M_PI*x[0]);

	case 2:
		return std::exp(x[0]);

	case 3:
		return x[0]*x[0];

	case 4:
		return x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3];

	case 5: // not square-integrable
		return 1.0 / std::sqrt(x[0]);

	default:
		return 0.0;
	}
}

BOOST_AUTO_TEST_CASE(plain_integration_constant)
{
	std::size_t steps = 10;
	hep::plain_result<double> result =
		hep::plain<double>(1, steps, 0, integrand, 0);

	BOOST_CHECK_EQUAL(result.value, 1.0);
	BOOST_CHECK_EQUAL(result.error, 0.0);
}

BOOST_AUTO_TEST_CASE(plain_integration_sine)
{
	std::size_t steps = 1000;
	hep::plain_result<double> result =
		hep::plain<double>(1, steps, 0, integrand, 1);

	BOOST_CHECK_CLOSE(result.value, 1.0, 100.0 / std::sqrt(steps));
	BOOST_CHECK_CLOSE(result.value, 1.0, 100.0 * result.error);
}

BOOST_AUTO_TEST_CASE(plain_integration_exponential)
{
	std::size_t steps = 1000;
	hep::plain_result<double> result =
		hep::plain<double>(1, steps, 0, integrand, 2);
	double reference = std::exp(1.0) - 1.0;

	BOOST_CHECK_CLOSE(result.value, reference, 100.0 / std::sqrt(steps));
	BOOST_CHECK_CLOSE(result.value, reference, 100.0 * result.error);
}

BOOST_AUTO_TEST_CASE(plain_integration_square)
{
	std::size_t steps = 10000;
	hep::plain_result<double> result =
		hep::plain<double>(1, steps, 0, integrand, 3);

	BOOST_CHECK_CLOSE(result.value, 1.0 / 3.0, 100.0 / std::sqrt(steps));
	BOOST_CHECK_CLOSE(result.value, 1.0 / 3.0, 100.0 * result.error);
}

BOOST_AUTO_TEST_CASE(plain_integration_square_4d)
{
	std::size_t steps = 10000;
	hep::plain_result<double> result =
		hep::plain<double>(4, steps, 0, integrand, 4);

	BOOST_CHECK_CLOSE(result.value, 4.0 / 3.0, 100.0 / std::sqrt(steps));
	BOOST_CHECK_CLOSE(result.value, 4.0 / 3.0, 100.0 * result.error);
}

BOOST_AUTO_TEST_CASE(plain_integration_one_over_sqrt)
{
	std::size_t steps = 1000;
	hep::plain_result<double> result =
		hep::plain<double>(1, steps, 0, integrand, 5);

	BOOST_CHECK_CLOSE(result.value, 2.0, 100.0 / std::sqrt(steps));
	// function is not square-integrable, error is unreliable
}
