#include "genz_integrand.hpp"

#include <algorithm>
#include <bitset>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <random>

namespace genz
{

template <typename T>
integrand<T>::integrand(
	integrand_type type,
	std::vector<T> const& a,
	std::vector<T> const& u
)
	: type(type)
	, a(a)
	, u(u)
{
	assert( a.size() == u.size() );
}

template <typename T>
T integrand<T>::operator()(hep::mc_point<T> const& x) const
{
	std::size_t const s = a.size();

	assert( x.point().size() == s );

	T ax = T();

	for (std::size_t i = 0; i != s; ++i)
	{
		ax += a[i] * x.point()[i];
	}

	T exponent = T();
	T result = T(1.0);

	switch (type)
	{
	case oscillatory:
		return std::cos(T(2.0) * std::acos(T(-1.0)) * u[0] + ax);

	case product_peak:
		for (std::size_t i = 0; i != s; ++i)
		{
			T const xmu = x.point()[i] - u[i];
			result /= T(1.0) / (a[i] * a[i]) + xmu * xmu;
		}

		return result;

	case corner_peak:
		// do not put sign of exponent into T( ... ) because s is unsigned!
		return std::pow(T(1.0) + ax, -T(s+1));

	case gaussian:
		for (std::size_t i = 0; i != s; ++i)
		{
			T const xmu = x.point()[i] - u[i];
			exponent += a[i] * a[i] * xmu * xmu;
		}

		return std::exp(-exponent);

	case c0_function:
		for (std::size_t i = 0; i != s; ++i)
		{
			exponent += a[i] * std::fabs(x.point()[i] - u[i]);
		}

		return std::exp(-exponent);

	case discontinuous:
		if ((x.point()[0] > u[0]) || ((x.point().size() > 1) &&
			(x.point()[1] > u[1])))
		{
			return T();
		}

		return std::exp(ax);

	default:
		assert( false );
	}
}

template <typename T>
T integrand<T>::reference_result() const
{
	std::size_t const s = a.size();

	T result = T(1.0);
	T product = T(1.0);
	T sum = T();

	std::bitset<std::numeric_limits<std::uintmax_t>::digits> w;

	switch (type)
	{
	case oscillatory:
		for (std::size_t i = 0; i != s; ++i)
		{
			T const tmp = T(0.5) * a[i];
			sum += tmp;
			product *= std::sin(tmp) / tmp;
		}

		return std::cos(T(2.0) * std::acos(T(-1.0)) * u[0] + sum) * product;

	case product_peak:
		for (std::size_t i = 0; i != s; ++i)
		{
			result *= a[i] * (std::atan(a[i] * (T(1.0) - u[i]))
				- std::atan(-a[i] * u[i]));
		}

		return result;

	case corner_peak:
		// make sure we can represent 2^s in the bitset
		assert( s < w.size() );

		result = T();

		for (std::size_t i = 0; i != s; ++i)
		{
			product *= a[i];
		}

		for (std::uintmax_t i = 0; i != (1u << s); ++i)
		{
			w = i;

			sum = T(1.0);
			for (std::size_t j = 0; j != s; ++j)
			{
				sum += w[j] ? a[j] : T();
			}

			T const sign = (w.count() % 2 == 0) ? T(1.0) : T(-1.0);
			result += sign / sum;
		}

		return result / std::tgamma(T(s+1)) / product;

	case gaussian:
		for (std::size_t i = 0; i != s; ++i)
		{
			result *= std::sqrt(std::acos(T(-1.0))) * T(0.5) / a[i] *
				(std::erf(a[i] * (T(1.0) - u[i])) - std::erf(a[i] * -u[i]));
		}

		return result;

	case c0_function:
		for (std::size_t i = 0; i != s; ++i)
		{
			result *= (T(2.0) - std::exp(-a[i] * u[i])
				- std::exp(a[i] * (u[i] - T(1.0)))) / a[i];
		}

		return result;

	case discontinuous:
		for (std::size_t i = 0; i != s; ++i)
		{
			T const tprime = i < 2 ? std::fmax(T(), std::fmin(T(1.0), u[i]))
				: T(1.0);
			result *= (std::exp(a[i] * tprime) - T(1.0)) / a[i];
		}

		return result;

	default:
		assert( false );
	}
}

//template class integrand<float>;
template class integrand<double>;
//template class integrand<long double>;

template <typename T>
parameters<T>::parameters(
	std::size_t dimension,
	T difficulty,
	T compatibility,
	T limit
)
	: a(dimension)
	, u(dimension)
{
	std::mt19937 rng;
	std::uniform_real_distribution<T> distribution(limit, T(1.0) - limit);

	std::generate(a.begin(), a.end(), [&]() {
		return distribution(rng);
	});

	T length = std::accumulate(a.begin(), a.end(), T());
	T factor = difficulty / std::pow(T(dimension), compatibility) / length;

	std::transform(a.begin(), a.end(), a.begin(), [&](T value) {
		return factor * value;
	});

	std::generate(u.begin(), u.end(), [&]() {
		return distribution(rng);
	});
}

template <typename T>
std::vector<T> parameters<T>::affective() const
{
	return a;
}

template <typename T>
std::vector<T> parameters<T>::unaffective() const
{
	return u;
}

//template class parameters<float>;
template class parameters<double>;
//template class parameters<long double>;

}
