#ifndef GENZ_INTEGRAND_HPP
#define GENZ_INTEGRAND_HPP

#include "hep/mc/mc_point.hpp"

#include <cstddef>
#include <vector>

namespace genz
{

enum integrand_type
{
	oscillatory,
	product_peak,
	corner_peak,
	gaussian,
	c0_function,
	discontinuous
};

template <typename T>
class integrand
{
public:
	integrand(
		integrand_type type,
		std::vector<T> const& a,
		std::vector<T> const& u
	);

	T operator()(hep::mc_point<T> const& x) const;

	T reference_result() const;

private:
	integrand_type type;
	std::vector<T> a;
	std::vector<T> u;
};

template <typename T>
class parameters
{
public:
	parameters(
		std::size_t dimension,
		T difficulty,
		T compatibility,
		T limit = T(0.05)
	);

	std::vector<T> affective() const;
	std::vector<T> unaffective() const;

private:
	std::vector<T> a;
	std::vector<T> u;
};

}

#endif
