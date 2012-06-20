#ifndef HEP_MC_VEGAS_RESULT_HPP
#define HEP_MC_VEGAS_RESULT_HPP

#include <cstddef>
#include <vector>

namespace hep
{

/**
 * The result of a VEGAS Monte Carlo integration.
 */
template <typename T>
class vegas_result
{
public:
	/**
	 * 
	 */
	vegas_result(std::vector<std::size_t> const& steps);

	/**
	 * 
	 */
	void add_iteration(
		std::size_t samples,
		T const& sum,
		T const& sum_of_squares
	);

	/**
	 * The computed approximations for the integral for each iteration.
	 */
	std::vector<T> values;

	/**
	 * The errors for each iteration.
	 */
	std::vector<T> errors;

	/**
	 * A std::vector containing the number of steps performed by VEGAS for each
	 * iteration.
	 */
	std::vector<std::size_t> steps;

private:
	/**
	 * 
	 */
	T sum_of_inv_variances;

	/**
	 * 
	 */
	T sum_of_averages;
};

}

#include <hep/mc/impl/vegas_result.hpp>

#endif
