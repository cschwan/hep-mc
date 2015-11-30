#include "hep/mc.hpp"

#include <iostream>
#include <vector>

int main()
{
	// vector containing the results of five (fictious) iterations
	std::vector<hep::mc_result<double>> results = {
		// create mc_results with the tuples (calls, estimate, error)
		hep::create_result(100000, 0.987449, 0.0165879),
		hep::create_result(100000, 0.988517, 0.00929261),
		hep::create_result(100000, 0.999915, 0.00825228),
		hep::create_result(100000, 1.0012,   0.0120161),
		hep::create_result(100000, 1.01968,  0.00926521)
	};

	for (std::size_t i = 0; i != results.size(); ++i)
	{
		// combine results from iterations [0, i]
		auto result = hep::cumulative_result0(results.begin(),
			results.begin() + (i+1));
		auto chi_square_dof = hep::chi_square_dof0(results.begin(),
			results.begin() + (i+1));

		// print results
		std::cout << ">> printing cumulative result for results 0 till " << i
			<< "\nmethod 0: E=" << result.value() << " +- " << result.error()
			<< " chi^2/dof=" << chi_square_dof << "\n";

		// calculate cumulative results again, this time with a different method
		result = hep::cumulative_result1(results.begin(),
			results.begin() + (i+1));
		chi_square_dof = hep::chi_square_dof1(results.begin(),
			results.begin() + (i+1));

		// print results
		std::cout << "method 1: E=" << result.value() << " +- "
			<< result.error() << " chi^2/dof=" << chi_square_dof << "\n\n";
	}

	return 0;
}
