#ifndef HEP_PS_DEFAULT_VEGAS_PARALLELIZER_HPP
#define HEP_PS_DEFAULT_VEGAS_PARALLELIZER_HPP

#include <vector>

namespace hep
{

/**
 * Default VEGAS parallelizer that does not parallelize.
 */
template <typename T>
class default_vegas_parallelizer
{
public:
	/**
	 * Returns zero.
	 */
	std::size_t rank() const
	{
		return 0;
	}

	/**
	 * Returns one.
	 */
	std::size_t world_size() const
	{
		return 1;
	}

	/**
	 * Does nothing.
	 */
	void all_reduce(std::vector<T> const& values) const
	{
	}
};

}

#endif
