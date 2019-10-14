#include "stream_rng.hpp"

#include <istream>

namespace hep
{

stream_rng::stream_rng() = default;

std::string const& stream_rng::state() const
{
    return state_;
}

void stream_rng::state(std::string const& state)
{
    state_ = state;
}

std::istream& operator>>(std::istream& in, stream_rng& rng)
{
    std::string state;
    std::getline(in, state);
    rng.state(state);

    return in;
}

std::ostream& operator<<(std::ostream& out, stream_rng const& rng)
{
    out << rng.state();

    return out;
}

}
