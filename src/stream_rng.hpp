#ifndef STREAM_RNG_HPP
#define STREAM_RNG_HPP

#include <iosfwd>
#include <string>

namespace hep
{

class stream_rng
{
public:
    stream_rng();

    std::string const& state() const;

    void state(std::string const& state);

private:
    std::string state_;
};

std::istream& operator>>(std::istream& in, stream_rng& rng);

std::ostream& operator<<(std::ostream& in, stream_rng const& rng);

}

#endif
