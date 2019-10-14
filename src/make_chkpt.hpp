#ifndef MAKE_CHKPT_HPP
#define MAKE_CHKPT_HPP

#include "hep/mc/chkpt.hpp"
#include "hep/mc/multi_channel_chkpt.hpp"
#include "hep/mc/plain_chkpt.hpp"
#include "hep/mc/vegas_chkpt.hpp"

namespace hep
{

template <typename Chkpt>
hep::chkpt_with_rng<hep::stream_rng, Chkpt> make_chkpt(std::istream& in);

template <>
hep::plain_chkpt_with_rng<hep::stream_rng, float> make_chkpt<hep::plain_chkpt<float>>(std::istream& in)
{
    return hep::make_plain_chkpt<float, hep::stream_rng>(in);
}

template <>
hep::plain_chkpt_with_rng<hep::stream_rng, double> make_chkpt<hep::plain_chkpt<double>>(std::istream& in)
{
    return hep::make_plain_chkpt<double, hep::stream_rng>(in);
}

template <>
hep::plain_chkpt_with_rng<hep::stream_rng, long double> make_chkpt<hep::plain_chkpt<long double>>(std::istream& in)
{
    return hep::make_plain_chkpt<long double, hep::stream_rng>(in);
}

template <>
hep::vegas_chkpt_with_rng<hep::stream_rng, float> make_chkpt<hep::vegas_chkpt<float>>(std::istream& in)
{
    return hep::make_vegas_chkpt<float, hep::stream_rng>(in);
}

template <>
hep::vegas_chkpt_with_rng<hep::stream_rng, double> make_chkpt<hep::vegas_chkpt<double>>(std::istream& in)
{
    return hep::make_vegas_chkpt<double, hep::stream_rng>(in);
}

template <>
hep::vegas_chkpt_with_rng<hep::stream_rng, long double> make_chkpt<hep::vegas_chkpt<long double>>(std::istream& in)
{
    return hep::make_vegas_chkpt<long double, hep::stream_rng>(in);
}

template <>
hep::multi_channel_chkpt_with_rng<hep::stream_rng, float> make_chkpt<hep::multi_channel_chkpt<float>>(std::istream& in)
{
    return hep::make_multi_channel_chkpt<float, hep::stream_rng>(in);
}

template <>
hep::multi_channel_chkpt_with_rng<hep::stream_rng, double> make_chkpt<hep::multi_channel_chkpt<double>>(std::istream& in)
{
    return hep::make_multi_channel_chkpt<double, hep::stream_rng>(in);
}

template <>
hep::multi_channel_chkpt_with_rng<hep::stream_rng, long double> make_chkpt<hep::multi_channel_chkpt<long double>>(std::istream& in)
{
    return hep::make_multi_channel_chkpt<long double, hep::stream_rng>(in);
}

}

#endif
