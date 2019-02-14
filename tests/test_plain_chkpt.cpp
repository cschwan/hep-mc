#include "gtest/gtest.h"

#include "hep/mc.hpp"

#include <cstddef>
#include <fstream>
#include <vector>

using T = double;

T linear_function(hep::mc_point<T> const& point, hep::projector<T>& projector)
{
    T const x = point.point().at(0);
    T const v = x;

    projector.add(0, x, v);

    return v;
}

TEST(TestPlainCheckpoint, Bla)
{
    // create a checkpoint using five iterations with PLAIN
    auto const chkpt1 = hep::plain(
        hep::make_integrand<T>(
            linear_function,
            1,
            hep::make_dist_params<T>(10, T(0.0), T(1.0), "distribution name")
        ),
        std::vector<std::size_t>(5, 1000)
    );

    // serialize checkpoint
    std::ostringstream stream1;
    chkpt1.serialize(stream1);

    // unserialize checkpoint
    auto in = std::istringstream(stream1.str());
    auto const chkpt2 = decltype (chkpt1) {in};

    // serialize previous checkpoint
    std::ostringstream stream2;
    chkpt2.serialize(stream2);

    // string representation should be the same
    ASSERT_STREQ( stream1.str().c_str() , stream2.str().c_str() );
}
