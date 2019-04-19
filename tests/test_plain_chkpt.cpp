#include "hep/mc/integrand.hpp"
#include "hep/mc/plain.hpp"

#include "catch2/catch.hpp"

#include <cstddef>
#include <sstream>
#include <vector>

template <typename T>
T linear_function(hep::mc_point<T> const& point, hep::projector<T>& projector)
{
    T const x = point.point().at(0);
    T const v = x;

    projector.add(0, x, v);

    return v;
}

TEMPLATE_TEST_CASE("plain_chkpt serialization", "", float, double, long double)
{
    using T = TestType;

    // create a checkpoint using five iterations with PLAIN
    auto const chkpt1 = hep::plain(
        hep::make_integrand<T>(
            linear_function<T>,
            1,
            hep::make_dist_params<T>(10, T(0.0), T(1.0), "distribution #1")
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
    CHECK( stream1.str() == stream2.str() );
}

TEMPLATE_TEST_CASE("plain_chkpt series", "", float, double, long double)
{
    using T = TestType;

    // check if we get the same checkpoint when we perform five iterations in one go, and when we
    // perfom first one iteration, then resume the resulting checkpoint with three iterations, and
    // finally resume it with one iteration

    auto const chkptA = hep::plain(
        hep::make_integrand<T>(
            linear_function<T>,
            1,
            hep::make_dist_params<T>(10, T(0.0), T(1.0), "distribution #1")
        ),
        std::vector<std::size_t>(5, 1000)
    );

    auto const tmp_chkpt1 = hep::plain(
        hep::make_integrand<T>(
            linear_function<T>,
            1,
            hep::make_dist_params<T>(10, T(0.0), T(1.0), "distribution #1")
        ),
        std::vector<std::size_t>(1, 1000)
    );

    auto const tmp_chkpt2 = hep::plain(
        hep::make_integrand<T>(
            linear_function<T>,
            1,
            hep::make_dist_params<T>(10, T(0.0), T(1.0), "distribution #1")
        ),
        std::vector<std::size_t>(3, 1000),
        tmp_chkpt1
    );

    auto const chkptB = hep::plain(
        hep::make_integrand<T>(
            linear_function<T>,
            1,
            hep::make_dist_params<T>(10, T(0.0), T(1.0), "distribution #1")
        ),
        std::vector<std::size_t>(1, 1000),
        tmp_chkpt2
    );

    // Check if both checkpoints are the same
    std::ostringstream stream1;
    std::ostringstream stream2;

    chkptA.serialize(stream1);
    chkptB.serialize(stream2);

    CHECK( stream1.str() == stream2.str() );
}
