#include "gtest/gtest.h"

#include <algorithm>
#include <cstddef>
//#include <iostream>
#include <limits>
#include <vector>

#ifndef HEP_USE_MPI
#include "hep/mc.hpp"
#else
#include "hep/mc-mpi.hpp"
#endif

template <typename T>
T function(hep::mc_point<T> const& x)
{
    return T(3.0) / T(2.0) * (x.point()[0] * x.point()[0] + x.point()[1] * x.point()[1]);
}

template <typename T>
T function_with_distribution(
    hep::mc_point<T> const& x,
    hep::projector<T>& projector
) {
    T const value = T(3.0) / T(2.0) * (x.point()[0] * x.point()[0] + x.point()[1] * x.point()[1]);

    projector.add(0, T(0.5), value);

    return value;
}

template <typename T>
std::vector<T> reference_results();

template <>
std::vector<float> reference_results() {
    return {
        // PLAIN
        9.980040789e-01f, 1.994426595e-03f,
        // VEGAS
        9.980040789e-01f, 1.994426595e-03f,
        1.001339316e+00f, 1.100822701e-03f,
        1.000011802e+00f, 7.872990682e-04f,
        9.995998144e-01f, 6.881667650e-04f,
        9.999216199e-01f, 6.571758422e-04f,
        // MULTI CHANNEL
        9.991801381e-01f, 1.994440332e-03f,
        1.000196576e+00f, 2.003321424e-03f,
        9.994604588e-01f, 1.999591244e-03f,
        9.999244809e-01f, 1.995231956e-03f,
        1.000989318e+00f, 1.998512540e-03f
    };
}

template <>
std::vector<double> reference_results() {
    return {
        // PLAIN
        9.99946966830306927e-01, 1.99451738543641663e-03,
        // VEGAS
        9.99946966830306927e-01, 1.99451738543641663e-03,
        1.00056702660240293e+00, 1.10050281971454359e-03,
        1.00009739823556631e+00, 7.84370395933593040e-04,
        9.99527435913331264e-01, 6.90829838859883150e-04,
        9.99963655759120829e-01, 6.62088512423425978e-04,
        // MULTI CHANNEL
        1.00189177191546164e+00, 2.00036054797089529e-03,
        1.00067313510539035e+00, 1.99717780871673479e-03,
        1.00199106786454206e+00, 2.01179040426654543e-03,
        1.00226985613340203e+00, 1.99751827870606581e-03,
        1.00237676222675454e+00, 2.00439702757806981e-03
    };
}

template <>
std::vector<long double> reference_results() {
    return {
        // PLAIN
        9.999469668303068833263e-01l, 1.994517385436427074006e-03l,
        // VEGAS
        9.999469668303068833263e-01l, 1.994517385436427074006e-03l,
        1.000567026602402857837e+00l, 1.100502819714536443759e-03l,
        1.000097398235566315537e+00l, 7.843703959335768267279e-04l,
        9.995274359133312137328e-01l, 6.908298388598802587929e-04l,
        9.999636557591208273232e-01l, 6.620885124234442225630e-04l,
        // MULTI CHANNEL
        1.001891771915461505746e+00l, 2.000360547970930202107e-03l,
        1.000673135105390259723e+00l, 1.997177808716723910491e-03l,
        1.001991067864541958866e+00l, 2.011790404266571298049e-03l,
        1.002269856133402095794e+00l, 1.997518278706081122539e-03l,
        1.002376762226754475585e+00l, 2.004397027578075676893e-03l
    };
}

typedef testing::Types<float, double, long double> MyT;
template <typename T> class NumericalResults : public testing::Test { };
TYPED_TEST_CASE(NumericalResults, MyT);

TYPED_TEST(NumericalResults, CheckPlainIntegration)
{
    using T = TypeParam;

    std::size_t const calls = 100000;

#ifndef HEP_USE_MPI
    auto const results = hep::plain(
#else
    auto const results = hep::mpi_plain(
        MPI_COMM_WORLD,
#endif
        hep::make_integrand<T>(function<T>, 2),
        std::vector<std::size_t>(1, calls)
    ).results();
    auto const reference = reference_results<T>();

#ifndef HEP_USE_MPI
    T const value_abs_error = T();
    T const error_abs_error = T();
#else
    T const value_abs_error = T(4.0) * std::numeric_limits<T>::epsilon();
    T const error_abs_error = T(4.0) * std::numeric_limits<T>::epsilon();
#endif

    EXPECT_NEAR( results.front().value() , reference[0] , value_abs_error );
    EXPECT_NEAR( results.front().error() , reference[1] , error_abs_error );

    EXPECT_EQ( results.front().calls() , calls );
    EXPECT_EQ( results.front().non_zero_calls() , calls );
    EXPECT_EQ( results.front().finite_calls() , calls );

//    std::cout.precision(std::numeric_limits<T>::max_digits10);
//    std::cout.setf(std::ios_base::scientific);
//    std::cout << result.value() << ", " << result.error() << ",\n";
}

TYPED_TEST(NumericalResults, CheckVegasIntegration)
{
    using T = TypeParam;

    std::size_t const calls = 100000;
    std::size_t const iterations = 5;

#ifndef HEP_USE_MPI
    auto const results = hep::vegas(
#else
    auto const results = hep::mpi_vegas(
        MPI_COMM_WORLD,
#endif
        hep::make_integrand<T>(function<T>, 2),
        std::vector<std::size_t>(iterations, calls)
    ).results();
    auto const reference = reference_results<T>();

#ifndef HEP_USE_MPI
    T const value_abs_error = T();
    T const error_abs_error = T();
#else
    T const value_abs_error = T(4.0) * std::numeric_limits<T>::epsilon();
    T const error_abs_error = T(4.0) * std::numeric_limits<T>::epsilon();
#endif

    ASSERT_EQ( results.size() , iterations );

    EXPECT_NEAR( results[0].value() , reference[ 2] , value_abs_error );
    EXPECT_NEAR( results[0].error() , reference[ 3] , error_abs_error );
    EXPECT_NEAR( results[1].value() , reference[ 4] , value_abs_error );
    EXPECT_NEAR( results[1].error() , reference[ 5] , error_abs_error );
    EXPECT_NEAR( results[2].value() , reference[ 6] , value_abs_error );
    EXPECT_NEAR( results[2].error() , reference[ 7] , error_abs_error );
    EXPECT_NEAR( results[3].value() , reference[ 8] , value_abs_error );
    EXPECT_NEAR( results[3].error() , reference[ 9] , error_abs_error );
    EXPECT_NEAR( results[4].value() , reference[10] , value_abs_error );
    EXPECT_NEAR( results[4].error() , reference[11] , error_abs_error );

    EXPECT_EQ( results[0].calls() , calls );
    EXPECT_EQ( results[1].calls() , calls );
    EXPECT_EQ( results[2].calls() , calls );
    EXPECT_EQ( results[3].calls() , calls );
    EXPECT_EQ( results[4].calls() , calls );

    EXPECT_EQ( results[0].non_zero_calls() , calls );
    EXPECT_EQ( results[1].non_zero_calls() , calls );
    EXPECT_EQ( results[2].non_zero_calls() , calls );
    EXPECT_EQ( results[3].non_zero_calls() , calls );
    EXPECT_EQ( results[4].non_zero_calls() , calls );

    EXPECT_EQ( results[0].finite_calls() , calls );
    EXPECT_EQ( results[1].finite_calls() , calls );
    EXPECT_EQ( results[2].finite_calls() , calls );
    EXPECT_EQ( results[3].finite_calls() , calls );
    EXPECT_EQ( results[4].finite_calls() , calls );

//    std::cout.precision(std::numeric_limits<T>::max_digits10);
//    std::cout.setf(std::ios_base::scientific);
//    for (auto i : results)
//    {
//        std::cout << i.value() << ", " << i.error() << ",\n";
//    }
}

TYPED_TEST(NumericalResults, CheckMultiChannelIntegration)
{
    using T = TypeParam;

    std::size_t const calls = 100000;
    std::size_t const iterations = 5;

    auto unit_densities = [](
        std::size_t,
        std::vector<T> const& random_numbers,
        std::vector<T>& coordinates,
        std::vector<std::size_t> const&,
        std::vector<T>& densities,
        hep::multi_channel_map action
    ) {
        if (action == hep::multi_channel_map::calculate_densities)
        {
            // we calculated the densities already
            return T(1.0);
        }

        std::copy(
            random_numbers.begin(),
            random_numbers.end(),
            coordinates.begin()
        );

        std::fill(densities.begin(), densities.end(), T(1.0));

        return T(1.0);
    };

#ifndef HEP_USE_MPI
    auto const results = hep::multi_channel(
#else
    auto const results = hep::mpi_multi_channel(
        MPI_COMM_WORLD,
#endif
        hep::make_multi_channel_integrand<T>(
            function<T>,
            2,
            unit_densities,
            2,
            2
        ),
        std::vector<std::size_t>(iterations, calls)
    ).results();
    auto const reference = reference_results<T>();

#ifndef HEP_USE_MPI
    T const value_abs_error = T();
    T const error_abs_error = T();
#else
    T const value_abs_error = T(4.0) * std::numeric_limits<T>::epsilon();
    T const error_abs_error = T(4.0) * std::numeric_limits<T>::epsilon();
#endif

    ASSERT_EQ( results.size() , iterations );

    EXPECT_NEAR( results[0].value() , reference[12] , value_abs_error );
    EXPECT_NEAR( results[0].error() , reference[13] , error_abs_error );
    EXPECT_NEAR( results[1].value() , reference[14] , value_abs_error );
    EXPECT_NEAR( results[1].error() , reference[15] , error_abs_error );
    EXPECT_NEAR( results[2].value() , reference[16] , value_abs_error );
    EXPECT_NEAR( results[2].error() , reference[17] , error_abs_error );
    EXPECT_NEAR( results[3].value() , reference[18] , value_abs_error );
    EXPECT_NEAR( results[3].error() , reference[19] , error_abs_error );
    EXPECT_NEAR( results[4].value() , reference[20] , value_abs_error );
    EXPECT_NEAR( results[4].error() , reference[21] , error_abs_error );

    EXPECT_EQ( results[0].calls() , calls );
    EXPECT_EQ( results[1].calls() , calls );
    EXPECT_EQ( results[2].calls() , calls );
    EXPECT_EQ( results[3].calls() , calls );
    EXPECT_EQ( results[4].calls() , calls );

    EXPECT_EQ( results[0].non_zero_calls() , calls );
    EXPECT_EQ( results[1].non_zero_calls() , calls );
    EXPECT_EQ( results[2].non_zero_calls() , calls );
    EXPECT_EQ( results[3].non_zero_calls() , calls );
    EXPECT_EQ( results[4].non_zero_calls() , calls );

    EXPECT_EQ( results[0].finite_calls() , calls );
    EXPECT_EQ( results[1].finite_calls() , calls );
    EXPECT_EQ( results[2].finite_calls() , calls );
    EXPECT_EQ( results[3].finite_calls() , calls );
    EXPECT_EQ( results[4].finite_calls() , calls );

//    std::cout.precision(std::numeric_limits<T>::max_digits10);
//    std::cout.setf(std::ios_base::scientific);
//    for (auto i : results)
//    {
//        std::cout << i.value() << ", " << i.error() << ",\n";
//    }

#ifndef HEP_USE_MPI
    auto const results2 = hep::multi_channel(
#else
    auto const results2 = hep::mpi_multi_channel(
        MPI_COMM_WORLD,
#endif
        hep::make_multi_channel_integrand<T>(
            function_with_distribution<T>,
            2,
            unit_densities,
            2,
            2,
            hep::make_dist_params<T>(1, T(), T(1.0))
        ),
        std::vector<std::size_t>(iterations, calls)
    ).results();

    ASSERT_EQ( results2.size() , iterations );
    ASSERT_EQ( results2[0].distributions().size() , 1u );
    ASSERT_EQ( results2[1].distributions().size() , 1u );
    ASSERT_EQ( results2[2].distributions().size() , 1u );
    ASSERT_EQ( results2[3].distributions().size() , 1u );
    ASSERT_EQ( results2[4].distributions().size() , 1u );
    ASSERT_EQ( results2[0].distributions()[0].results().size() , 1u );
    ASSERT_EQ( results2[1].distributions()[0].results().size() , 1u );
    ASSERT_EQ( results2[2].distributions()[0].results().size() , 1u );
    ASSERT_EQ( results2[3].distributions()[0].results().size() , 1u );
    ASSERT_EQ( results2[4].distributions()[0].results().size() , 1u );

    EXPECT_EQ( results2[0].value() , results[0].value() );
    EXPECT_EQ( results2[0].error() , results[0].error() );
    EXPECT_EQ( results2[0].calls() , results[0].calls() );
    EXPECT_EQ( results2[1].value() , results[1].value() );
    EXPECT_EQ( results2[1].error() , results[1].error() );
    EXPECT_EQ( results2[1].calls() , results[1].calls() );
    EXPECT_EQ( results2[2].value() , results[2].value() );
    EXPECT_EQ( results2[2].error() , results[2].error() );
    EXPECT_EQ( results2[2].calls() , results[2].calls() );
    EXPECT_EQ( results2[3].value() , results[3].value() );
    EXPECT_EQ( results2[3].error() , results[3].error() );
    EXPECT_EQ( results2[3].calls() , results[3].calls() );
    EXPECT_EQ( results2[4].value() , results[4].value() );
    EXPECT_EQ( results2[4].error() , results[4].error() );
    EXPECT_EQ( results2[4].calls() , results[4].calls() );

    EXPECT_EQ( results2[0].non_zero_calls() , results[0].calls() );
    EXPECT_EQ( results2[1].non_zero_calls() , results[1].calls() );
    EXPECT_EQ( results2[2].non_zero_calls() , results[2].calls() );
    EXPECT_EQ( results2[3].non_zero_calls() , results[3].calls() );
    EXPECT_EQ( results2[4].non_zero_calls() , results[4].calls() );

    EXPECT_EQ( results2[0].finite_calls() , results[0].calls() );
    EXPECT_EQ( results2[1].finite_calls() , results[1].calls() );
    EXPECT_EQ( results2[2].finite_calls() , results[2].calls() );
    EXPECT_EQ( results2[3].finite_calls() , results[3].calls() );
    EXPECT_EQ( results2[4].finite_calls() , results[4].calls() );

    EXPECT_EQ( mid_points_x(results2[0].distributions()[0])[0] , T(0.5) );
    EXPECT_EQ( results2[0].distributions()[0].results()[0].value() , results[0].value() );
    EXPECT_EQ( results2[0].distributions()[0].results()[0].error() , results[0].error() );
    EXPECT_EQ( results2[0].distributions()[0].results()[0].calls() , results[0].calls() );
    EXPECT_EQ( results2[0].distributions()[0].results()[0].non_zero_calls() , results[0].calls() );
    EXPECT_EQ( results2[0].distributions()[0].results()[0].finite_calls() , results[0].calls() );

    EXPECT_EQ( mid_points_x(results2[1].distributions()[0])[0] , T(0.5) );
    EXPECT_EQ( results2[1].distributions()[0].results()[0].value() , results[1].value() );
    EXPECT_EQ( results2[1].distributions()[0].results()[0].error() , results[1].error() );
    EXPECT_EQ( results2[1].distributions()[0].results()[0].calls() , results[1].calls() );
    EXPECT_EQ( results2[1].distributions()[0].results()[0].non_zero_calls() , results[1].calls() );
    EXPECT_EQ( results2[1].distributions()[0].results()[0].finite_calls() , results[1].calls() );

    EXPECT_EQ( mid_points_x(results2[2].distributions()[0])[0] , T(0.5) );
    EXPECT_EQ( results2[2].distributions()[0].results()[0].value() , results[2].value() );
    EXPECT_EQ( results2[2].distributions()[0].results()[0].error() , results[2].error() );
    EXPECT_EQ( results2[2].distributions()[0].results()[0].calls() , results[2].calls() );
    EXPECT_EQ( results2[2].distributions()[0].results()[0].non_zero_calls() , results[2].calls() );
    EXPECT_EQ( results2[2].distributions()[0].results()[0].finite_calls() , results[2].calls() );

    EXPECT_EQ( mid_points_x(results2[3].distributions()[0])[0] , T(0.5) );
    EXPECT_EQ( results2[3].distributions()[0].results()[0].value() , results[3].value() );
    EXPECT_EQ( results2[3].distributions()[0].results()[0].error() , results[3].error() );
    EXPECT_EQ( results2[3].distributions()[0].results()[0].calls() , results[3].calls() );
    EXPECT_EQ( results2[3].distributions()[0].results()[0].non_zero_calls() , results[3].calls() );
    EXPECT_EQ( results2[3].distributions()[0].results()[0].finite_calls() , results[3].calls() );

    EXPECT_EQ( mid_points_x(results2[4].distributions()[0])[0] , T(0.5) );
    EXPECT_EQ( results2[4].distributions()[0].results()[0].value() , results[4].value() );
    EXPECT_EQ( results2[4].distributions()[0].results()[0].error() , results[4].error() );
    EXPECT_EQ( results2[4].distributions()[0].results()[0].calls() , results[4].calls() );
    EXPECT_EQ( results2[4].distributions()[0].results()[0].non_zero_calls() , results[4].calls() );
    EXPECT_EQ( results2[4].distributions()[0].results()[0].finite_calls() , results[4].calls() );
}
