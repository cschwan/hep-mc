#include "hep/mc/vegas_pdf.hpp"

#include <catch2/catch.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

TEMPLATE_TEST_CASE("vegas_pdf construction", "", float, double, long double)
{
    using T = TestType;

    hep::vegas_pdf<T> pdf{2, 4};

    REQUIRE( pdf.bins()       == 4u );
    REQUIRE( pdf.dimensions() == 2u );

    CHECK( pdf.bin_left(0, 0) == T(0.00) );
    CHECK( pdf.bin_left(0, 1) == T(0.25) );
    CHECK( pdf.bin_left(0, 2) == T(0.50) );
    CHECK( pdf.bin_left(0, 3) == T(0.75) );
    CHECK( pdf.bin_left(0, 4) == T(1.00) );
    CHECK( pdf.bin_left(1, 0) == T(0.00) );
    CHECK( pdf.bin_left(1, 1) == T(0.25) );
    CHECK( pdf.bin_left(1, 2) == T(0.50) );
    CHECK( pdf.bin_left(1, 3) == T(0.75) );
    CHECK( pdf.bin_left(1, 4) == T(1.00) );
}

TEMPLATE_TEST_CASE("vegas_pdf setter", "", float, double, long double)
{
    using T = TestType;

    hep::vegas_pdf<T> pdf(2, 3);

    pdf.set_bin_left(0, 0, T(0.0));
    pdf.set_bin_left(0, 1, T(0.1));
    pdf.set_bin_left(0, 2, T(0.2));
    pdf.set_bin_left(0, 3, T(1.0));
    pdf.set_bin_left(1, 0, T(0.0));
    pdf.set_bin_left(1, 1, T(0.4));
    pdf.set_bin_left(1, 2, T(0.5));
    pdf.set_bin_left(1, 3, T(1.0));

    CHECK( pdf.bin_left(0, 0) == T(0.0) );
    CHECK( pdf.bin_left(0, 1) == T(0.1) );
    CHECK( pdf.bin_left(0, 2) == T(0.2) );
    CHECK( pdf.bin_left(0, 3) == T(1.0) );
    CHECK( pdf.bin_left(1, 0) == T(0.0) );
    CHECK( pdf.bin_left(1, 1) == T(0.4) );
    CHECK( pdf.bin_left(1, 2) == T(0.5) );
    CHECK( pdf.bin_left(1, 3) == T(1.0) );
}

TEMPLATE_TEST_CASE("vegas_icdf uniform", "", float, double /*, long double*/)
{
    using T = TestType;

    std::vector<T> random_numbers(1);
    std::vector<std::size_t> bin(1);

    // check a uniform pdf
    hep::vegas_pdf<T> pdf{1, 5};

    for (std::size_t i = 0; i != 99; ++i)
    {
        INFO( "i=" << i );

        T const number = T(1+i) / T(100.0);
        std::size_t const bin_index = number * pdf.bins();
        random_numbers[0] = number;

        // weight for a uniform pdf is 1
        CHECK_THAT( vegas_icdf(pdf, random_numbers, bin), Catch::WithinULP(T(1.0), 4)  );

        // the bins are equal-sized
        CHECK( bin[0] == bin_index );

        // "random" number must be unmodified
        CHECK_THAT( random_numbers[0] , Catch::WithinULP(number, 4) );
    }
}

TEMPLATE_TEST_CASE("vegas_icdf non uniform", "", float, double /*, long double*/)
{
    using T = TestType;

    std::vector<T> random_numbers(1);
    std::vector<std::size_t> bin(1);

    // check a uniform pdf
    hep::vegas_pdf<T> pdf{1, 5};

    // check a pdf in which bin 0 is twice as large as bin 1 ...
    pdf.set_bin_left(0, 0, T( 0.0)          );
    pdf.set_bin_left(0, 1, T( 2.0) / T( 3.0));
    pdf.set_bin_left(0, 2, T( 8.0) / T( 9.0));
    pdf.set_bin_left(0, 3, T(26.0) / T(27.0));
    pdf.set_bin_left(0, 4, T(80.0) / T(81.0));
    pdf.set_bin_left(0, 5, T( 1.0)          );

    std::array<T, 5> const weight = {
        T(10.0) / T(  3.0),
        T(10.0) / T(  9.0),
        T(10.0) / T( 27.0),
        T(10.0) / T( 81.0),
        T(15.0) / T(243.0)
    };

    std::array<T, 5> const offset = {
        T(),
        T( 2.0) / T(5.0),
        T(10.0) / T(5.0),
        T(36.0) / T(5.0),
        T(76.0) / T(5.0)
    };

    for (std::size_t i = 0; i != 99; ++i)
    {
        INFO( "i=" << i );

        T const number = T(1+i) / T(100.0);
        std::size_t const bin_index = number * pdf.bins();
        random_numbers[0] = number;

        CHECK_THAT( hep::vegas_icdf(pdf, random_numbers, bin) ,
            Catch::WithinULP(weight.at(bin_index), 64) );
        CHECK( bin[0]  == bin_index );
        CHECK_THAT( random_numbers[0] ,
            Catch::WithinULP((number + offset.at(bin_index)) * weight.at(bin_index), 4) );
    }
}

TEMPLATE_TEST_CASE("", "", float, double /*, long double*/)
{
    using std::exp;
    using T = TestType;

    // uniform distribution
    hep::vegas_pdf<T> old_pdf{1, 128};

    std::vector<T> adjustment_data(128 + 2);
    for (std::size_t i = 0; i != 128; ++i)
    {
        // gaussian data, but could also be noise
        T const power = (T(0.5) - i / T(128.0)) / T(0.1);
        adjustment_data.at(i) = T(1000.0) * exp(-power * power);
    }

    // refine with `alpha = 0` so that the PDF doesn't change
    auto new_pdf = hep::vegas_refine_pdf(old_pdf, T(), adjustment_data);

    for (std::size_t i = 0; i != old_pdf.bins() + 1; ++i)
    {
        CHECK_THAT( new_pdf.bin_left(0, i), Catch::WithinULP(old_pdf.bin_left(0, i), 4) );
    }
}
