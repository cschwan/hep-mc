#include "hep/mc.hpp"

#include <cstddef>
#include <fstream>
#include <vector>

template <typename T>
T function(hep::mc_point<T> const& point, hep::projector<T>& projector)
{
    T const x = point.point().at(0);
    T const v = x;

    projector.add(0, x, v);

    return v;
}

void create_chkpt(std::string const& file)
{
    // create a checkpoint using five iterations with PLAIN and create a distribution
    auto const chkpt = hep::plain(
        hep::make_integrand<double>(
            function<double>,
            1,
            hep::make_dist_params<double>(10, 0.0, 1.0, "distribution #1")
        ),
        std::vector<std::size_t>(5, 1000),
        hep::make_plain_chkpt<double>(std::mt19937()) // default argument
    );

    // write the checkpoint to a file
    std::ofstream out(file);
    chkpt.serialize(out);
}

int main()
{
    // create a checkpoint and save it to disk
    create_chkpt("chkpt");

    // IMPORTANT: here we need to specify the random number generator type we've used before
    std::ifstream in("chkpt");
    auto const chkpt_from_disk = hep::make_plain_chkpt<double, std::mt19937>(in);

    // resume the integration with 5 additional iterations
    auto const chkpt = hep::plain(
        hep::make_integrand<double>(
            function<double>,
            1,
            hep::make_dist_params<double>(10, 0.0, 1.0, "distribution #1")
        ),
        std::vector<std::size_t>(5, 1000),
        chkpt_from_disk
    );

    std::cout << "result using five iterations from a chkpt and five addtional iterations:\n    ";

    // print the result using enough precision
    std::cout << std::scientific << std::setprecision(16)
        << chkpt.results().back().value() << " +- " << chkpt.results().back().error() << '\n';

    // perform the same integration in one go; the result is the same (even numerical)
    auto const result = hep::plain(
        hep::make_integrand<double>(
            function<double>,
            1,
            hep::make_dist_params<double>(10, 0.0, 1.0, "distribution #1")
        ),
        std::vector<std::size_t>(10, 1000)
    ).results().back();

    std::cout << "result using a single integration with ten iterations:\n    ";

    // print the result using enough precision
    std::cout << std::scientific << std::setprecision(16)
        << chkpt.results().back().value() << " +- " << chkpt.results().back().error() << '\n';
}
