#include "hep/mc.hpp"

#include <cstddef>
#include <fstream>
#include <random>
#include <vector>

template <typename T>
T function(hep::mc_point<T> const& point, hep::projector<T>& projector)
{
    T const x = point.point().at(0);
    T const v = T(2.0) * x;

    projector.add(0, x, v);

    return v;
}

int main()
{
    std::cout << ">>> loading checkpoint from `chkpt2` ...\n\n";

    std::ifstream in("chkpt2");
    auto const chkpt_from_disk = hep::make_plain_chkpt<double, std::mt19937>(in);
    in.close();

    // resume the integration with 5 additional iterations
    auto const chkpt = hep::plain(
        hep::make_integrand<double>(
            function<double>,
            1,
            hep::make_dist_params<double>(10, 0.0, 1.0, "distribution #1")
        ),
        std::vector<std::size_t>(5, 100000),
        chkpt_from_disk
    );

    std::cout << ">>> overwriting checkpoint ...\n";

    std::ofstream out("chkpt2");
    chkpt.serialize(out);

    std::cout << ">>> - restart this program as often as you like to improve the result\n";
    std::cout << ">>> - remove the file `chkpt2` to reset\n";
}
