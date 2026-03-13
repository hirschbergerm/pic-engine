#include "world.hpp" 
#include "species.hpp"
#include <fstream>

namespace Output {

    void fields_output(World& world, std::vector<Species>& species);
    void screen_output(World &world, std::vector<Species>& species);
    void diagnostic_output(World &world, std::vector<Species>& species);

    std::ofstream f_diag;
}