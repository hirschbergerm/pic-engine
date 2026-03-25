#pragma once
#ifndef OUTPUT_HPP
#define OUTPUT_HPP

#include <fstream>
#include "world.hpp"
#include "species.hpp"

namespace Output {

    void fields_output(World& world, std::vector<Species*>& species);
    void screen_output(World &world, std::vector<Species*>& species);
    void diagnostic_output(World &world, std::vector<Species*>& species);

    extern std::ofstream f_diag;
}

#endif // OUTPUT_HPP