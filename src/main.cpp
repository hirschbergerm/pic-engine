#include <iostream>
#include <string>
#include <vector>

#include "constants.hpp"
#include "world.hpp"
#include "field.hpp"
#include "species.hpp"
#include "potential_solver.hpp"
#include "output.hpp"

int main(int argc, char* argv[]) {
    
    // Initialize simulation domain
    World world(21, 21, 21);
    world.set_extents(-0.1, -0.1, 0, 0.1, 0.1, 0.2);
    world.compute_node_volumes(); // Compute node volumes for the domain, needed for charge density calculation and output
    world.set_time(2e-10, 200);
    std::cout << "Created World" << std::endl;

    // Set up particle species
    std::vector<Species*> species;
    species.emplace_back(new Species("O+", Const::mp, Const::e, world)); // Proton species
    species.emplace_back(new Species("e-", Const::me, -Const::e, world)); // Electron species
    std::cout << "Created Species" << std::endl;

    // Initialize potential solver and solve initial potential
    PotentialSolver p_solver(world, 10000, 1e-4);
    p_solver.solve();

    // Compute initial electric field from potential
    p_solver.compute_electric_field();

    // Load particles specifying number of simulation particles in each dimension
    species[0]->load_particles_box_quiet_start(world.get_origin(), world.get_xmax(), 1e11, {21, 21, 21});
    species[1]->load_particles_box_quiet_start(world.get_origin(), world.get_xcenter(), 1e11, {41, 41, 41});

    // Compute initial number density fields for each species
    for (auto& sp : species) {
        sp->compute_number_density();
    }

    // Main simulation loop
    while (world.advance_time()) {
        // Compute charge density on the grid based on current particle positions
        world.compute_charge_density(species);

        // Solve potential based on this charge density
        p_solver.solve();

        // Compute electric field from potential
        p_solver.compute_electric_field();

        // Push paticles based on electric field
        for (auto& sp : species) {
            sp->push_particles();
            sp->compute_number_density();
        }

        // Screen and file output
        Output::screen_output(world, species);
        Output::diagnostic_output(world, species);

        // periodically write output files for visualization
        if (world.get_timestep() % 10 == 0 || world.is_last_timestep()) {
            Output::fields_output(world, species);
        }
    }

    std::cout << "Simulation took " << world.get_wall_time() << "seconds"<< std::endl;

    return 0;

}