#include "world.hpp"
#include "species.hpp"

/**
 * 
 */
explicit Species::Species(std::string name, double mass, double charge, World& world) : 
    name(name), 
    mass(mass), 
    charge(charge),
    _x(0),
    _y(0),
    _z(0),
    _vx(0),
    _vy(0),
    _vz(0),
    _world(world),
    den(_world._ni, _world._nj, _world._nk)
{}

void Species::addParticle(const Eigen::Vector3d& pos, Eigen::Vector3d& vel, double mpwt) {
    
}

void Species::loadParticlesBox(Eigen::Vector3d& box_min, Eigen::Vector3d box_max, double num_density, int num_sim_particles) {
    
    // Calculate the box volume
    double box_volume = (box_max[0] - box_min[0]) * (box_max[1] - box_min[1]) * (box_max[2] - box_min[2]);
    double num_real = num_density * box_volume; // Calculate the number of real particles that would be in the domain
    double mpw = num_sim_particles / num_real; // Calculate the macroparticle weight

    // Resize the storage vectors
    resize_storage(num_sim_particles);

    // Initialize particles 
    for (int p = 0; p < num_sim_particles; p++) {
        // sample position
    }

}

/**
 * @brief Resize the internal storage vectors to accommodate the passed number (new_size) of particles
 */
void Species::resize_storage(const int& new_size) {
    _x.conservativeResize(new_size);
    _y.conservativeResize(new_size);
    _z.conservativeResize(new_size);

    _vx.conservativeResize(new_size);
    _vy.conservativeResize(new_size);
    _vz.conservativeResize(new_size);

    _mpwt.conservativeResize(new_size);
}