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

/**
 * @brief Load particles uniformly in a box.
 * @remarks Particles are loaded with zero velocity. This function modifies the internal storage of the Species object.
 * 
 * @param box_min Minimum corner of the box.
 * @param box_max Maximum corner of the box.
 * @param num_density Number density of particles to load (in real units).
 * @param num_particles Number of simulation particles to load.
 * 
 * @returns void
 */
void Species::load_particles_box(Eigen::Vector3d& box_min, Eigen::Vector3d box_max, double num_density, int num_sim_particles) {
    
    // Calculate the box volume
    double box_volume = (box_max[0] - box_min[0]) * (box_max[1] - box_min[1]) * (box_max[2] - box_min[2]);
    double num_real = num_density * box_volume; // Calculate the number of real particles that would be in the domain
    double mpw = num_sim_particles / num_real; // Calculate the macroparticle weight

    // Resize the storage vectors
    resize_storage(num_sim_particles);

    // Initialize particles 
    for (int p = 0; p < num_sim_particles; p++) {
        // sample position
        _x[p] = box_min[0] + global_rnd() * (box_max[0] - box_min[0]);
        _y[p] = box_min[1] + global_rnd() * (box_max[1] - box_min[1]);
        _z[p] = box_min[2] + global_rnd() * (box_max[2] - box_min[2]);

        // sample velocity
        _vx[p] = 0.0;
        _vy[p] = 0.0;
        _vz[p] = 0.0;

        // set macroparticle weight
        _mpwt[p] = mpw;
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