#include "species.hpp"
#include "world.hpp"

/**
 * 
 */
explicit Species::Species(std::string name, double mass, double charge, World& world) : 
    _name(name), 
    _mass(mass), 
    _charge(charge),
    _x(0),
    _y(0),
    _z(0),
    _vx(0),
    _vy(0),
    _vz(0),
    _world(world),
    _den(_world._ni, _world._nj, _world._nk)
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
 * @brief Loading algorithm that ensures that the initial carge density is nearly zero throughout the grid. This algorithm minimizes the non-physical,
 * non-nuetrality caused by random particle loading.
 * 
 * @param
 */
void Species::load_particles_box_quiet_start(Eigen::Vector3d& box_min, Eigen::Vector3d box_max, double num_density, Eigen::Vector3i& num_sim_particles) {

    double box_volume = (box_max[0] - box_min[0]) * (box_max[1] - box_min[1]) * (box_max[2] - box_min[2]);
    int total_num_sim_particles = (num_sim_particles[0]-1) * (num_sim_particles[1]-1) * (num_sim_particles[2]-1);
    double num_real = num_density * box_volume; // Calculate the number of real particles that would be present
    double mpw = total_num_sim_particles / num_real; // Calculate the macroparticle weight

    // Compute spacing between particles in each dimension
    // You can imagine a box with num_sim_particles particles each a distance di apart.
    double di = (box_max[0] - box_min[0]) / (num_sim_particles[0]-1);
    double dj = (box_max[1] - box_min[1]) / (num_sim_particles[1]-1);
    double dk = (box_max[2] - box_min[2]) / (num_sim_particles[2]-1);

    // Resize the storage vectors
    resize_storage(total_num_sim_particles);

    // Begin loading particles on equally spaced grid
    for (int i = 0; i < num_sim_particles[0]; i++) {
        for (int j = 0; j < num_sim_particles[1]; j++) {
            for (int k = 0; k < num_sim_particles[2]; k++) {
                int p = i * num_sim_particles[0] + j * num_sim_particles[1] + k; // Generate a particle index. Order doesn't matter

                // Load particles in an equally spaced grid
                _x[p] = box_min[0] + i * di;
                _y[p] = box_min[1] + j * dj;
                _z[p] = box_min[2] + k * dk;

                // Shift particles on max faces back into the domain
                // Domain is up to but not including max faces x < [x_0, x_max)
                if (_x[p]==box_max[0]) _x[p]-=1e-4*di;
                if (_y[p]==box_max[1]) _y[p]-=1e-4*dj;
                if (_z[p]==box_max[2]) _z[p]-=1e-4*dk;

                // Assign weights
                // We scale the particles on the domain edges by 0.5 to account for that fact that 
                // edge node volumes are half of an internal node
                double w = 1; // relative weight 
                if (i==0 || i==num_sim_particles[0] - 1) w*= 0.5;
                if (j==0 || j==num_sim_particles[1] - 1) w*= 0.5;
                if (k==0 || k==num_sim_particles[2] - 1) w*= 0.5;

                // Sample velocities
                _vx[p] = 0.0;
                _vy[p] = 0.0;
                _vz[p] = 0.0;

                // Set macroparticle weight (mpw * relative scaling value)
                _mpwt[p] = mpw * w; 
            }
        }
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

/**
 * @brief Computes the number density field of the species based on the current particle positions and macroparticle weights. Modifies the internal 
 * density Field member of the Species object.
 * 
 * @return void
 */
void Species::compute_number_density() { 

    _den = 0.0; // 

    // Loop over all sim particles
    const int num_sim_particles = _x.size();

    for (int p = 0; p < num_sim_particles; p++) {
        // Convert to logical coordinates
        Eigen::Vector3d l = _world.XtoL({_x[p], _y[p], _z[p]});

        // Scatter the density value to the grid
        _den.scatter(l, _mpwt[p]);
        
    }

    _den /= _world.get_node_volumes(); // Divide by the node volumes to get number density

}