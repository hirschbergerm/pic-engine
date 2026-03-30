#pragma once
#ifndef SPECIES_HPP
#define SPECIES_HPP

#include <string>
#include "field.hpp"

class Species {
    public:

        // Constructors and Destructor
        explicit Species(std::string name, double mass, double charge, class World& world);
        ~Species() = default;

        // Public Methods
        void load_particles_box(Eigen::Vector3d& box_min, Eigen::Vector3d box_max, double num_density, int num_particles); 
        void load_particles_box_quiet_start(Eigen::Vector3d box_min, Eigen::Vector3d box_max, double num_density, Eigen::Vector3i num_sim_particles);
        void compute_number_density();

        void push_particles();

        void get_kinetic_energy(double& ke) const;
        void get_momentum(double& px, double& py, double& pz) const;
        void get_macroparticle_count(size_t& mp_count) const;
        void get_real_count(double& real_count) const;
        const std::string _name;
        const double _charge; 
        const double _mass;

    private:

        void resize_storage(const int& new_size);

        // Storage for the particle velocities and positions
        // We are using a Data-Oriented Design approach, so NO arrays of structs
        std::vector<double> _x; // x coordinates of all the particles
        std::vector<double> _y; // y coordinates of all the particles
        std::vector<double> _z; // z coordinates of all the particles

        std::vector<double> _vx; // x velocities of all the particles
        std::vector<double> _vy; // y velocities of all the particles
        std::vector<double> _vz; // z velocities of all the particles

        std::vector<double> _mpwt; // macroparticle weights of all the particles

        class World& _world; // Reference to the world object the species belongs to
        // Using a reference here to avoid circular dependencies and because Species should not exist without a World

    public:
        Field<double> _den; // Number density field of the species, defined on the grid nodes. This is public for now for easy access in output routines, but we may want to make this private and add a getter function later.
};

#endif // SPECIES_HPP