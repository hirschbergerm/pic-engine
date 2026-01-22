#include <string>
#include "field.hpp"

#ifndef SPECIES_HPP
#define SPECIES_HPP

class Species {
    public:

        // Constructors and Destructor
        explicit Species(std::string name, double mass, double charge, class World& world);
        ~Species() = default;
        
        //Species(const Species& other) = default; // Copy Constructor

        // Public Methods
        void addParticle(const Eigen::Vector3d& pos, Eigen::Vector3d& vel, double mpwt);
        void loadParticlesBox(Eigen::Vector3d& box_min, Eigen::Vector3d box_max, double num_density, int num_particles); 

        const std::string name;
        const double charge; 
        const double mass;

        Field den; // Number density field of the species

    private:

        void resize_storage(const int& new_size);

        // Storage for the particle velocities and positions
        // We are using a Data-Oriented Design approach, so NO arrays of structs
        Eigen::VectorXd _x; // x coordinates of all the particles
        Eigen::VectorXd _y; // y coordinates of all the particles
        Eigen::VectorXd _z; // z coordinates of all the particles

        Eigen::VectorXd _vx; // x velocities of all the particles
        Eigen::VectorXd _vy; // y velocities of all the particles
        Eigen::VectorXd _vz; // z velocities of all the particles

        Eigen::VectorXd _mpwt; // macroparticle weights of all the particles

        class World& _world; // Reference to the world object the species belongs to
        // Using a reference here to avoid circular dependencies and because Species should not exist without a World
};

#endif // SPECIES_HPP