#include <Eigen/Core>
#include <vector>
#include "field.hpp"

#ifndef WORLD_HPP
#define WORLD_HPP

// Not sure if I want this to be a pure singleton yet
class World {
    public: 
        
        explicit World(const int& ni, const int& nj, const int& nk); // Explicit constructor
        ~World(); // Destructor

        World(const World& other) = delete; // Delete copy consturctor
        World& operator=(const World& other) = delete; // Delete copy assignment operator 
        World(World&& other) = delete;// Delete the move constructor

        void setExtents(const double& x1, const double& y1, const double& z1,
                        const double& x2, const double& y2, const double& z2);

        const Eigen::Vector3i _nn; // Number of cells as an Eigen vector
        const int _ni, _nj, _nk; // Number of cells along each dimension

    private:

        // The origin is always in the bottom corner of the simulation domain. It is NOT centered in a box.
        Eigen::Vector3d _x0; // Origin
        Eigen::Vector3d _dh; // Cell size along each dimension [dx, dy, dz]
        Eigen::Vector3d _xmax; // Maximum coordinate along each dimension
        Eigen::Vector3d _xc; // Domain centroid

        // Fields
        Field _phi; // Electric potential
        Field _rho; // Charge density
        Field3 _E; // Electric field

        std::vector<class Species*> _species;

};

#endif 