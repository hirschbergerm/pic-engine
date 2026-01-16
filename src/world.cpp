#include "world.hpp"

World::World(const int& ni, const int& nj, const int& nk) : 
    _ni(ni), 
    _nj(nj), 
    _nk(nk), 
    nn{ni, nj, nk},
    _phi(ni, nj, nk),
    _rho(ni, nj, nk),
    _E(ni, nj, nk) 
    {}

World::~World() {

}

void World::setExtents(const double& x1, const double& y1, const double& z1, const double& x2, const double& y2, const double& z2) {

    // Set the origin coordinates
    x0[0] = x1;
    x0[1] = y1; 
    x0[2] = z1;

    // Set the maximum coordiantes
    xmax[0] = x2;
    xmax[1] = y2;
    xmax[2] = z2;

    // Compute the cell spacings
    dh[0] = (xmax[0] - x0[0]) / _ni; // You will need static_cast here
    dh[1] = (xmax[1] - x0[1]) / _nj;
    dh[2] = (xmax[2] - x0[2]) / _nk;

    // Compute the domain centroid
    // This is solely needed for loading the electron population
    xc[0] = 0.5 * (xmax[0] + x0[0]);
    xc[1] = 0.5 * (xmax[1] + x0[1]);
    xc[2] = 0.5 * (xmax[2] + x0[2]);
    
}
