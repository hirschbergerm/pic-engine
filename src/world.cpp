#include "world.hpp"
#include "species.hpp"
#include "field.hpp"

/**
 * @brief Explicit constructor for the World class.
 */
explicit World::World(const int& ni, const int& nj, const int& nk) : 
    _ni(ni), 
    _nj(nj), 
    _nk(nk), 
    _nn{ni, nj, nk},
    _phi(ni, nj, nk),
    _rho(ni, nj, nk),
    _E(ni, nj, nk),
    _node_vol(ni, nj, nk)
    {}

World::~World() {

}

/**
 * @brief Converts from physical coordinates to logical (cell-centered) coordinates.
 * 
 * @param x Vector of physical coordinates [x,y,z].
 * @return const Eigen::Vector3d of logical coordinates [lx, ly, lz].
 */
const Eigen::Vector3d World::XtoL(const Eigen::Vector3d& x) {
    Eigen::Vector3d l;
    l[0] = (x[0] - _x0[0]) / _dh[0];
    l[1] = (x[1] - _x0[1]) / _dh[1];
    l[2] = (x[2] - _x0[2]) / _dh[2];
    return l;
}

const Eigen::Vector3d World::get_dh() {
    return _dh;
}

const Eigen::Vector3d World::get_origin() {
    return _x0;
}

const Field& World::get_node_volumes() {
    return _node_vol;
}

void World::set_extents(const double& x1, const double& y1, const double& z1, const double& x2, const double& y2, const double& z2) {

    // Set the origin coordinates
    _x0[0] = x1;
    _x0[1] = y1; 
    _x0[2] = z1;

    // Set the maximum coordiantes
    _xmax[0] = x2;
    _xmax[1] = y2;
    _xmax[2] = z2;

    // Compute the cell spacings
    _dh[0] = (_xmax[0] - _x0[0]) / _ni; // You will need static_cast here
    _dh[1] = (_xmax[1] - _x0[1]) / _nj;
    _dh[2] = (_xmax[2] - _x0[2]) / _nk;

    // Compute the domain centroid
    // This is solely needed for loading the electron population
    _xc[0] = 0.5 * (_xmax[0] + _x0[0]);
    _xc[1] = 0.5 * (_xmax[1] + _x0[1]);
    _xc[2] = 0.5 * (_xmax[2] + _x0[2]);
    
}

/**
 * @brief Computes the volume associated with each grid node and stores it in the internal node volume field.
 * 
 * @return void
 */
void World::compute_node_volumes() {

    for (int i = 0; i < _ni; i++) {
        for (int j = 0; j < _nj; j++) {
            for (int k = 0; k < _nk; k++) {
                double vol = _dh[0] * _dh[1] * _dh[2];
                if (i == 0 || i == _ni-1) vol *= 0.5;
                if (j == 0 || j == _nj-1) vol *= 0.5;
                if (k == 0 || k == _nk-1) vol *= 0.5;
                _node_vol(i, j, k) = vol;
            }
        }
    }
}

void World::compute_charge_density() {

    _rho = 0.0; // Initialize the charge density field to zero

    for (const auto& species : _species) {
        if (species->_charge == 0.0) continue; // Skip neutral species
        _rho += species->_den * species->_charge; // Accumulate charge density from each species

    }

}
