#include "world.hpp"
#include "species.hpp"
#include "constants.hpp"

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
    _node_vol(ni, nj, nk),
    _start_time(std::chrono::high_resolution_clock::now())
    {}

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

Eigen::Vector3d World::get_dh() const {
    return _dh;
}

Eigen::Vector3d World::get_origin() const {
    return _x0;
}

Eigen::Vector3d World::get_xmax() const {
    return _xmax;
}

const Field<double>& World::get_node_volumes() const {
    return _node_vol;
}

double World::get_wall_time() const {
    auto now = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = now - _start_time;
    return elapsed.count();
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

void World::compute_charge_density(const std::vector<Species*>& species) {

    _rho = 0.0; // Initialize the charge density field to zero

    for (const auto& sp : species) {
        if (sp->_charge == 0.0) continue; // Skip neutral species

        sp->_den *= sp->_charge; // Scale density by charge to get charge density contribution from this species
        _rho += sp->_den; // Add this species' contribution to the total charge density

    }

} 

/**
 * @brief Advances the simulation time by one timestep. If the simulation has reached the end, it returns false to signal that the simulation is over.
 */
bool::World::advance_time() {
    
    // Check that the user doesn't advance time beyond the number of timesteps specified. 
    if (_current_timestep >= _num_timesteps) {
        return false; 
    }

    _current_time += _dt;
    _current_timestep++;

    return _current_timestep < _num_timesteps;

}

/**
 * @brief Calculates the total potential energy of the system
 */
void World::get_potential_energy(double& pe) const {

    pe = 0.0;
    Eigen::Vector3d ef{0.0};

    for (int i = 0; i < _ni; i++) {
        for (int j = 0; j < _nj; j++) {
            for (int k = 0; k < _nk; k++) {
                ef = _E(i, j, k); // Get field value at node i,j,k 
                double ef2 = ef[0]* ef[0] + ef[1]*ef[1] + ef[2]*ef[2]; // Take L2 norm of field at this point
                pe += ef2*_node_vol(i, j, k); // Add contribution to potential energy using the node volume as a weight
            }
        }
    }

    pe *= 0.5 * Const::EPS_0; // Scale by proper physical factor

}

