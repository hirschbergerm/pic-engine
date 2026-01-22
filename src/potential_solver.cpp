#include "potential_solver.hpp"
#include "constants.hpp"

PotentialSolver::PotentialSolver(const unsigned int& max_solver_it, const double& tolerance) :
 _max_solver_it(max_solver_it),
 _tolerance(tolerance) 
 {} 

PotentialSolver::~PotentialSolver() {}

/**
 * @brief Solves Poisson's equation using the Gauss-Siedel algorithm with successive over-relaxation (SOR).
 * 
 * @param phi The electric potential field to be solved (modified in place).
 */
auto PotentialSolver::solve(const Field& phi, const Field& rho, const Eigen::Vector3i& nn, const Eigen::Vector3d& dh) -> bool {
    
    // precompute (1/dx^2)
    double idx2 = 1.0/(dh[0]*dh[0]);
    double idy2 = 1.0/(dh[1]*dh[1]);
    double idz2 = 1.0/(dh[2]*dh[2]);
    
    double L2 = 0; // L2 norm
    bool converged = false;

    // main solver loop
    for (unsigned int it = 0; it < _max_solver_it; it++) {
        for(int i=1; i < nn[0]-1; i++) {
            for(int j=1; j < nn[1]-1; j++) {
                for(int k=1; k < nn[2]-1; k++) {
                    // update step for internal node
                    double phi_new = (rho(i,j,k))
                }
            }
        }
    }





    return true; // Placeholder return value
}