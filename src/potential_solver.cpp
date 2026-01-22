#include "potential_solver.hpp"
#include "constants.hpp"
#include <iostream>

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
auto PotentialSolver::solve(Field& phi, Field& rho, const Eigen::Vector3i& nn, const Eigen::Vector3d& dh) -> bool {
    
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
                    double phi_new = (rho(i,j,k)/Const::EPS_0 +
                                      (phi(i+1,j,k) + phi(i-1,j,k)) * idx2 +
                                      (phi(i,j+1,k) + phi(i,j-1,k)) * idy2 +
                                      (phi(i,j,k+1) + phi(i,j,k-1)) * idz2) /
                                     (2.0 * (idx2 + idy2 + idz2));

                    // SOR step
                    phi(i,j,k) = phi(i,j,k) + 1.5 * (phi_new - phi(i,j,k));

                



                }
            }
        }

        if (it % 25 == 0) {
            double sum = 0;
            for(int i=1; i < nn[0]-1; i++) {
                for(int j=1; j < nn[1]-1; j++) {
                    for(int k=1; k < nn[2]-1; k++) {
                        // Calculate the residual
                        double R = -phi(i,j,k)*(2*idx2 + 2*idy2 + 2*idz2) +
                                    rho(i,j,k)/Const::EPS_0 + 
                                    idx2*(phi(i-1,j,k) + phi(i+1,j,k)) +
                                    idy2*(phi(i,j-1,k) + phi(i,j+1,k)) +
                                    idz2*(phi(i,j,k-1) + phi(i,j,k+1));

                        sum += R * R;
                    }
                }
            }

            L2 = std::sqrt(sum / ((nn[0])*(nn[1])*(nn[2])));
            if (L2 < _tolerance) {
                converged = true;
                break;
            }
        }
    }

    if (!converged) {
        std::cerr << "Gauss-Siedel Solver failed to converge. L2=" << L2 << " after " << _max_solver_it << " iterations." << std::endl;
    }

    return converged;
}