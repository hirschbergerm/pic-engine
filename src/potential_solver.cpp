#include "potential_solver.hpp"
#include "constants.hpp"
#include <iostream>

PotentialSolver::PotentialSolver(World& world, const unsigned int& max_solver_it, const double& tolerance) :
 _world(world),
 _max_solver_it(max_solver_it),
 _tolerance(tolerance) 
 {} 

PotentialSolver::~PotentialSolver() {}

/**
 * @brief Solves Poisson's equation using the Gauss-Siedel algorithm with successive over-relaxation (SOR).
 * 
 * @param phi The electric potential field to be solved (modified in place).
 */
auto PotentialSolver::solve() -> bool {

    Field<double>& phi = _world._phi;
    Field<double>& rho = _world._rho;
    const Eigen::Vector3i& nn = _world._nn;
    const Eigen::Vector3d& dh = _world.get_dh();
    
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

/** * @brief Computes the electric field from the electric potential using central differences.
 * 
 * @param phi The electric potential field (input).
 * @param E The electric field to be computed (modified in place).
*/
auto PotentialSolver::computeElectricField() -> void {
 
    Field<double>& phi = _world._phi;
    Field3& E = _world._E;

    const Eigen::Vector3d& dh = _world.get_dh();

    for (int i=0; i<_world._nn[0]; i++) {
        for (int j=0; j<_world._nn[1]; j++) {
            for (int k=0; k<_world._nn[2]; k++) {
                // Get references to field components to update
                double& Ex = E.x(i,j,k);
                double& Ey = E.y(i,j,k);
                double& Ez = E.z(i,j,k);

                // x component
                if (i==0) {
                    Ex = -(-3*phi(i,j,k) + 4*phi(i+1,j,k) - phi(i+2,j,k)) / (2.0 * dh[0]); // forward difference
                } else if (i==_world._nn[0]-1) {
                    Ex = -(phi(i-2,j,k) - 4*phi(i-1,j,k) + 3*phi(i,j,k)) / (2.0 * dh[0]); // backward difference
                } else {
                    Ex = -(phi(i+1,j,k) - phi(i-1,j,k)) / (2.0 * dh[0]); // central difference
                }

                if (j==0) {
                    Ey = -(-3*phi(i,j,k) + 4*phi(i,j+1,k) - phi(i,j+2,k)) / (2.0 * dh[1]); // forward difference
                } else if (j==_world._nn[1]-1) {
                    Ey = -(phi(i,j-2,k) - 4*phi(i,j-1,k) + 3*phi(i,j,k)) / (2.0 * dh[1]); // backward difference
                } else {
                    Ey = -(phi(i,j+1,k) - phi(i,j-1,k)) / (2.0 * dh[1]); // central difference
                }

                if (k==0) {
                    Ez = -(-3*phi(i,j,k) + 4*phi(i,j,k+1) - phi(i,j,k+2)) / (2.0 * dh[2]); // forward difference
                } else if (k==_world._nn[2]-1) {
                    Ez = -(phi(i,j,k-2) - 4*phi(i,j,k-1) + 3*phi(i,j,k)) / (2.0 * dh[2]); // backward difference
                } else {
                    Ez = -(phi(i,j,k+1) - phi(i,j,k-1)) / (2.0 * dh[2]); // central difference
                }

            }
        }
    }   
}