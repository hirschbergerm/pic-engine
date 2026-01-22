#include "field.hpp"
#include <Eigen/Core>

#ifndef POTENTIAL_SOLVER_HPP
#define POTENTIAL_SOLVER_HPP

class PotentialSolver {
    public:
        explicit PotentialSolver(const unsigned int& max_solver_it, const double& tolerance);
        ~PotentialSolver();

        auto solve(const Field& phi, const Field& rho, const Eigen::Vector3i& nn, const Eigen::Vector3d& dh) -> bool; // Solve Poisson's Equation

        auto computeElectricField(const Field& phi) -> void; // Calculates static EF = -grad(phi)

    private:
        // I disagree with the idea of making world a private member
        unsigned int _max_solver_it;
        double _tolerance;

};

#endif // POTENTIAL_SOLVER_HPP