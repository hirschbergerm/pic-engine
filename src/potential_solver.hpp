#include "world.hpp"
#include "field.hpp"
#include <Eigen/Core>

#ifndef POTENTIAL_SOLVER_HPP
#define POTENTIAL_SOLVER_HPP

class PotentialSolver {
    public:
        explicit PotentialSolver(World& world, const unsigned int& max_solver_it, const double& tolerance);
        ~PotentialSolver();

        auto solve() -> bool; // Solve Poisson's Equation

        auto computeElectricField() -> void; // Calculates static EF = -grad(phi)

    private:
        // I disagree with the idea of making world a private member
        unsigned int _max_solver_it;
        double _tolerance;
        World& _world;

};

#endif // POTENTIAL_SOLVER_HPP