#include <Eigen/Core>
#include <vector>
#include <random>

#ifndef WORLD_HPP
#define WORLD_HPP

// Not sure if I want this to be a pure singleton yet
class World {
    public: 
        
        explicit World(const int& ni, const int& nj, const int& nk); // Explicit constructor
        ~World() = default; // Destructor

        World(const World& other) = delete; // Delete copy consturctor
        World& operator=(const World& other) = delete; // Delete copy assignment operator 
        World(World&& other) = delete;// Delete move constructor

        // Getters
        const Eigen::Vector3d XtoL(const Eigen::Vector3d& x);
        const Eigen::Vector3d get_dh();
        const Eigen::Vector3d get_origin();
        const Field<double>& get_node_volumes();

        // Setters
        void set_extents(const double& x1, const double& y1, const double& z1,
                        const double& x2, const double& y2, const double& z2);

        void set_time(double dt, int num_ts);

        void compute_node_volumes();
        void compute_charge_density();

        void particle_push();

        const Eigen::Vector3i _nn; // Number of cells as an Eigen vector
        const int _ni, _nj, _nk; // Number of cells along each dimension

        // Fields
        class Field<double> _phi; // Electric potential
        class Field<double> _rho; // Charge density
        class Field3 _E; // Electric field

    protected:

        // The origin is always in the bottom corner of the simulation domain. It is NOT centered in a box.
        Eigen::Vector3d _x0; // Origin
        Eigen::Vector3d _dh; // Cell size along each dimension [dx, dy, dz]
        Eigen::Vector3d _xmax; // Maximum coordinate along each dimension
        Eigen::Vector3d _xc; // Domain centroid
        class Field<double> _node_vol; // Volume associated with each grid node

        std::vector<class Species*> _species;

        double _dt;
        int _num_timesteps;

};

// Making a choice here to have Rnd managed by the World header. Might revisit this later.
// I do want this to be a singleton.
class Rnd {
    public:
        Rnd() : _mt_gen{std::random_device()()}, _rnd_dist{0, 1.0} {};
        ~Rnd();

        Rnd(const Rnd& other) = delete; // Delete copy constructor
        Rnd& operator=(const Rnd& other) = delete; // Delete copy assignment operator
        Rnd(Rnd&& other) = delete; // Delete move constructor
        
        double operator()() {
            return _rnd_dist(_mt_gen);
        }; 

    private:
        std::mt19937 _mt_gen; // Random number generator
        std::uniform_real_distribution<double> _rnd_dist; // Uniform distribution between 0 and 1
};

extern Rnd global_rnd; 

#endif 