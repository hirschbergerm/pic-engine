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

    private:

        const int nn[3]; // Number of cells along each dimension [ni, nj, nk]
        const int ni, nj, nk;

        // The origin is always in the bottom corner of the simulation domain. It is NOT centered in a box.
        double x0[3]; // Origin
        double dh[3]; // Cell size along each dimension [dx, dy, dz]
        double xmax[3]; // Maximum coordinate along each dimension
        double xc[3]; // Domain centroid
};

#endif 