#pragma once
#ifndef FIELD_HPP
#define FIELD_HPP

#include <Eigen/Core>
#include <ostream>

template <typename T>
class Field {
    public:

        Field<T>(size_t ni, size_t nj, size_t nk);
        ~Field<T>() = default;
        Field<T>(Field<T>&& other) noexcept; // Move constructor
        Field<T>& operator=(Field<T>&& other) noexcept; // Move assignment operator

        // Operator overloads
        void operator=(const T& value); // Set all field values to the assigned value

        // Elementwise Field operations. The majority of these only modify Field in-place for performance.
        Field<T>& operator*=(const T& value); // Multiply all field values by the passed value
        Field<T>& operator/=(const T& value); // Divide all field values by the passed value
        Field<T>& operator+=(const Field<T>& other); // Elementwise addition with another Field
        Field<T>& operator-=(const Field<T>& other); // Elementwise subtraction with another Field
        
        Field<T>& operator/=(const Field<T>& other); // Elementwise division by another Field
        
        // Read-Write Field access operator for value at (i,j,k)
        inline T& operator()(size_t i, size_t j, size_t k) {
            return _data[index(i,j,k)];
        } 

        // Read-Only Field access operator for value at (i,j,k)
        inline const T& operator()(size_t i, size_t j, size_t k) const { 
            return _data[index(i,j,k)];
        }

        void scatter(const Eigen::Vector3d& l, const T& value); // Scatter a field value at the logical coordinate l
        void gather(const Eigen::Vector3d& l, T& value); // Gather a field value at the logical coordinate l

        bool in_bound_logical(const Eigen::Vector3d& l) const; // Check if the logical coordinate l is within bounds of the field

        // Const component access operators.
        inline const size_t ni() const{
            return _ni;
        }

        inline const size_t nj() const{
            return _nj;
        }

        inline const size_t nk() const{
            return _nk;
        }
        
    private:
        std::vector<T> _data; // 3D field data

    protected:

        inline size_t index(const size_t i, const size_t j, const size_t k) const {   
            assert(i < _ni && j < _nj && k < _nk);
            return i + _ni * (j + _nj * k);
        }

        const size_t _ni, _nj, _nk; // Number of nodes in each direction

};

template class Field<double>;

// Stream output operator declaration
template <typename T>
std::ostream& operator<<(std::ostream& out, const Field<T>& field);

class Field3 {
    public: 
        Field3(size_t ni, size_t nj, size_t nk);
        ~Field3() = default;

    
        bool in_bound_logical(const Eigen::Vector3d& l) const;

        void scatter(const Eigen::Vector3d& l, const Eigen::Vector3d& value);
        void gather(const Eigen::Vector3d& l, Eigen::Vector3d& value);

        // Mutable access operator 
        auto operator()(size_t i, size_t j, size_t k) {
            return VecProxy{x(i,j,k), y(i,j,k), z(i,j,k)};
        }

        // Const access operator
        Eigen::Vector3d operator()(size_t i, size_t j, size_t k) const {
            return Eigen::Vector3d(_dataX[index(i,j,k)], _dataY[index(i,j,k)], _dataZ[index(i,j,k)]);
        }

        // Mutable component access operators.
        inline double& x(size_t i, size_t j, size_t k) { 
            return _dataX[index(i,j,k)];
        };

        inline double& y(size_t i, size_t j, size_t k) { 
            return _dataY[index(i,j,k)];
        };

        inline double& z(size_t i, size_t j, size_t k) { 
            return _dataZ[index(i,j,k)];
        };

        // Const component access operators.
        inline const size_t ni() const{
            return _ni;
        }

        inline const size_t nj() const{
            return _nj;
        }

        inline const size_t nk() const{
            return _nk;
        }

    private:
        std::vector<double> _dataX;
        std::vector<double> _dataY;
        std::vector<double> _dataZ;

    protected:

        // VecProxy struct is private so that users never see it. Access operators handle automatic conversion to Eigen::Vector3d
        struct VecProxy {
            double& x;
            double& y;
            double& z;

            // This is a type-casting operator basically
            // When you try to assign a VecProxy object to an Eigen::Vector3d, this function will be called
            operator Eigen::Vector3d() const {
                return Eigen::Vector3d(x, y, z);
            };

            // Return a mutable vector proxy that 
            VecProxy& operator=(const Eigen::Vector3d& value) {
                x = value[0];
                y = value[1];
                z = value[2];
                return *this;
            }
        };

        inline size_t index(const size_t i, const size_t j, const size_t k) const {
            assert(i < _ni && j < _nj && k < _nk);
            return i + _ni * (j + _nj * k);
        }

        const size_t _ni, _nj, _nk; // Number of nodes in each direction
}; 

// Declare output stream operator for Field3
std::ostream& operator<<(std::ostream& out, const Field3& field);

#endif // FIELD_HPP