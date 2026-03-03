#include <unsupported/Eigen/CXX11/Tensor>

#ifndef FIELD_HPP
#define FIELD_HPP

template <typename T>
class Field {
    public:

        Field<T>(const int& ni, const int& nj, const int& nk);
        ~Field<T>() = default;
        Field<T>(Field<T>&& other); // Move constructor
        Field<T>& operator=(Field<T>&& other) noexcept; // Move assignment operator

        // Operator overloads
        void operator=(const T& value); // Set all field values to the assigned value
        void operator*(const T& value); // Multiply all field values by the passed value
        Field<T>& operator*=(const T& value); // Multiply all field values by the passed value
        Field<T>& operator/=(const T& value); // Divide all field values by the passed value
        Field<T>& operator+=(const Field<T>& other); // Elementwise addition with another Field
        Field<T>& operator/=(const Field<T>& other);; // Elementwise division by another Field
        
        T& operator()(const int& i, const int& j, const int& k); // Read-Write Field access operator for value at (i,j,k)

        void scatter(const Eigen::Vector3d& l, const T& value); // Scatter a field value at the logical coordinate l
        T& gather(const Eigen::Vector3d& l, T& value); // Gather a field value at the logical coordinate l

        bool in_bound_logical(const Eigen::Vector3d& l) const; // Check if the logical coordinate l is within bounds of the field

    private:
        Eigen::Tensor<T, 3> _data; // 3D field data

    protected:
        const int _ni, _nj, _nk; // Number of nodes in each direction

};


class Field3 {
    public: 
        Field3(const int& ni, const int& nj, const int& nk);
        ~Field3() = default;

        bool in_bound_logical(const Eigen::Vector3d& l) const;

        void scatter(const Eigen::Vector3d& l, const Eigen::Vector3d& value);
        void gather(const Eigen::Vector3d& l, Eigen::Vector3d& value);

        Eigen::Vector3d& operator()(const int& i, const int& j, const int& k); // Read-Write Field access operator for value at (i,j,k)

    private:
        Eigen::Tensor<double, 3> _dataX;
        Eigen::Tensor<double, 3> _dataY;
        Eigen::Tensor<double, 3> _dataZ;
        
    protected:
        const int _ni, _nj, _nk; // Number of nodes in each direction
}; 

#endif // FIELD_HPP