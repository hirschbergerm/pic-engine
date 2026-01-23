#include <unsupported/Eigen/CXX11/Tensor>

#ifndef FIELD_HPP
#define FIELD_HPP

class Field {
    public:
        Field(const int& ni, const int& nj, const int& nk);
        ~Field();

        // Operator overloads
        void operator=(const double& value); // Set all field values to the assigned value
        void operator/=(const double& value); // Divide all field values by the passed value
        void operator/=(const Field& other);; // Elementwise division by another Field
        void operator+=(const Field& other); // Elementwise addition with another Field
        double& operator()(const int& i, const int& j, const int& k); // Read-Write Field access operator for value at (i,j,k)

        void scatter(const Eigen::Vector3d& l, const double& value); // Scatter a field value at the logical coordinate l
        void gather(const Eigen::Vector3d& l, double& value); // Gather a field value at the logical coordinate l

    private:
        Eigen::Tensor<double, 3> _data; // 3D field data

    protected:
        const int _ni, _nj, _nk; // Number of nodes in each direction

};


class Field3 : Field {
    public: 
        Field3(const int& ni, const int& nj, const int& nk);
        ~Field3();

    private:
        Eigen::Tensor<double, 3> _dataX;
        Eigen::Tensor<double, 3> _dataY;
        Eigen::Tensor<double, 3> _dataZ;
        
}; 

#endif // FIELD_HPP