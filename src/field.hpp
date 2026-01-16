#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

#ifndef FIELD_HPP
#define FIELD_HPP

class Field {
    public:
        Field(const int& ni, const int& nj, const int& nk);
        ~Field();

        // Create math operators for field
        void operator /=(const Field& other);
        Field operator+=(const Field& other); // Add two fields 

    private:
        Eigen::Tensor<double, 3> _data; // 3D field data

    protected:
        const int _ni, _nj, _nk;
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