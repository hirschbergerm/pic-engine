#include "field.hpp"

Field::Field(const int& ni, const int& nj, const int& nk) : _ni(ni), _nj(nj), _nk(nk) {

    _data = Eigen::Tensor<double, 3>(ni, nj, nk);
    _data.setZero(); // Set all of the data to 0


}

Field::~Field() {

    //delete &(_data); Free the allocated memory for data
    
}

double& Field::operator()(const int& i, const int& j, const int& k) {
    return _data(i, j, k);
}

