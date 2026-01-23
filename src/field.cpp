#include "field.hpp"
#include <Eigen/Core>


Field::Field(const int& ni, const int& nj, const int& nk) : _ni(ni), _nj(nj), _nk(nk) {

    _data = Eigen::Tensor<double, 3>(ni, nj, nk); // Initialize the data tensor to the correct size
    _data.setZero(); // Set all of the data to 0

}

Field::~Field() {

    //delete &(_data); Free the allocated memory for data
    
}

/**
 * @brief Sets all field values to the assigned value.
 * 
 * @param value The value to set all field elements to.
 * @return void
 */
void Field::operator=(const double& value) {

    _data.setConstant(value); // Set all of the Eigen tensor values to value

}

/**
 * @brief Divides all field values by the passed value.
 * 
 * @param value The value to divide all field elements by.
 * @return void
 */
void Field::operator/=(const double& value) {

    _data = _data / value; // Divide all field values by the passed value

}

/**
 * @brief Elementwise division of this field by another field.
 * 
 * @param other The other field to divide by.
 * @return void
 */
void Field::operator/=(const Field& other) {

    _data = _data / other._data;
}

void Field::operator+=(const Field& other) {

    _data = _data + other._data;
}

double& Field::operator()(const int& i, const int& j, const int& k) {

    return _data(i, j, k);

}

/**
 * @brief Scatters a field value at logical coordinate (l) to the adjacent cell nodes. Modifies the internal data of the Field.
 * @attention This is an ADDITIVE operation - the value is added to the existing field values at the nodes. You may need to reset your
 * field values before calling this method.
 * 
 * @param l The logical coordinate to scatter to (3D vector).
 * @param value The field value to scatter.
 * @return void
 */
void Field::scatter(const Eigen::Vector3d& l, const double& value) {

    // Check that the logical coordinates are within bounds
    if (l(0) < 0 || l(0) > _ni-1 || l(1) < 0 || l(1) > _nj-1 || l(2) < 0 || l(2) > _nk-1) {
        return; 
    }
    
    // Get the node indices
    int i = static_cast<int>(std::floor(l(0)));
    int j = static_cast<int>(std::floor(l(1)));
    int k = static_cast<int>(std::floor(l(2)));

    // Get the fractional distance to each node
    double di = l(0) - i;
    double dj = l(1) - j;
    double dk = l(2) - k;

    // Scatter data to the surrounding nodes
    // Bottom four nodes
    _data(i,j,k) += value * (1-di) * (1-dj) * (1-dk);
    _data(i+1,j,k) += value * (di) * (1-dj) * (1-dk);
    _data(i,j+1,k) += value * (1-di) * (dj) * (1-dk); 
    _data(i+1,j+1,k) += value * (di) * (dj) * (1-dk);
    // Top four nodes
    _data(i,j,k+1) += value * (1-di) * (1-dj) * (dk);
    _data(i+1,j,k+1) += value * (di) * (1-dj) * (dk);
    _data(i,j+1,k+1) += value * (1-di) * (dj) * (dk);
    _data(i+1,j+1,k+1) += value * (di) * (dj) * (dk);

}
