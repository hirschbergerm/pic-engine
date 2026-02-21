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

void Field::operator*(const double& value) {

    _data = _data * value; // Multiply all field values by the passed value
}

double& Field::operator()(const int& i, const int& j, const int& k) {

    return _data(i, j, k);

}

/**
 * @brief Checks if the logical coordinate l is within bounds of the field.
 */
bool Field::in_bound_logical(const Eigen::Vector3d& l) const {

    return (l(0) >= 0 && l(0) <= _ni-1 && l(1) >= 0 && l(1) <= _nj-1 && l(2) >= 0 && l(2) <= _nk-1);

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
    
    // Logical coordinate is an index with a fractional part
    // 1.75 -> i=1, di=0.75

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

/*
 * @brief Gathers a field value at logical coordinate (l) from the adjacent cell nodes. Modifies the passed value.
 * 
 * @param l The logical coordinate to gather from (3D vector).
 * @param value The field value to gather into.
 * @return void
 */
void Field::gather(const Eigen::Vector3d& l, double& value) {

    // Check that the logical coordinates are within bounds
    if (!in_bound_logical(l)) {
        return;
    }

    // Logical coordinate is an index with a fractional part
    // 1.75 -> i=1, di=0.75

    // Get the node indices
    int i = static_cast<int>(std::floor(l(0)));
    int j = static_cast<int>(std::floor(l(1)));
    int k = static_cast<int>(std::floor(l(2)));

    // Get the fractional distances
    double di = l(0) - i;
    double dj = l(1) - j;
    double dk = l(2) - k;

    // Interpolate data from surrounding nodes
    // We add the contributions from each of the surrounding nodes weighted by their distance to the logical coordinate
    value = (1-di)*(1-dj)*(1-dk)*_data(i,j,k) +
            (di)*(1-dj)*(1-dk)*_data(i+1,j,k) +
            (1-di)*(dj)*(1-dk)*_data(i,j+1,k) + 
            (di)*(dj)*(1-dk)*_data(i+1,j+1,k) +
            (1-di)*(1-dj)*(dk)*_data(i,j,k+1) +
            (di)*(1-dj)*(dk)*_data(i+1,j,k+1) +
            (1-di)*(dj)*(dk)*_data(i,j+1,k+1) +
            (di)*(dj)*(dk)*_data(i+1,j+1,k+1);

    return; 
}