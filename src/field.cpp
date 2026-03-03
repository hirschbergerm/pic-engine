#include "field.hpp"
#include <Eigen/Core>

template <typename T>
Field<T>::Field(const int& ni, const int& nj, const int& nk) : _ni(ni), _nj(nj), _nk(nk) {

    _data = Eigen::Tensor<T, 3>(ni, nj, nk); // Initialize the data tensor to the correct size
    _data.setZero(); // Set all of the data to 0

}

/** 
* @brief Move constructor for the Field class. Moves the data from the other field to this field.
* 
* @param other The other field to move from.
*/
template <typename T>
Field<T>::Field(Field<T>&& other) : 
    _ni{other._ni},
    _nj{other._nj},
    _nk{other._nk} {

    _data = std::move(other._data); // Move the data from the other field to this field

}

/**
 * @brief Move assignment operator for the Field class. Moves the data from the other field to this field.
 * 
 * @param other The other field to move from.
 * @return Field& A reference to this field after the move assignment.
 */
template <typename T>
Field<T>& Field<T>::operator=(Field<T>&& other) noexcept {

    _data = std::move(other._data); // Move the data from the other field to this field

    return *this;
}

/**
 * @brief Value assignment operator. Sets all field values to the assigned value.
 * 
 * @param value The value to set all field elements to.
 * @return void
 */
template <typename T>
void Field<T>::operator=(const T& value) {

    _data.setConstant(value); // Set all of the Eigen tensor values to value

}

/**
 * @brief Divides all field values by the passed value.
 * 
 * @param value The value to divide all field elements by.
 * @return void
 */
template <typename T>
Field<T>& Field<T>::operator/=(const T& value) {

    _data = _data / value; // Divide all field values by the passed value
    return *this;

}

/**
 * @brief Elementwise division of this field by another field.
 * 
 * @param other The other field to divide by.
 * @return void
 */
template <typename T>
Field<T>& Field<T>::operator/=(const Field<T>& other) {

    _data = _data / other._data;
    return *this;

}

/**
 * @brief Addition chain operator.Allows for elementwise addition of this field with another field.
 * 
 * @param other The other field to add to this field.
 * @return *this A reference to this field after the addition.
 */
template <typename T>
Field<T>& Field<T>::operator+=(const Field<T>& other) {

    _data = _data + other._data;
    return *this;
}

/**<
 * @brief Multiplication chain operator. Allows for elementwise multiplication of this field with a scalar value.
 * 
 * @param value The scalar value to multiply all field elements by.
 * @return *this A reference to this field after the multiplication.
 */
template <typename T>
Field<T>& Field<T>::operator*=(const T& value) {

    _data = value * _data; // Multiply all field values by the passed value
    return *this;
}

/**
 * @brief Data access operator overload. Allows for read-write access to field values at the specified indices.
 * 
 * @param i Tensor index in the x direction.
 * @param j Tensor index in the y direction.
 * @param k Tensor index in the z direction.
 */
template <typename T>
T& Field<T>::operator()(const int& i, const int& j, const int& k) {

    return _data(i, j, k);

}

/**
 * @brief Checks if the logical coordinate l is within bounds of the field.
 * 
 * @param l The logical coordinate to check (3D vector).
 * @return true if the logical coordinate is within bounds, false otherwise.
 */
template <typename T>
bool Field<T>::in_bound_logical(const Eigen::Vector3d& l) const {

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
template <typename T>
void Field<T>::scatter(const Eigen::Vector3d& l, const T& value) {

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
    _data(i,j,k) += value * (1-di) * (1-dj) * (1-dk);
    _data(i+1,j,k) += value * (di) * (1-dj) * (1-dk);
    _data(i,j+1,k) += value * (1-di) * (dj) * (1-dk); 
    _data(i+1,j+1,k) += value * (di) * (dj) * (1-dk);
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
template <typename T>
T& Field<T>::gather(const Eigen::Vector3d& l, T& value) {

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

    return value;
}

// Field method explicit instantiations
template double& Field<double>::gather(const Eigen::Vector3d& l, double& value);
template void Field<double>::scatter(const Eigen::Vector3d& l, const double& value);

// Field 3 methods

Field3::Field3(const int& ni, const int& nj, const int& nk) : 
    _ni(ni),
    _nj(nj),
    _nk(nk),
    _dataX(ni, nj, nk), 
    _dataY(ni, nj, nk), 
    _dataZ(ni, nj, nk) {

    _dataX.setZero();
    _dataY.setZero();
    _dataZ.setZero();

}

bool Field3::in_bound_logical(const Eigen::Vector3d& l) const {
    
    return (l(0) >= 0 && l(0) <= _ni-1 && l(1) >= 0 && l(1) <= _nj-1 && l(2) >= 0 && l(2) <= _nk-1);
}

void Field3::scatter(const Eigen::Vector3d& l, const Eigen::Vector3d& value) {
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

    // Get the fractional distance to each node
    double di = l(0) - i;
    double dj = l(1) - j;
    double dk = l(2) - k;

    // Scatter data to the surrounding nodes for each component of the field
    _dataX(i,j,k) += value[0] * (1-di) * (1-dj) * (1-dk);
    _dataX(i+1,j,k) += value[0] * (di) * (1-dj) * (1-dk);
    _dataX(i,j+1,k) += value[0] * (1-di) * (dj) * (1-dk); 
    _dataX(i+1,j+1,k) += value[0] * (di) * (dj) * (1-dk);
    _dataX(i,j,k+1) += value[0] * (1-di) * (1-dj) * (dk);
    _dataX(i+1,j,k+1) += value[0] * (di) * (1-dj) * (dk);
    _dataX(i,j+1,k+1) += value[0] * (1-di) * (dj) * (dk);
    _dataX(i+1,j+1,k+1) += value[0] * (di) * (dj) * (dk);

    _dataY(i,j,k) += value[1] * (1-di) * (1-dj) * (1-dk);
    _dataY(i+1,j,k) += value[1] * (di) * (1-dj) * (1-dk);
    _dataY(i,j+1,k) += value[1] * (1-di) * (dj) * (1-dk); 
    _dataY(i+1,j+1,k) += value[1] * (di) * (dj) * (1-dk);
    _dataY(i,j,k+1) += value[1] * (1-di) * (1-dj) * (dk);
    _dataY(i+1,j,k+1) += value[1] * (di) * (1-dj) * (dk);
    _dataY(i,j+1,k+1) += value[1] * (1-di) * (dj) * (dk);
    _dataY(i+1,j+1,k+1) += value[1] * (di) * (dj) * (dk);

    _dataZ(i,j,k) += value[2] * (1-di) * (1-dj) * (1-dk);
    _dataZ(i+1,j,k) += value[2] * (di) * (1-dj) * (1-dk);
    _dataZ(i,j+1,k) += value[2] * (1-di) * (dj) * (1-dk); 
    _dataZ(i+1,j+1,k) += value[2] * (di) * (dj) * (1-dk);
    _dataZ(i,j,k+1) += value[2] * (1-di) * (1-dj) * (dk);
    _dataZ(i+1,j,k+1) += value[2] * (di) * (1-dj) * (dk);
    _dataZ(i,j+1,k+1) += value[2] * (1-di) * (dj) * (dk);
    _dataZ(i+1,j+1,k+1) += value[2] * (di) * (dj) * (dk);

    return;
}


void Field3::gather(const Eigen::Vector3d& l, Eigen::Vector3d& value) {

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
    
    // Interpolate data from surrounding nodes for each component of the field
    value[0] = (1-di)*(1-dj)*(1-dk)*_dataX(i,j,k) +
            (di)*(1-dj)*(1-dk)*_dataX(i+1,j,k) +
            (1-di)*(dj)*(1-dk)*_dataX(i,j+1,k) + 
            (di)*(dj)*(1-dk)*_dataX(i+1,j+1,k) +
            (1-di)*(1-dj)*(dk)*_dataX(i,j,k+1) +
            (di)*(1-dj)*(dk)*_dataX(i+1,j,k+1) +
            (1-di)*(dj)*(dk)*_dataX(i,j+1,k+1) +
            (di)*(dj)*(dk)*_dataX(i+1,j+1,k+1);

    value[1] = (1-di)*(1-dj)*(1-dk)*_dataY(i,j,k) +
            (di)*(1-dj)*(1-dk)*_dataY(i+1,j,k) +
            (1-di)*(dj)*(1-dk)*_dataY(i,j+1,k) + 
            (di)*(dj)*(1-dk)*_dataY(i+1,j+1,k) +
            (1-di)*(1-dj)*(dk)*_dataY(i,j,k+1) +
            (di)*(1-dj)*(dk)*_dataY(i+1,j,k+1) +
            (1-di)*(dj)*(dk)*_dataY(i,j+1,k+1) +
            (di)*(dj)*(dk)*_dataY(i+1,j+1,k+1);

    value[2] = (1-di)*(1-dj)*(1-dk)*_dataZ(i,j,k) +
            (di)*(1-dj)*(1-dk)*_dataZ(i+1,j,k) +
            (1-di)*(dj)*(1-dk)*_dataZ(i,j+1,k) + 
            (di)*(dj)*(1-dk)*_dataZ(i+1,j+1,k) +
            (1-di)*(1-dj)*(dk)*_dataZ(i,j,k+1) +
            (di)*(1-dj)*(dk)*_dataZ(i+1,j,k+1) +
            (1-di)*(dj)*(dk)*_dataZ(i,j+1,k+1) +
            (di)*(dj)*(dk)*_dataZ(i+1,j+1,k+1);

    return; 
}

Eigen::Vector3d& Field3::operator()(const int& i, const int& j, const int& k) {

    static Eigen::Vector3d value; // Static variable to hold the value to return a reference to

    value[0] = _dataX(i,j,k);
    value[1] = _dataY(i,j,k);
    value[2] = _dataZ(i,j,k);

    return value;