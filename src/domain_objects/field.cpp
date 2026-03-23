#include "field.hpp"
#include <algorithm>

template <typename T>
Field<T>::Field(size_t ni, size_t nj, size_t nk) : _ni(ni), _nj(nj), _nk(nk) {

    _data.resize(ni * nj * nk); // Resize the data vector to hold all field values
    _data.assign(_data.size(), T{}); // Initialize all field values to 0.0

}

/** 
* @brief Move constructor for the Field class. Moves the data from the other field to this field.
* 
* @param other The other field to move from.
*/
template <typename T>
Field<T>::Field(Field<T>&& other) noexcept : 
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

    _data.assign(_data.size(), value); // Set all of the Eigen tensor values to value

}

/**
 * @brief Multiplication chain operator. Allows for elementwise multiplication of this field with a scalar value.
 * 
 * @param value The scalar value to multiply all field elements by.
 * @return *this A reference to this field after the multiplication.
 */
template <typename T>
Field<T>& Field<T>::operator*=(const T& value) {

    std::transform(_data.begin(), _data.end(), _data.begin(), [&value](T& element) {
        return element * value; // Multiply each element by the passed value
    });

    return *this;
}

/**
 * @brief Divides all field values by the passed value.
 * 
 * @param value The value to divide all field elements by.
 * @return void
 */
template <typename T>
Field<T>& Field<T>::operator/=(const T& value) {

    assert(value != 0); // Make sure we are not dividing by zero

    // Modify in place using lambda function
    std::transform(_data.begin(), _data.end(), _data.begin(), [&value](T& element) {
        return element / value; // Divide each element by the passed value
    });

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

    std::transform(_data.begin(), _data.end(), other._data.begin(), _data.begin(), std::plus<T>());
    return *this;
}

/**
 * @brief Subtraction chain operator. Allows for elementwise subtraction of this field with another field.
 * 
 * @param other The other field to subtract from this field.
 * @return *this A reference to this field after the subtraction.
 */
template <typename T>
Field<T>& Field<T>::operator-=(const Field<T>& other) {

    std::transform(_data.begin(), _data.end(), other._data.begin(), _data.begin(), std::minus<T>());
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

    // Make sure that the dimensions of the two fields match up, otherwise something has gone really wrong
    assert(_data.size() == other._data.size());

    std::transform(_data.begin(), _data.end(), other._data.begin(), _data.begin(), std::divides<T>());
    return *this;

}

/**
 * @brief Output stream operator for Field class. 
 */
template <typename T>
std::ostream& operator<<(std::ostream& out, const Field<T>& field) {
    
    for (size_t k=0; k < field.nk(); k++,out<<"\n") { // Write a new line each time we move to a new k plane
        for (size_t j=0; j < field.nj(); j++) {
            for (size_t i=0; i< field.ni(); i++) {
                out<<field(i,j,k)<<" "; // Write the field value at (i,j,k) followed by a space
            }
        }
    }
    return out;
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
        // TODO: Throw exception here and handle it.
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
    _data[index(i,j,k)] += value * (1-di) * (1-dj) * (1-dk);
    _data[index(i+1,j,k)] += value * (di) * (1-dj) * (1-dk);
    _data[index(i,j+1,k)] += value * (1-di) * (dj) * (1-dk); 
    _data[index(i+1,j+1,k)] += value * (di) * (dj) * (1-dk);
    _data[index(i,j,k+1)] += value * (1-di) * (1-dj) * (dk);
    _data[index(i+1,j,k+1)] += value * (di) * (1-dj) * (dk);
    _data[index(i,j+1,k+1)] += value * (1-di) * (dj) * (dk);
    _data[index(i+1,j+1,k+1)] += value * (di) * (dj) * (dk);

}

/*
 * @brief Gathers a field value at logical coordinate (l) from the adjacent cell nodes. Modifies the passed value.
 * 
 * @param l The logical coordinate to gather from (3D vector).
 * @param value The field value to gather into.
 * @return void
 */
template <typename T>
void Field<T>::gather(const Eigen::Vector3d& l, T& value) {

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
    value = (1-di) * (1-dj) * (1-dk) * _data[index(i,j,k)] +
            (di) * (1-dj) * (1-dk) * _data[index(i+1,j,k)] +
            (1-di) * (dj) * (1-dk) * _data[index(i,j+1,k)] + 
            (di) * (dj) * (1-dk) * _data[index(i+1,j+1,k)] +
            (1-di) * (1-dj) * (dk) * _data[index(i,j,k+1)] +
            (di) * (1-dj) * (dk) * _data[index(i+1,j,k+1)] +
            (1-di) * (dj) * (dk) * _data[index(i,j+1,k+1)] +
            (di) * (dj) * (dk) * _data[index(i+1,j+1,k+1)];

}

// Explicit instantiation of the stream operator for Field<double>
template std::ostream& operator<< <double>(std::ostream& out, const Field<double>& field);

// Field 3 methods
Field3::Field3(size_t ni, size_t nj, size_t nk) : 
    _ni(ni),
    _nj(nj),
    _nk(nk) {

    // Resize the data vector to hold all field values
    _dataX.resize(ni * nj * nk); 
    _dataY.resize(ni * nj * nk);
    _dataZ.resize(ni * nj * nk);

    // Initialize all field values to 0.0
    _dataX.assign(_dataX.size(),0.0); 
    _dataY.assign(_dataY.size(), 0.0);
    _dataZ.assign(_dataZ.size(), 0.0);

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
    x(i,j,k) += value[0] * (1-di) * (1-dj) * (1-dk);
    x(i+1,j,k) += value[0] * (di) * (1-dj) * (1-dk);
    x(i,j+1,k) += value[0] * (1-di) * (dj) * (1-dk); 
    x(i+1,j+1,k) += value[0] * (di) * (dj) * (1-dk);
    x(i,j,k+1) += value[0] * (1-di) * (1-dj) * (dk);
    x(i+1,j,k+1) += value[0] * (di) * (1-dj) * (dk);
    x(i,j+1,k+1) += value[0] * (1-di) * (dj) * (dk);
    x(i+1,j+1,k+1) += value[0] * (di) * (dj) * (dk);

    y(i,j,k) += value[1] * (1-di) * (1-dj) * (1-dk);
    y(i+1,j,k) += value[1] * (di) * (1-dj) * (1-dk);
    y(i,j+1,k) += value[1] * (1-di) * (dj) * (1-dk); 
    y(i+1,j+1,k) += value[1] * (di) * (dj) * (1-dk);
    y(i,j,k+1) += value[1] * (1-di) * (1-dj) * (dk);
    y(i+1,j,k+1) += value[1] * (di) * (1-dj) * (dk);
    y(i,j+1,k+1) += value[1] * (1-di) * (dj) * (dk);
    y(i+1,j+1,k+1) += value[1] * (di) * (dj) * (dk);

    z(i,j,k) += value[2] * (1-di) * (1-dj) * (1-dk);
    z(i+1,j,k) += value[2] * (di) * (1-dj) * (1-dk);
    z(i,j+1,k) += value[2] * (1-di) * (dj) * (1-dk); 
    z(i+1,j+1,k) += value[2] * (di) * (dj) * (1-dk);
    z(i,j,k+1) += value[2] * (1-di) * (1-dj) * (dk);
    z(i+1,j,k+1) += value[2] * (di) * (1-dj) * (dk);
    z(i,j+1,k+1) += value[2] * (1-di) * (dj) * (dk);
    z(i+1,j+1,k+1) += value[2] * (di) * (dj) * (dk);

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
    value[0] = (1-di) * (1-dj) * (1-dk) * x(i,j,k) +
            (di) * (1-dj) * (1-dk) * x(i+1,j,k) +
            (1-di) * (dj) * (1-dk) * x(i,j+1,k) + 
            (di) * (dj) * (1-dk) * x(i+1,j+1,k) +
            (1-di) * (1-dj) * (dk) * x(i,j,k+1) +
            (di) * (1-dj) * (dk) * x(i+1,j,k+1) +
            (1-di) * (dj) * (dk) * x(i,j+1,k+1) +
            (di) * (dj) * (dk) * x(i+1,j+1,k+1);

    value[1] = (1-di) * (1-dj) * (1-dk) * y(i,j,k) +
            (di) * (1-dj) * (1-dk) * y(i+1,j,k) +
            (1-di) * (dj) * (1-dk) * y(i,j+1,k) + 
            (di) * (dj) * (1-dk) * y(i+1,j+1,k) +
            (1-di) * (1-dj) * (dk) * y(i,j,k+1) +
            (di) * (1-dj) * (dk) * y(i+1,j,k+1) +
            (1-di) * (dj) * (dk) * y(i,j+1,k+1) +
            (di) * (dj) * (dk) * y(i+1,j+1,k+1);

    value[2] = (1-di) * (1-dj) * (1-dk) * z(i,j,k) +
            (di) * (1-dj) * (1-dk) * z(i+1,j,k) +
            (1-di) * (dj) * (1-dk) * z(i,j+1,k) + 
            (di) * (dj) * (1-dk) * z(i+1,j+1,k) +
            (1-di) * (1-dj) * (dk) * z(i,j,k+1) +
            (di) * (1-dj) * (dk) * z(i+1,j,k+1) +
            (1-di) * (dj) * (dk) * z(i,j+1,k+1) +
            (di) * (dj) * (dk) * z(i+1,j+1,k+1);

    return; 
}

std::ostream& operator<<(std::ostream& out, const Field3& field) {

    Eigen::Vector3d vec({0,0,0});
    for (size_t k=0; k < field.nk(); k++,out<<"\n") { // Write a new line each time we move to a new k plane
        for (size_t j=0; j < field.nj(); j++) {
            for (size_t i=0; i< field.ni(); i++) {
                vec = field(i,j,k);
                out<<vec[0]<<" "<<vec[1]<<" "<<vec[2]<<" ";
            }
        }
    }
    return out;

}
