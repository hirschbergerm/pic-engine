#include "field.hpp"

Field::Field(const int& ni, const int& nj, const int& nk) : _ni(ni), _nj(nj), _nk(nk) {

    _data = Eigen::Tensor<double, 3>(ni, nj, nk);
    _data.setZero(); // Set all of the data to 0


}

Field::~Field() {

    //delete &(_data); Free the allocated memory for data
    
}


/*
* @brief Overload the division operator for Field class
* @param other:
*/
void Field::operator/=(const Field& other) {


}

Field Field::operator+=(const Field& other) {

    Field result(_ni, _nj, _nk);
    result._data = this->_data + other._data;
    return result;

}

