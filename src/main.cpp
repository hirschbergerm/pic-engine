#include <iostream>
#include <string>
#include <Eigen/Dense>

int main(int argc, char* argv[]) {


    Eigen::MatrixXd mat(2,2);
    mat(0,0) = 1;
    mat(0,1) = 2;
    mat(1,0) = 3;
    mat(1,1) = 4;

    std::cout << mat << std::endl;

    std::cout << "Hello, World!" << std::endl;
    return 0;

}