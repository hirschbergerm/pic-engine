#include "output.hpp"
#include <iostream>
#include <fstream>

void Output::fields(World& world) {
    // Build filename
    std::stringstream name;
    name<<"fields.vti";

    // open output file
    std::ofstream file(name.str());
    if (!file.is_open()) {
        std::cerr<<"Error: Could not open file "<<name.str()<<" for writing."<<std::endl;
        return;
    }

    Eigen::Vector3d x_min = world.get_origin();
    Eigen::Vector3d dh = world.get_dh();
    const Field<double>& node_volumes = world.get_node_volumes(); 

    // Header tags
    file<<"VTKFile type=\"ImageData\">/n";
    file<<"ImageData Origin=\""<<x_min(0)<<" "<<x_min(1)<<" "<<x_min(2)<<"\"";
    file<<"Spacing=\""<<dh(0)<<" "<<dh(1)<<" "<<dh(2)<<"\"";
    file<<"WholeExtent=\"0 "<<world._ni-1<<" 0 "<<world._nj-1<<" 0 "<<world._nk-1<<"\">/n";

    // output data stored on nodes 
    file<<"<PointData>\n";

    // Node volumes
    file<<"<DataArray Name=\"NodeVolumes\" type=\"Float64\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    file<<node_volumes;
    file<<"</DataArray>\n";

    // Potential field
    file<<"<DataArray Name=\"Potential\" type=\"Float64\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    file<<world._phi;
    file<<"</DataArray>\n";

    // Charge density field
    file<<"<DataArray Name=\"ChargeDensity\" type=\"Float64\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    file<<world._rho;
    file<<"</DataArray>\n";

    // Electric Field
    file<<"<DataArray Name=\"ElectricField\" type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    file<<world._E;
    file<<"</DataArray>\n";

    // Closing tags
    file<<"</PointData>/n";
    file<<"</ImageData>/n";
    file<<"</VTKFile>/n";

    // Close the file explicitly
    file.close();
}