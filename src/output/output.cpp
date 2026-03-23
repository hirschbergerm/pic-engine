#include "output.hpp"
#include <iostream>
#include <iomanip>

std::ofstream Output::f_diag; // Define the static member variable for diagnostics output file

void Output::fields_output(World& world, std::vector<Species*>& species) {
    // Build filename
    std::stringstream name;
    name<<"/results/fields_"<<std::setfill('0')<<std::setw(5)<<world.get_timestep()<<".vti";

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

    // Species number density
    for (auto& sp: species) {
        file<<"<DataArray Name=\"NumberDensity."<<sp->_name<<"\" type=\"Float64\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        file<<sp->_den;
        file<<"</DataArray>\n";
    }

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

void Output::screen_output(World &world, std::vector<Species*>& species) {
    std::cout<<"time: "<<world.get_time();
    for (auto& sp : species) {
        std::cout<<std::setprecision(3)<<"\t"<<sp->_name<<":"<<0;
    }
    std::cout<<std::endl;
}

void Output::diagnostic_output(World& world, std::vector<Species*>& species) {
    
    if (!f_diag.is_open()){ // If file isn't open we create it
        f_diag.open("runtime_diags.csv"); // Create the file
        f_diag<<"ts,time,wall_time"; // make the csv columns
        for (auto& sp : species) {
            f_diag<<",mp_count."<<sp->_name<<",real_count."<<sp->_name<<",KE."
                <<sp->_name<<",px."<<sp->_name<<",py."<<sp->_name<<",pz."<<sp->_name; // make the csv columns for each species
        }
        f_diag<<",PE,total_E"<<std::endl; // finish header row
    }

    f_diag<<world.get_timestep()<<","<<world.get_time()<<","<<world.get_wall_time(); // Write timestep, sim time, and wall time to file

    double tot_ke = 0.0;
    double ke, px, py, pz, real_count, mp_count;

    for (auto& sp : species) {
        
        sp->get_kinetic_energy(ke);
        sp->get_momentum(px, py, pz);
        sp->get_real_count(real_count);
        sp->get_macro_particle_count(mp_count);
        tot_ke += ke;

        f_diag<<","<<mp_count<<","<<real_count<<","<<ke<<","<<px<<","<<py<<","<<pz; // Write species diagnostics to file
    }

    double pe; 
    world.get_potential_energy(pe);
    f_diag<<","<<pe<<","<<tot_ke + pe; // Write potential energy and total energy to file

    f_diag<<"\n"; // Use \n to avoid flush to disk
    if (world.get_timestep()%25==0) {
        f_diag.flush(); // Flush to disk every 25 timesteps to avoid data loss in case of a crash
    }
    
}