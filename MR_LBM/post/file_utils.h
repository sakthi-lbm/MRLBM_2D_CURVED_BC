#ifndef FILE_UTILS_H
#define FILE_UTILS_H

#include<string>
#include<fstream>
#include<sstream>
#include<iostream>
#include<string>
#include <vector>


struct InputParameters
{
    double Re;
    double ncy;
    double Dcy;
    double uo;
    double rho_infty;
};


[[nodiscard]] inline int count_lines_in_file(const std::string& filename){

    std::ifstream count_file(filename);
    if(!count_file.is_open()){
        std::cerr << "Failed to open the file" << std::endl;
        return 1;
    }

    int line_count = 0;
    std::string dummy;
    while(std::getline(count_file, dummy)){
        line_count++;
    }
    count_file.close();

   return (line_count);
}

[[nodiscard]] inline int get_npoints_on_cylinder(const std::string& filename){
    std::ifstream theta_file(filename);
    if(!theta_file.is_open()){
        std::cout << "Error opening pressure file!!" <<std::endl;
        return 1;
    }
    
    int npoints;
    std::string line;
    std::getline(theta_file, line);
    std::istringstream iss(line);
    iss >> npoints;

    theta_file.close();


    return(npoints);
}

// Function to read input from file
[[nodiscard]] inline InputParameters read_input_info(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Cannot open input file.");
    }

    std::string line;
    InputParameters params;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string key;
        double value;

        // Match lines like "Dcy = 32.9848440"
        if (line.find("Re") != std::string::npos && iss >> key >> key >> value) {
            params.Re = value;
        } else if (line.find("uo") != std::string::npos && iss >> key >> key >> value) {
            params.uo = value;
        } else if (line.find("ncy") != std::string::npos && iss >> key >> key >> value) {
            params.ncy = value;
        } else if (line.find("Dcy") != std::string::npos && iss >> key >> key >> value) {
            params.Dcy = value;
        } else if (line.find("rho_infty") != std::string::npos && iss >> key >> key >> value) {
            params.rho_infty = value;
        }
    }

    return params;
}

[[nodiscard]] inline int read_theta(const std::string& filename, const int npoints, std::vector<double>& theta) {
    std::ifstream theta_file(filename);
    if(!theta_file.is_open()){
        std::cout << "Error opening pressure file!!" <<std::endl;
        return 1;
    }
    
    int dummy;
    std::string line;
    std::getline(theta_file, line);
    std::istringstream iss(line);

    iss >> dummy;  // read the first integer (npoints) and ignore
    for (int j = 0; j < npoints; j++) {
        if (!(iss >> theta[j])) {
            std::cout << "Error reading theta value at index " << j << std::endl;
            return 1;
        }
    }

    return 0;

}

[[nodiscard]] inline int read_forces(const std::string& filename, const int nsteps, std::vector<int>& time, std::vector<double>& F_drag, std::vector<double>& F_lift){
    std::ifstream force_file(filename);
    if(!force_file.is_open()){
        std::cout << "Failed to open the force file!" << std::endl;
        return 1;
    }

    int step;
    double fx,fy;

    for(int i=0; i<nsteps; i++){
        force_file >> step >> fx >> fy;
        time[i] = step;
        F_drag[i] = fx;
        F_lift[i] = fy;
    }
    force_file.close();


    return 0;
}

[[nodiscard]] inline int read_pressure(const std::string& filename, const int nsteps, const int npoints, std::vector<std::vector<double>>& ps){
    std::ifstream pressure_file(filename);
   
    if(!pressure_file.is_open()){
        std::cout << "Error opening pressure file!!" <<std::endl;
        return 1;
    }

    int i = 0;
    std::string line;
    while (std::getline(pressure_file, line) && i < nsteps) {
        std::istringstream iss(line);
        int dummy;
        iss >> dummy;
        
        for (int j = 0; j < npoints; j++) {
            iss >> ps[i][j];
        }
        i++;
    }
    pressure_file.close();

    return 0;
}

[[nodiscard]] inline int read_pressure_binary(const std::string& filename, const int nsteps, const int npoints, std::vector<std::vector<double>>& ps) {
    std::ifstream pressure_file(filename, std::ios::binary);
    if (!pressure_file.is_open()) {
        std::cout << "Error opening pressure file!!" << std::endl;
        return 1;
    }

    for (int i = 0; i < nsteps; ++i) {

        double step_dummy;
        pressure_file.read(reinterpret_cast<char*>(&step_dummy), sizeof(double));

        pressure_file.read(reinterpret_cast<char*>(ps[i].data()), npoints * sizeof(double));
        if (!pressure_file) {
            std::cout << "Error reading at step " << i << std::endl;
            return 1;
        }
    }

    pressure_file.close();
    return 0;
}

[[nodiscard]] inline int read_cp(const std::string& filename, double *cps, double *cpb){

    std::ifstream count_file(filename);
    if (!count_file.is_open()) {
        std::cout << "Error opening cp file!!" << std::endl;
        return 1;
    }
    std::string line;
    double theta_first, pressure_first;
    double theta_last, pressure_last;

    // Read first line
    if (std::getline(count_file, line)) {
        std::istringstream iss(line);
        iss >> theta_first >> *cpb;
    }

    // Read last line
    while (std::getline(count_file, line)) {
        std::istringstream iss(line);
        iss >> theta_last >> *cps;
    }
    count_file.close();
    return 0;
}





















#endif  //FILE_UTILS_H