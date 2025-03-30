#pragma once

#include <iostream>
#include <fstream>

static void save_universe(std::filesystem::path file_path, Universe& universe){
    // std::cout << "Saving universe to: " << file_path << std::endl;
    std::ofstream universe_file;
    universe_file.open(file_path);
    
    // store settings
    universe_file << "### Bodies" << std::endl;
    universe_file << universe.num_bodies << std::endl;

    // store positions
    universe_file << "### Positions" << std::endl;
    for(auto position: universe.positions){
        universe_file << std::to_string(position[0]) << " " << std::to_string(position[1]) << std::endl;
    }

    // store weights
    universe_file << "### Weights" << std::endl;
    for(auto weight: universe.weights){
        universe_file << std::to_string(weight) << std::endl;
    }

    // store velocities
    universe_file << "### Velocities" << std::endl;
    for(auto velocity: universe.velocities){
        universe_file << std::to_string(velocity[0]) << " " << std::to_string(velocity[1]) << std::endl;
    }

    // store forces
    universe_file << "### Forces" << std::endl;
    for(auto force: universe.forces){
        universe_file << std::to_string(force[0]) << " " << std::to_string(force[1]) << std::endl;
    }

    universe_file.close();
}


