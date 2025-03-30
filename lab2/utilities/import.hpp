#pragma once
#include <fstream>
#include <iostream>
#include <string>


static void load_universe(std::filesystem::path load_universe_path, Universe& universe){
    // important: ordering of the elements has to be preserverd to ensure correct unpacking!
    std::ifstream universe_file(load_universe_path);
    std::string line;

    if(! universe_file.is_open()){
        throw std::invalid_argument("Could not load universe from given file!");
    }

    // unpack values
    std::int32_t line_count = 0;

    // get body count
    getline(universe_file, line);  // ignore comment line
    getline(universe_file, line);
    line_count = 2;

    std::uint32_t num_bodies = std::stoi(line);
    universe.num_bodies = num_bodies;

    universe.weights.resize(num_bodies);
    universe.velocities.resize(num_bodies);
    universe.positions.resize(num_bodies);
    universe.forces.resize(num_bodies);

    // unpack positions
    getline(universe_file, line); // ignore comment line
    for(int i = 0; i < num_bodies; i++){
        getline(universe_file, line);
        auto delimiter_position = line.find(" ");
        std::string value_1 = line.substr(0, delimiter_position);
        std::string value_2 = line.substr(delimiter_position+1);
        universe.positions[i] = Vector2d(std::stod(value_1), std::stod(value_2));
    }

    // unpack weights
    getline(universe_file, line); // ignore comment line
    for(int i = 0; i < num_bodies; i++){
        getline(universe_file, line);
        universe.weights[i] = std::stod(line);
    }

    // unpack velocities
    getline(universe_file, line); // ignore comment line
    for(int i = 0; i < num_bodies; i++){
        getline(universe_file, line);
        auto delimiter_position = line.find(" ");
        std::string value_1 = line.substr(0, delimiter_position);
        std::string value_2 = line.substr(delimiter_position+1);
        universe.velocities[i] = Vector2d(std::stod(value_1), std::stod(value_2));
    }

    // unpack forces
    getline(universe_file, line); // ignore comment line
    for(int i = 0; i < num_bodies; i++){
        getline(universe_file, line);
        auto delimiter_position = line.find(" ");
        std::string value_1 = line.substr(0, delimiter_position);
        std::string value_2 = line.substr(delimiter_position+1);
        universe.forces[i] = Vector2d(std::stod(value_1), std::stod(value_2));
    }

    

}
