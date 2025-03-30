#pragma once

#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "../structures/universe.h"

// Universe is const, as we do not want to change any universe values
static void save_universe(std::filesystem::path save_universe_path, const Universe& universe) {

	std::ofstream file(save_universe_path);

	if (!file) {
		std::cerr << "save_universe called, file not present" << std::endl;
	}

	file << "### Bodies" << std::endl;
	file << universe.num_bodies << std::endl;

	file << "### Positions" << std::endl;
	//repeat for numbodies
	for (int i1 = 0; i1 < universe.num_bodies; i1++) {
	
		//0 for x coordinate, 1 for y coordinate
		file << universe.positions[i1][0] << " " << universe.positions[i1][1] << std::endl;

	}

	file << "### Weights" << std::endl;
	for (int i2 = 0; i2 < universe.num_bodies; i2++) {

		file << universe.weights[i2] << std::endl;
		
	}

	file << "### Velocities" << std::endl;
	for (int i3 = 0; i3 < universe.num_bodies; i3++) {

		file << universe.velocities[i3][0] << " " << universe.velocities[i3][1] << std::endl;

	}

	file << "### Forces" << std::endl;
	for (int i4 = 0; i4 < universe.num_bodies; i4++) {

		file << universe.forces[i4][0] << " " << universe.forces[i4][1] << std::endl;

	}

}