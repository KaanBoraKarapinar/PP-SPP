#pragma once
#include <cstdint> // Included to call specific integer types.
#include <iostream>
#include <vector>
#include <limits>

#include <structures/vector2d.h>
#include <structures/bounding_box.h>

class Universe {

    public:

        // Universe Attributes:
	    std::uint32_t num_bodies; // Total amount of celestial bodies.
	    std::uint32_t current_simulation_epoch; // Current simulation epoch.

        /*
        Standart Constructor:
        Number of bodies: 0 & at 0th Simulation Epoch
        */
        Universe() : num_bodies(0), current_simulation_epoch(0) {
	        std::cout << "A universe with 0 bodies and 0th epoch created by standart constructor" << std::endl;
        }

	    //Vectors containing RESPECTIVE DETAILS OF EVERY body:
	    std::vector<double> weights; // kg
	    std::vector<Vector2d<double>> forces; // N
	    std::vector<Vector2d<double>> velocities; // m/s
	    std::vector<Vector2d<double>> positions; // m

        BoundingBox get_bounding_box() {
            // Empty Universe:
            if (positions.empty()) {
                return BoundingBox(0.0, 0.0, 0.0, 0.0);
            }
            
            // Initial max and min values. These will immediately get updated.
            double x_min = std::numeric_limits<double>::max();
            double x_max = std::numeric_limits<double>::lowest();
            double y_min = std::numeric_limits<double>::max();
            double y_max = std::numeric_limits<double>::lowest();

            // For each body in this universe, the smallest min and largest max values will be found.
            for (const auto& pos : positions) {
                if (pos[0] < x_min) x_min = pos[0];
                if (pos[0] > x_max) x_max = pos[0];
                if (pos[1] < y_min) y_min = pos[1];
                if (pos[1] > y_max) y_max = pos[1];
            }
      
            return BoundingBox(x_min, x_max, y_min, y_max);
            
        }

};