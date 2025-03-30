#pragma once
#include <vector>
#include <iostream>

#include <limits>
#include <filesystem>

#include "structures/vector2d.h"
#include "structures/bounding_box.h"
#include "image/bitmap_image.h"

class Universe {
public:
    Universe(){
        num_bodies = 0;
        current_simulation_epoch = 0;
    }
    void print_bodies_to_console();
    BoundingBox get_bounding_box();
    BoundingBox parallel_cpu_get_bounding_box();


    std::uint32_t num_bodies;
    std::vector<double> weights;  // in kg
    std::vector<Vector2d<double>> forces; // in N
    std::vector<Vector2d<double>> velocities;  // in m/s
    std::vector<Vector2d<double>> positions;  // in m
    std::uint32_t current_simulation_epoch;

};