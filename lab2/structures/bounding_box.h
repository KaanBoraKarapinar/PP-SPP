#pragma once

#include <cmath>
#include "structures/vector2d.h"
#include <iostream>

class BoundingBox{
public:
    BoundingBox(): x_min(0), x_max(0), y_min(0), y_max(0){};
    BoundingBox(double arg_x_min, double arg_x_max, double arg_y_min, double arg_y_max): x_min(arg_x_min), x_max(arg_x_max), y_min(arg_y_min), y_max(arg_y_max){}

    [[nodiscard]] bool contains(Vector2d<double> position){
        if((x_min <= position[0]) && (position[0] <= x_max) && (y_min <= position[1]) && (position[1] <= y_max)){
            return true;
        }
        return false;
    }

    [[nodiscard]] BoundingBox get_quadrant(std::uint8_t quadrant_id){
        switch (quadrant_id)
        {
        case 0:
            return BoundingBox(x_min, x_min + ((x_max-x_min)/2), y_min+((y_max-y_min)/2), y_max);   
        case 1:
            return BoundingBox(x_min + ((x_max-x_min)/2), x_max, y_min+((y_max-y_min)/2), y_max);
        case 2: 
            return BoundingBox(x_min, x_min + ((x_max-x_min)/2), y_min, y_min + ((y_max-y_min)/2));
        case 3:
            return BoundingBox(x_min + ((x_max-x_min)/2), x_max, y_min, y_min + ((y_max-y_min)/2));
        default:
            throw std::invalid_argument("invalid quadrant id");
        }
    }

    [[nodiscard]] std::string get_string();
    [[nodiscard]] double get_diagonal();
    void plotting_sanity_check();
    [[nodiscard]] BoundingBox get_scaled(std::uint32_t scaling_factor);

    double x_min, x_max, y_min, y_max;
};