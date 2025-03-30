#pragma once

#include <string>
#include "vector2d.h"
#include <cstdint>  // for std::uint32_t

#ifndef BOUNDING_BOX_H
#define BOUNDING_BOX_H

class BoundingBox {

    public:
        // Member-Variables
        double x_min;
        double x_max;
        double y_min;
        double y_max;

        BoundingBox() : x_min(0.0), x_max(0.0), y_min(0.0), y_max(0.0) {} // Default Constructor

        BoundingBox(double x_min_val, double x_max_val, double y_min_val, double y_max_val) : x_min(x_min_val), x_max(x_max_val), y_min(y_min_val), y_max(y_max_val) {};

        std::string get_string();

        double get_diagonal();

        void plotting_sanity_check();

        BoundingBox get_scaled(std::uint32_t scaling_factor);

        bool contains(const Vector2d<double>& position);

        BoundingBox get_quadrant(std::uint8_t index);

};
#endif