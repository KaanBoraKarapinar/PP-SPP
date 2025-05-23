#include "structures/bounding_box.h"
#include <math.h>
#include <stdexcept>

std::string BoundingBox::get_string(){
    std::string res = "" + std::to_string(x_min) + " " + std::to_string(x_max) + "  -  " + std::to_string(y_min) + " " + std::to_string(y_max);
    return res;
}

double BoundingBox::get_diagonal(){
    // NOTE: berechnungsformel auf grund von max(1.0 , ..) vorgeben
    return std::max(1.0, sqrt(pow((x_max-x_min), 2) + pow((y_max-y_min), 2)));
}

void BoundingBox::plotting_sanity_check(){
    // if one direction has a size of 0, convert to a square
    double x_size = x_max - x_min;
    x_size = x_size < 0 ? -x_size : x_size;
    double y_size = y_max - y_min;
    y_size = y_size < 0 ? -y_size : y_size;

    if(x_size < 1e-10 && y_size > 1e-10){
        // fix x size
        //std::cout << "BoundingBox.sanity_check: inflating bounding box in x direction." << std::endl;
        x_min = x_min - (y_size / 2);
        x_max = x_max + (y_size / 2);
        // update x_size
        x_size = x_max - x_min;
        x_size = x_size < 0 ? -x_size : x_size;
    }
    if(y_size < 1e-10 && x_size > 1e-10){
        // fix y size
        //std::cout << "BoundingBox.sanity_check: inflating bounding box in y direction." << std::endl;
        y_min = y_min - (x_size / 2);
        y_max = y_max + (x_size / 2);
    }
    if(x_size < 1e-10 && y_size < 1e-10){
        throw std::invalid_argument("x and y size of bounding box are 0, thus can not be fixed by the sanity check.");
    }
}

BoundingBox BoundingBox::get_scaled(std::uint32_t scaling_factor){
    // calculate edge lenghts and multiply by scaling factor
    double x_middle = (x_min + x_max) / 2;
    double y_middle = (y_min + y_max) / 2;

    double x_size = x_max - x_min;
    double y_size = y_max - y_min;

    return BoundingBox(x_middle - x_size*scaling_factor, x_middle + x_size*scaling_factor, y_middle - y_size*scaling_factor, y_middle + y_size*scaling_factor);
}

bool BoundingBox::contains(const Vector2d<double>& position_vector) {

    double x = position_vector[0];
    double y = position_vector[1];

    if (x <= x_max && x >= x_min && y <= y_max && y >= y_min) {
        return true;
    }
    else
        return false;
}

BoundingBox BoundingBox::get_quadrant(std::uint8_t index) {

    double x_middle = (x_min + x_max) / 2;
    double y_middle = (y_min + y_max) / 2;

    switch (index) {

        case 0: // Top-left Quadrant: 0
            return BoundingBox(x_min, x_middle, y_middle, y_max);
        case 1: // Top-right Quadrant: 1
            return BoundingBox(x_middle, x_max, y_middle, y_max);
        case 2: // Bottom-left Quadrant: 2
            return BoundingBox(x_min, x_middle, y_min, y_middle);
        case 3: // Bottom-right Quadrant: 3
            return BoundingBox(x_middle, x_max, y_min, y_middle);
    
        default:
            throw std::out_of_range("Invalid access!");
    }
}