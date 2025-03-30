#include "plotting/plotter.h"
#include "io/image_parser.h"

#include <exception>
#include <cmath>


void Plotter::write_and_clear(){
    // create plot serial number string
    std::string serial_number_string = std::to_string(image_serial_number);
    while(serial_number_string.length() < 9){
        serial_number_string = "0" + serial_number_string;
    }

    std::string file_name = filename_prefix + "_" + serial_number_string + ".bmp";
    ImageParser::write_bitmap(output_folder_path / file_name, image);
    clear_image();
    image_serial_number += 1;
}

BitmapImage::BitmapPixel Plotter::get_pixel(std::uint32_t x, std::uint32_t y){
    return image.get_pixel(y, x);
}

void Plotter::mark_pixel(std::uint32_t x, std::uint32_t y, std::uint8_t red, std::uint8_t green, std::uint8_t blue) {
    
    if (x >= plot_width || y >= plot_height) {
        throw std::out_of_range("Pixel is out of range!");
    }

    BitmapImage::BitmapPixel color = {red, green, blue};

    image.set_pixel(y, x, color);
   
}

void Plotter::mark_position(Vector2d<double> position, std::uint8_t red, std::uint8_t green, std::uint8_t blue) {

    bool out_of_bounding_box = 
        position[0] < plot_bounding_box.x_min || 
        position[0] > plot_bounding_box.x_max ||
        position[1] < plot_bounding_box.y_min || 
        position[1] > plot_bounding_box.y_max;
    
    if (out_of_bounding_box) {
        std::cout << "Position is out of bounds!" << std::endl;
    }
    else
    {   
        // Values of The Pixel To Be Painted
        std::uint32_t pixel_x;
        std::uint32_t pixel_y;

        // CHECK EDGE CASES:
        if (position[0] == plot_bounding_box.x_max) pixel_x = plot_width - 1;
        //else if (position[0] == plot_bounding_box.x_min) pixel_x = 0;
        else
        {
            // Normalize position to [0, 1] range
            //double normalized_x = (position[0] - plot_bounding_box.x_min) / (plot_bounding_box.x_max - plot_bounding_box.x_min);

            // Map to pixel indices
            //pixel_x = static_cast<u_int32_t>(std::floor(normalized_x * plot_width));
            
            pixel_x = static_cast<std::uint32_t>((position[0] - plot_bounding_box.x_min) / 
            (plot_bounding_box.x_max - plot_bounding_box.x_min) * (plot_width-1));
        }
        if (position[1] == plot_bounding_box.y_max) pixel_y = plot_height - 1;
        //else if (position[1] == plot_bounding_box.y_min) pixel_y = 0;
        else
        {
            //double normalized_y = (position[1] - plot_bounding_box.y_min) / (plot_bounding_box.y_max - plot_bounding_box.y_min);
            //pixel_y = static_cast<uint32_t>(std::floor(normalized_y * plot_height));
            pixel_y = static_cast<std::uint32_t>((position[1] - plot_bounding_box.y_min) / 
            (plot_bounding_box.y_max - plot_bounding_box.y_min) * (plot_height-1));
        }

        mark_pixel(pixel_x, pixel_y, red, green, blue);
    }   
}


void Plotter::highlight_position(Vector2d<double> position, std::uint8_t red, std::uint8_t green, std::uint8_t blue) {
    
    bool out_of_bounding_box = 
        position[0] < plot_bounding_box.x_min || 
        position[0] > plot_bounding_box.x_max ||
        position[1] < plot_bounding_box.y_min || 
        position[1] > plot_bounding_box.y_max;
    
    if (!out_of_bounding_box)
    {
        for (uint32_t i = 0; i < plot_height; ++i)
        {
            Vector2d<double> vertical(position[0], i);
            mark_position(vertical, red, green, blue);
        }

        for (uint32_t i = 0; i < plot_width; ++i)
        {
            Vector2d<double> horizontal(i, position[1]);
            mark_position(horizontal, red, green, blue);
        }
        //Vector2d<double> up_cross(position[0], position[1]+1);
        ///Vector2d<double> down_cross(position[0], position[1]-1);
        //Vector2d<double> right_cross(position[0]+1, position[1]);
        //Vector2d<double> left_cross(position[0]-1, position[1]);

        //mark_position(position, red, green, blue);
        //mark_position(up_cross, red, green, blue);
        //mark_position(down_cross, red, green, blue);
        //mark_position(right_cross, red, green, blue);
        //mark_position(left_cross, red, green, blue);
    }
}

void Plotter::add_bodies_to_image(Universe& universe) {

    for (std::size_t i = 0; i < universe.num_bodies; ++i)
    {
        bool out_of_bounding_box = 
        universe.positions[i][0] < plot_bounding_box.x_min || 
        universe.positions[i][0] > plot_bounding_box.x_max ||
        universe.positions[i][1] < plot_bounding_box.y_min || 
        universe.positions[i][1] > plot_bounding_box.y_max;

        if (!out_of_bounding_box)
        {
            mark_position(universe.positions[i], 255,255,255);
        }
    }
}