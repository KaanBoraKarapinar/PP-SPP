#pragma once
#include "structures/bounding_box.h"
#include "image/bitmap_image.h"
#include "image/pixel.h"
#include "structures/universe.h"
#include <cstdint>
#include <set>
#include <filesystem>
#include "structures/vector2d.h"

// VORLAGE START:
class Plotter{
public:
    Plotter(BoundingBox bb, const std::filesystem::path & arg_output_folder_path, std::uint32_t plot_width_arg, std::uint32_t plot_height_arg): plot_bounding_box(bb), output_folder_path(arg_output_folder_path), plot_width(plot_width_arg), plot_height(plot_height_arg), image(BitmapImage(plot_height_arg, plot_width_arg)), image_serial_number(0){
        // set default filename prefix
        filename_prefix = "plot";
    }

    void add_bodies_to_image(Universe& universe);
    void highlight_position(Vector2d<double> position, std::uint8_t red, std::uint8_t green, std::uint8_t blue);
    
    void set_plot_bounding_box(BoundingBox bb){
        plot_bounding_box = bb;
    }

    void write_and_clear();
    
    void clear_image(){
        image = BitmapImage(plot_height, plot_width);
    }

    void set_filename_prefix(std::string prefix){
        filename_prefix = prefix;
    }

    void mark_position(Vector2d<double> position, std::uint8_t red, std::uint8_t green, std::uint8_t blue);
    void mark_pixel(std::uint32_t x, std::uint32_t y, std::uint8_t red, std::uint8_t green, std::uint8_t blue);

    BitmapImage::BitmapPixel get_pixel(std::uint32_t x, std::uint32_t y);

    std::uint32_t get_next_image_serial_number(){
        return image_serial_number;
    }

private:
    std::string filename_prefix;
    std::uint32_t image_serial_number;
    BitmapImage image;
    BoundingBox plot_bounding_box;
    std::uint32_t plot_width, plot_height;
    std::filesystem::path output_folder_path;
};
// VORLAGE END