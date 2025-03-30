#pragma once

#include "common.cuh"

#include "structures/universe.h"
#include "quadtree/quadtree.h"
#include "plotting/plotter.h"


class NaiveCudaSimulation{
public:
    static void allocate_device_memory(Universe& universe, void** d_weights, void** d_forces, void** d_velocities, void** d_positions);
    static void free_device_memory(void** d_weights, void** d_forces, void** d_velocities, void** d_positions);

    static void copy_data_to_device(Universe& universe, void* d_weights, void* d_forces, void* d_velocities, void* d_positions);
    static void copy_data_from_device(Universe& universe, void* d_weights, void* d_forces, void* d_velocities, void* d_positions);

    static void calculate_forces(Universe& universe, void* d_positions, void* d_weights, void* d_forces);
    static void calculate_velocities(Universe& universe, void* d_forces, void* d_weights, void* d_velocities);
    static void calculate_positions(Universe& universe, void* d_velocities, void* d_positions);

    static void simulate_epochs(Plotter& plotter, Universe& universe, std::uint32_t num_epochs, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs);
    static void simulate_epoch(Plotter& plotter, Universe& universe, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs, void* d_weights, void* d_forces, void* d_velocities, void* d_positions);

    static std::vector<std::uint8_t> get_pixels(std::uint32_t plot_width, std::uint32_t plot_height, BoundingBox plot_bounding_box, void* d_positions, std::uint32_t num_bodies);
    static void compress_pixels(std::vector<std::uint8_t>& raw_pixels, std::vector<std::uint8_t>& compressed_pixels);

    // adapter functions used by automatic tests
    static void calculate_forces_kernel_test_adapter(std::uint32_t grid_dim, std::uint32_t block_dim, std::uint32_t num_bodies, void* d_positions, void* d_weights, void* d_forces);
    static void calculate_velocities_kernel_test_adapter(std::uint32_t grid_dim, std::uint32_t block_dim, std::uint32_t num_bodies, void* d_forces, void* d_weights, void* d_velocities);
    static void calculate_positions_kernel_test_adapter(std::uint32_t grid_dim, std::uint32_t block_dim, std::uint32_t num_bodies, void* d_velocities, void* d_positions);
};
