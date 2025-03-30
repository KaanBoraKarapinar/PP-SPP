#pragma once

#include "simulation/barnes_hut_simulation.h"

class BarnesHutSimulationWithCollisions : BarnesHutSimulation {
public:
    static void simulate_epochs(Plotter& plotter, Universe& universe, std::uint32_t num_epochs, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs);
    static void simulate_epoch(Plotter& plotter, Universe& universe, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs);

    static void find_collisions(Universe& universe);
    static void find_collisions_parallel(Universe& universe);
};