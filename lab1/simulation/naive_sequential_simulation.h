#pragma once

#include "structures/universe.h"
#include "plotting/plotter.h"

class NaiveSequentialSimulation{
public:
    static void simulate_epochs(Plotter& plotter, Universe& universe, std::uint32_t num_epochs, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs);
    static void simulate_epoch(Plotter& plotter, Universe& universe, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs);
    static void calculate_forces(Universe& universe);
    static void calculate_velocities(Universe& universe);
    static void calculate_positions(Universe& universe);

};