#pragma once


#include "structures/universe.h"
#include "quadtree/quadtree.h"
#include "plotting/plotter.h"

class BarnesHutSimulation{
public:
    static void simulate_epochs(Plotter& plotter, Universe& universe, std::uint32_t num_epochs, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs);
    static void simulate_epoch(Plotter& plotter, Universe& universe, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs);
    static void calculate_forces(Universe& universe, Quadtree& quadtree);
    static void get_relevant_nodes(Universe& universe, Quadtree& quadtree, std::vector<QuadtreeNode*>& relevant_nodes, Vector2d<double>& body_position, std::int32_t body_index, double threshold_theta);
    static void get_relevant_nodes_recursively(Universe& universe, QuadtreeNode* node, std::vector<QuadtreeNode*>& relevant_nodes, Vector2d<double>& body_position, std::int32_t body_index, double treshold_theta);
};