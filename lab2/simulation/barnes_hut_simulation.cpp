#include "simulation/barnes_hut_simulation.h"
#include "simulation/naive_parallel_simulation.h"
#include "physics/gravitation.h"
#include "physics/mechanics.h"
#include "structures/universe.h"
#include "quadtree/quadtreeNode.h"
#include <vector>
#include <iostream>

#include <cmath>

void BarnesHutSimulation::simulate_epochs(Plotter& plotter, Universe& universe, std::uint32_t num_epochs, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs){
    for(int i = 0; i < num_epochs; i++){
        simulate_epoch(plotter, universe, create_intermediate_plots, plot_intermediate_epochs);
    }
}

void BarnesHutSimulation::simulate_epoch(Plotter& plotter, Universe& universe, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs){
    
    BoundingBox bounding_box = universe.get_bounding_box(); // bounding box of universe
    std::int8_t construct_mode; 
    Quadtree* quadtree = nullptr; // quadtree object

    //build quadtree depending on construct_mode
    //serial structure
    if (construct_mode == 0) {
        quadtree = new Quadtree(universe, universe.get_bounding_box(), 0);
    }
    //parallel structure without Cut-Off
    if (construct_mode == 1) {
        quadtree = new Quadtree(universe, universe.get_bounding_box(), 1);
    }
    //parallel structure with Cut-Off
    if (construct_mode == 2) {
        quadtree = new Quadtree(universe, universe.get_bounding_box(), 2);
    }
    else
    {
        throw std::invalid_argument("Invalid construct_mode.");
    }
    
    //calculate center of mass & cumulative mass of quadtree
    quadtree->calculate_center_of_mass();
    quadtree->calculate_cumulative_masses();

    //calculate forces of quadtree in universe
    calculate_forces(universe, *quadtree);

    //calculate velocities and positions
    NaiveParallelSimulation::calculate_velocities(universe);
    NaiveParallelSimulation::calculate_positions(universe);

    //increment epoch counter and plott if needed
    universe.current_simulation_epoch++;
    if (create_intermediate_plots) {
        if ((universe.current_simulation_epoch % plot_intermediate_epochs) == 0) {
            plotter.add_bodies_to_image(universe);
            plotter.write_and_clear();
        }
    }
}

void BarnesHutSimulation::get_relevant_nodes(Universe& universe, Quadtree& quadtree, std::vector<QuadtreeNode*>& relevant_nodes, Vector2d<double>& body_position, std::int32_t body_index, double threshold_theta){
    
   QuadtreeNode* root = quadtree.root;

   get_relevant_nodes_recursively(universe, root, relevant_nodes, body_position, body_index, threshold_theta);
       
}

//Helper method
void BarnesHutSimulation::get_relevant_nodes_recursively(Universe& universe, QuadtreeNode* node, std::vector<QuadtreeNode*>& relevant_nodes, Vector2d<double>& body_position, std::int32_t body_index, double treshold_theta) {

    Vector2d<double> center_of_mass = node->calculate_node_center_of_mass();

    Vector2d<double> relative_vector = (body_position - center_of_mass);// vector in between two points in universe

    double distance = std::sqrt((relative_vector[0] * relative_vector[0]) + (relative_vector[1] * relative_vector[1])); //distance between body K & center of mass of quadrant

    if (distance == 0) {
        throw std::invalid_argument("Distance can't be zero when we need it for calculating theta.");
    }
    
    double diameter = node->bounding_box.get_diagonal(); // diameter of quadrant

    double theta = diameter / distance;
    
    //if quadrant contains body K, then divide it
    if (node->body_identifier == body_index) {
        
        if (node->children.empty()) { // check subquadrants; if true, then no subquadrant
            return; 
        }

        for (QuadtreeNode* child : node->children){ // check child nodes of current node recursively
            get_relevant_nodes_recursively(universe, node, relevant_nodes, body_position, body_index, treshold_theta);
        }
    }

    else {
        
        if (theta < treshold_theta) { // if theta < treshold, then push current node to the relevant nodes
            relevant_nodes.push_back(node);
        }

        else {
            if (!node->children.empty()) { // if current quadrant is irrelevant, but has children, check subquadrants
                for (QuadtreeNode* child : node->children) {
                    get_relevant_nodes_recursively(universe, node, relevant_nodes, body_position, body_index, treshold_theta);
                }
            
            }
        }
    }
}

void BarnesHutSimulation::calculate_forces(Universe& universe, Quadtree& quadtree){

    universe.forces.assign(universe.num_bodies, Vector2d<double>(0.0, 0.0)); // set all forces to zero

    #pragma omp parallel for // parallel calculation of forces for every body 
        for (std::int32_t body_index = 0; body_index < universe.num_bodies; ++body_index) {
            
            std::vector<QuadtreeNode*> relevant_nodes; // vector of relevant quadrants

            Vector2d<double>& body_position = universe.positions[body_index]; // position of current body

            get_relevant_nodes(universe, quadtree, relevant_nodes, body_position, body_index, 0.2); // determine relevant quadrants

            // calculate forces based on relevant quadrants
            
            Vector2d<double> final_force(0.0, 0.0);
            for (QuadtreeNode* node : relevant_nodes ) {
                if (node->body_identifier != 1) { // leaf node, calculate forces directly
                    Vector2d<double> relative_vector = universe.positions[node->body_identifier] - body_position;
                    double distance = std::sqrt(relative_vector[0] * relative_vector[0] + relative_vector[1] * relative_vector[1]);

                    if (distance != 0) {
                        double force_size = gravitational_force(universe.weights[body_index], universe.weights[node->body_identifier], distance);
                        Vector2d<double> scaled_vector = relative_vector * (force_size / distance);
                        final_force = final_force +  scaled_vector;
                    }
                }

                else {// inner node; use center of mass & cumulative mass
                    Vector2d<double> relative_vector = node->center_of_mass - body_position;
                    double distance = std::sqrt(relative_vector[0] * relative_vector[0] + relative_vector[1] * relative_vector[1]);

                    if (distance != 0) {
                        double force_size = gravitational_force(universe.weights[body_index], node->cumulative_mass, distance);
                        Vector2d<double> scaled_vector = relative_vector * (force_size / distance);
                        final_force = final_force + scaled_vector; // normalized direction vector
                    }
                
                }
            
            } // update total forces of body
            universe.forces[body_index] = final_force;
        }

    return;
}