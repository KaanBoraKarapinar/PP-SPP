#include "simulation/naive_sequential_simulation.h"
#include "simulation/constants.h"
#include "physics/gravitation.h"
#include "physics/mechanics.h"

#include <cmath>



void NaiveSequentialSimulation::simulate_epochs(Plotter& plotter, Universe& universe, std::uint32_t num_epochs, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs){
    for(int i = 0; i < num_epochs; i++){
        simulate_epoch(plotter, universe, create_intermediate_plots, plot_intermediate_epochs);
    }
}

void NaiveSequentialSimulation::simulate_epoch(Plotter& plotter, Universe& universe, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs){
    calculate_forces(universe);
    calculate_velocities(universe);
    calculate_positions(universe);
    universe.current_simulation_epoch++;
    if(create_intermediate_plots){
        if((universe.current_simulation_epoch % plot_intermediate_epochs) == 0){
            plotter.add_bodies_to_image(universe);
            plotter.write_and_clear();
        }
    }
}


void NaiveSequentialSimulation::calculate_forces(Universe& universe){
    for(int body_idx = 0; body_idx < universe.num_bodies; body_idx++){
        // get body positions
        Vector2d<double> body_position = universe.positions[body_idx];

        // get body mass
        double body_mass = universe.weights[body_idx];

        Vector2d<double> applied_force_vector;

        for(int distant_body_idx = 0; distant_body_idx < universe.num_bodies; distant_body_idx++){
            if(body_idx == distant_body_idx){
                continue;
            }
            // get distant body positions
            Vector2d<double> distant_body_position = universe.positions[distant_body_idx];

            // calculate vector between bodies to get the direction of the gravitational force
            Vector2d<double> direction_vector = distant_body_position - body_position;

            // calculate the distance between the bodies
            double distance = sqrt(pow(direction_vector[0], 2) + pow(direction_vector[1], 2));

            // calculate gravitational force between the bodies
            double force = gravitational_force(body_mass, universe.weights[distant_body_idx], distance);

            // create the force vector
            Vector2d<double> force_vector = direction_vector * (force / distance);

            // sum forces applied to body
            applied_force_vector = applied_force_vector + force_vector;
         }

        // store applied force 
        universe.forces[body_idx] = applied_force_vector;
    }
}


void NaiveSequentialSimulation::calculate_velocities(Universe& universe){
    // calculate velocity due to applied force
    for(int body_idx = 0; body_idx < universe.num_bodies; body_idx++){
        auto acceleration = calculate_acceleration(universe.forces[body_idx], universe.weights[body_idx]);
        universe.velocities[body_idx] = calculate_velocity(universe.velocities[body_idx], acceleration, epoch_in_seconds);
    }
}


void NaiveSequentialSimulation::calculate_positions(Universe& universe){
    for(int body_idx = 0; body_idx < universe.num_bodies; body_idx++){
        // calculate movement
        // s = v * t
        Vector2d<double> movement = universe.velocities[body_idx] * epoch_in_seconds;

        // calculate new position
        // p` = p0 + s 
        Vector2d<double> new_position = universe.positions[body_idx] + movement;

        // update position
        universe.positions[body_idx] = new_position;
    }
}
