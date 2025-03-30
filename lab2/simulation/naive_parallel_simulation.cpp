#include "simulation/naive_parallel_simulation.h"
#include "physics/gravitation.h"
#include "physics/mechanics.h"

#include <cmath>

#include "simulation/constants.h"

void NaiveParallelSimulation::simulate_epochs(Plotter& plotter, Universe& universe, std::uint32_t num_epochs, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs){
    for(int i = 0; i < num_epochs; i++){
        simulate_epoch(plotter, universe, create_intermediate_plots, plot_intermediate_epochs);
    }
}

void NaiveParallelSimulation::simulate_epoch(Plotter& plotter, Universe& universe, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs){
    calculate_forces(universe);
    calculate_velocities(universe);
    calculate_positions(universe);
    universe.current_simulation_epoch++;
    if(create_intermediate_plots){
        if(universe.current_simulation_epoch % plot_intermediate_epochs == 0){
            plotter.add_bodies_to_image(universe);
            plotter.write_and_clear();
        }
    }
}



// 3 HELPER METHODS below

/*
    Helper method for calculating distances between two bodies.
    Formula: d = sqrt( (x2 - x1)^2 + (y2 - y1)^2)
*/
double calculate_distance(Universe& universe, std::size_t body1, std::size_t body2) {
    // The displacement vector. (Difference of two vectors.)
    Vector2d<double> disv = universe.positions[body2] - universe.positions[body1];

    double distance = std::sqrt(disv[0] * disv[0] + disv[1] * disv[1]);
    if (distance == 0) throw std::invalid_argument("Error: zero vector/distance in get_normalized_vector at naive sequential simulation");
    return distance;
}

/*
    Helper method for calculating the force vector.
    This method uses force magnitude and distance to calculate the vector.
    The method for calculating gravitational force only gives its magnitude. But universe needs the vector.
    FORMULA: (F * ((x2 - x1) / distance), F * ((y1 - y2) / distance))
*/
Vector2d<double> force_vector(Universe& universe, std::size_t body1, std::size_t body2) {

    // The displacement vector. (Difference of two vectors.)
    Vector2d<double> disv = universe.positions[body2] - universe.positions[body1];

    double distance = calculate_distance(universe, body1, body2);

    double force_magnitude = gravitational_force(universe.weights[body1], universe.weights[body2], distance);

    double x_F = force_magnitude * (disv[0] / distance);
    double y_f = force_magnitude * (disv[1] / distance);

    return Vector2d<double>(x_F, y_f);

}

/*
    Helper method for calculating forces.
    This method calculates the force applied on a single body.
*/

Vector2d<double> calculate_force(Universe& universe, std::size_t body) {
    // Constructor call to initialize force vector.
    Vector2d<double> force(0, 0);

    // Calculate applied forces to body: Part 1
    for (std::size_t j = 0; j < body; j++)
    {
        force = force + force_vector(universe, body, j);
    }

    // Calculate applied forces to body: Part 2
    for (std::size_t i = body + 1; i < universe.num_bodies; i++)
    {
        force = force + force_vector(universe, body, i);
    }
    return force;

}




// IMPLEMENTATION

void NaiveParallelSimulation::calculate_forces(Universe& universe){

    universe.forces.assign(universe.num_bodies, Vector2d<double>(0.0,0.0));

    #pragma omp parallel for
        for (int i = 0; i < universe.num_bodies; ++i ) {
            universe.forces[i] = calculate_force(universe, i);
        }

    return;
}

void NaiveParallelSimulation::calculate_velocities(Universe& universe){

    #pragma omp parallel for
        for (int body_idx = 0; body_idx < universe.num_bodies; ++body_idx) {
            auto acceleration = calculate_acceleration(universe.forces[body_idx], universe.weights[body_idx]);
            universe.velocities[body_idx] = calculate_velocity(universe.velocities[body_idx], acceleration, epoch_in_seconds);
        }

    return;
}

void NaiveParallelSimulation::calculate_positions(Universe& universe){

    #pragma omp parallel for
        for (int i = 0; i < universe.num_bodies;++i) {
            auto s = universe.velocities[i] * epoch_in_seconds;
            universe.positions[i] = universe.positions[i] + s;
        }

    return;
}