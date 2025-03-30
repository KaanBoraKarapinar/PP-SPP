#include "naive_cuda_simulation.cuh"
#include "physics/gravitation.h"
#include "physics/mechanics.h"
#include "simulation/constants.h"
#include <cuda.h>
#include <cuda_runtime_api.h>
#include "cuda_wrappers.cuh"


// DEVICE FRIENDLY HILFSMETHODEN:


__constant__ double gravitational_constant_device = 6.67430e-11; // (m^3)/(kg*s^2)


__host__ __device__ inline
double compute_gravitational_force_device(double mass_1, double mass_2, double distance) {
    return gravitational_constant_device * ((mass_1 * mass_2) / (distance * distance));
}


constexpr double epoch_time_in_seconds_device = 2.628e+6; 




void NaiveCudaSimulation::allocate_device_memory(Universe& universe, void** d_weights, void** d_forces, void** d_velocities, void** d_positions){

 size_t num_bodies = universe.num_bodies;

//-> num bodies times double
cudaMalloc(d_weights, num_bodies * sizeof(double));

//-> VECTOR 2D IN PREV MODULES IS NOW DOUBLE2
    cudaMalloc(d_forces, num_bodies * sizeof(double2));
    cudaMalloc(d_velocities, num_bodies * sizeof(double2));
    cudaMalloc(d_positions, num_bodies * sizeof(double2));


}

void NaiveCudaSimulation::free_device_memory(void** d_weights, void** d_forces, void** d_velocities, void** d_positions){

cudaFree(*d_weights);
cudaFree(*d_forces);
cudaFree(*d_velocities);
cudaFree(*d_positions);

*d_weights = nullptr;
    *d_forces = nullptr;
    *d_velocities = nullptr;
    *d_positions = nullptr;


}


/*
implementation decision
-> no explicit performance criteria but instead of copying entire vectors, recieving constant references
better
-> check later if any problems

*/
std::vector<double2> vector2d_to_double2_translator(const std::vector<Vector2d<double>>& vector2d){

    size_t vector_size = vector2d.size(); //ideally num bodies if debug needed
    std::vector<double2> vector_to_return(vector_size);

    for (size_t i = 0; i < vector_size; i++) {
        //0 = x 1= y
        vector_to_return[i] = make_double2(vector2d[i][0], vector2d[i][1]);
    }

    return vector_to_return;
}


void NaiveCudaSimulation::copy_data_to_device(Universe& universe, void* d_weights, void* d_forces, void* d_velocities, void* d_positions){

/*

-> weight in plain double type

// DOES IT WORK???
     expected = an pointer stating the start of the array, like c type
     .data() returns a pointer for c++ vector types
     another implementation (doesnt match with verbindliche anforderung but implement if error)
    for loop for cudamemcpy or c type of arrays


*/
size_t num_bodies = universe.num_bodies;



//weights are in plain dobule vector type, no translation needed, loading directly
cudaMemcpy(d_weights, universe.weights.data(), num_bodies * sizeof(double), cudaMemcpyHostToDevice);

 std::vector<double2> forces = vector2d_to_double2_translator(universe.forces);
    std::vector<double2> velocities = vector2d_to_double2_translator(universe.velocities);
    std::vector<double2> positions = vector2d_to_double2_translator(universe.positions);



cudaMemcpy(d_forces, forces.data(), num_bodies * sizeof(double2), cudaMemcpyHostToDevice);
    cudaMemcpy(d_velocities, velocities.data(), num_bodies * sizeof(double2), cudaMemcpyHostToDevice);
    cudaMemcpy(d_positions, positions.data(), num_bodies * sizeof(double2), cudaMemcpyHostToDevice);

}

std::vector<Vector2d<double>> double2_to_vector2d_translator(const std::vector<double2>& cuda_vector){

    size_t vector_size = cuda_vector.size(); //ideally num bodies if debug needed
    std::vector<Vector2d<double>> vector_to_return(vector_size);

    for (size_t i = 0; i < vector_size; ++i) {
        //0 = x 1= y
        vector_to_return[i] = Vector2d{cuda_vector[i].x,cuda_vector[i].y};
    }

    return vector_to_return;
}


/*
implementation decision/gedankengang:
    -> saving values into universe object, as no explicit host variables are given.


*/
void NaiveCudaSimulation::copy_data_from_device(Universe& universe, void* d_weights, void* d_forces, void* d_velocities, void* d_positions){
 
     size_t num_bodies = universe.num_bodies;

     cudaMemcpy(universe.weights.data(), d_weights, num_bodies * sizeof(double), cudaMemcpyDeviceToHost);


    //pre translation data types
 std::vector<double2> forces(num_bodies);
    std::vector<double2> velocities(num_bodies);
    std::vector<double2> positions(num_bodies);

 
 
 cudaMemcpy(forces.data(), d_forces, num_bodies * sizeof(double2), cudaMemcpyDeviceToHost);
    cudaMemcpy(velocities.data(), d_velocities, num_bodies * sizeof(double2), cudaMemcpyDeviceToHost);
    cudaMemcpy(positions.data(), d_positions, num_bodies * sizeof(double2), cudaMemcpyDeviceToHost);


//setting universe variables to translated data types

universe.forces = double2_to_vector2d_translator(forces);
    universe.velocities = double2_to_vector2d_translator(velocities);
    universe.positions = double2_to_vector2d_translator(positions);


}



/*
> implementation decision: just convert one iteration of naive sequential
to cuda kernel + add schelifenbedingung as an if.


ERROR: GRAVITIONAL FORCE IS NOT WELL DEFINED IN GPU USE! -> changed


 */
__global__
void calculate_forces_kernel(std::uint32_t num_bodies, double2* d_positions, double* d_weights, double2* d_forces){

//at the same time himmelskörper index, as
    //body_id_x in sequential function
int thread_himmelskorper_index = blockIdx.x * blockDim.x + threadIdx.x;

//last legal move= num_bodies -1, right?
if (thread_himmelskorper_index < num_bodies){

    double2 body_position = d_positions[thread_himmelskorper_index];
    double body_mass = d_weights[thread_himmelskorper_index];

     double2 applied_force_vector = make_double2(0.0, 0.0);

     for (int distant_body_idx = 0; distant_body_idx < num_bodies; ++distant_body_idx) {

        //skip current korper:


        if (thread_himmelskorper_index == distant_body_idx){
            continue;
        }


          double2 distant_body_position = d_positions[distant_body_idx];
        double distant_body_mass = d_weights[distant_body_idx];

        double2 direction_vector = make_double2(
                distant_body_position.x - body_position.x,
                distant_body_position.y - body_position.y
            );



        double distance = sqrt(pow(direction_vector.x, 2) + pow(direction_vector.y, 2));

        double force_magnitude = compute_gravitational_force_device(body_mass, distant_body_mass, distance);


        double2 force_vector;
        force_vector.x = (direction_vector.x / distance) * force_magnitude;
        force_vector.y = (direction_vector.y / distance) * force_magnitude;


        applied_force_vector.x += force_vector.x;
        applied_force_vector.y += force_vector.y;

     }//end of for for distant bodies

    //save results
     d_forces[thread_himmelskorper_index] = applied_force_vector;

}//end of schleifenbedingung checker if

}


/*
implementation decision: how to choose block/grid/tile etc. sizes?
    in vorlesung = tile_width was an pre-defined macro and we haven't talked about it much

    reseaarch: a sweet spoot between 128-512 (must be dividable by warp size 32)
    -> occupancy api: too complez, we have no performance criteria, turn back if comes in upcoming aufgaben

    -> naive and sequential shouldnt return different results -> should we do one by one???? TEST!!!!


    -> DO WE NEED POINTER CASTING???

 */
void NaiveCudaSimulation::calculate_forces(Universe& universe, void* d_positions, void* d_weights, void* d_forces){


    //each block would calculate 256 results, 256 bodies
     int thread_amount_pro_block = 256;

     //how many blocks are needed given result array size:
    int amount_blocks = (universe.num_bodies + thread_amount_pro_block - 1) / thread_amount_pro_block;



     calculate_forces_kernel<<<amount_blocks, thread_amount_pro_block>>>(
        universe.num_bodies, 
        (double2*)d_positions, 
        (double*)d_weights,     
        (double2*)d_forces
    );

        cudaDeviceSynchronize();


}


__device__
double2 calculate_acceleration_in_double2(double2 applied_force, double mass) {
    
    // calculate acceleration 
    // a = F / m
    return make_double2(applied_force.x / mass, applied_force.y / mass);
}


__device__
double2 calculate_velocity_in_double2(double2 base_velocity, double2 acceleration, double time_in_seconds) {
    
    // v = v0 + a * t
    return make_double2(
        base_velocity.x + acceleration.x * time_in_seconds,
        base_velocity.y + acceleration.y * time_in_seconds
    );
}


/*

BIG TODO: DO WE INJECT "static const double epoch_in_seconds = 2.628e+6;" HERE?
CHECK IN A ENVIROMENT WITH MORE STABLE INTELLISENSE!!!!!!!!!!
-> currently white!!!!!!!!!!! -> added constant.h as header, do i miss sth?


*/
__global__
void calculate_velocities_kernel(std::uint32_t num_bodies, double2* d_forces, double* d_weights, double2* d_velocities){

    int thread_himmelskorper_index = blockIdx.x * blockDim.x + threadIdx.x;
    

    if(thread_himmelskorper_index< num_bodies){

    double2 force = d_forces[thread_himmelskorper_index];
    double mass = d_weights[thread_himmelskorper_index];

    
    double2 acceleration = calculate_acceleration_in_double2(force, mass);
    d_velocities[thread_himmelskorper_index] = calculate_velocity_in_double2(d_velocities[thread_himmelskorper_index], acceleration, epoch_time_in_seconds_device);



    }//end of schleifenbedingung if 


}

void NaiveCudaSimulation::calculate_velocities(Universe& universe, void* d_forces, void* d_weights, void* d_velocities){

  //each block would calculate 256 results, 256 bodies
     int thread_amount_pro_block = 256;

     //how many blocks are needed given result array size:
    int amount_blocks = (universe.num_bodies + thread_amount_pro_block - 1) / thread_amount_pro_block;


 calculate_velocities_kernel<<<amount_blocks , thread_amount_pro_block>>>(
        universe.num_bodies, 
        (double2*)d_forces,  // Use the passed argument
        (double*)d_weights,  // Use the passed argument
        (double2*)d_velocities  // Use the passed argument
    );

    cudaDeviceSynchronize();




}

__global__
void calculate_positions_kernel(std::uint32_t num_bodies, double2* d_velocities, double2* d_positions){

int thread_himmelskorper_index = blockIdx.x * blockDim.x + threadIdx.x;


  if(thread_himmelskorper_index< num_bodies){


    double2 velocity = d_velocities[thread_himmelskorper_index];
    double2 position = d_positions[thread_himmelskorper_index];

// calculate movement
        // s = v * t
     double2 movement;
    movement.x = velocity.x * epoch_time_in_seconds_device;
    movement.y = velocity.y * epoch_time_in_seconds_device;

// calculate new position
        // p` = p0 + s 
     double2 new_position;
    new_position.x = position.x + movement.x;
    new_position.y = position.y + movement.y;


     d_positions[thread_himmelskorper_index] = new_position;


  }//end of schleifenbedingung if


}

void NaiveCudaSimulation::calculate_positions(Universe& universe, void* d_velocities, void* d_positions){


//each block would calculate 256 results, 256 bodies
     int thread_amount_pro_block = 256;

     //how many blocks are needed given result array size:
    int amount_blocks = (universe.num_bodies + thread_amount_pro_block - 1) / thread_amount_pro_block;


 calculate_positions_kernel<<<amount_blocks, thread_amount_pro_block>>>(
        universe.num_bodies, 
        (double2*)d_velocities, 
        (double2*)d_positions
    );
    cudaDeviceSynchronize();


}



/* implementation decisions:

ANFORDERUNG: gpu speicher allocation + freeing

sequantial workflow for ONE EPOCH

-> calculate forces
->c velocities
-> c positions
-> +1 epoch

if create intermediate plots:,
     if((universe.current_simulation_epoch % plot_intermediate_epochs) == 0){
            plotter.add_bodies_to_image(universe);
            plotter.write_and_clear();
        }


epochS : repeat given times


-> where allocating memory

-> where freeing



-> memory allocated, pointer of pointer type of pointers MUST BE CASTED TO VOID!
    -> should we make the initial pointers null pointer?
        -> i did, it wasnt like that in vorlesung example but research says better practice
            -> TODO IF ERROR: make them unitialized




 */

// this function would be called by CPU, 
void NaiveCudaSimulation::simulate_epochs(Plotter& plotter, Universe& universe, std::uint32_t num_epochs, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs){

//ASSIGN MEMORY

double* d_weights = nullptr;
double2* d_forces = nullptr;
double2* d_velocities = nullptr;
double2* d_positions = nullptr;

NaiveCudaSimulation::allocate_device_memory(universe, 
    (void**)&d_weights, 
    (void**)&d_forces, 
    (void**)&d_velocities, 
    (void**)&d_positions
    );

//COPY DATA FROM HOST TO DEVICE

NaiveCudaSimulation::copy_data_to_device(universe, d_weights, d_forces, d_velocities, d_positions);

// DO THE SIMULATIONS

 for(int i = 0; i < num_epochs; i++){
        NaiveCudaSimulation::simulate_epoch(plotter, universe, create_intermediate_plots, plot_intermediate_epochs, d_weights, d_forces, d_velocities, d_positions);
    }


//copy data back

NaiveCudaSimulation::copy_data_from_device(universe, d_weights, d_forces, d_velocities, d_positions);

//free memory
NaiveCudaSimulation::free_device_memory((void**)&d_weights, (void**)&d_forces, (void**)&d_velocities, (void**)&d_positions);

}





__global__
void get_pixels_kernel(std::uint32_t num_bodies, double2* d_positions, std::uint8_t* d_pixels, std::uint32_t plot_width, std::uint32_t plot_height, double plot_bounding_box_x_min, double plot_bounding_box_x_max, double plot_bounding_box_y_min, double plot_bounding_box_y_max){
    
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= num_bodies) return;

    double2 pos = d_positions[i];

    if (pos.x >= plot_bounding_box_x_min && pos.x <= plot_bounding_box_x_max &&
        pos.y >= plot_bounding_box_y_min && pos.y <= plot_bounding_box_y_max) {
        
        int pixel_x = static_cast<int>(((pos.x - plot_bounding_box_x_min) / 
                      (plot_bounding_box_x_max - plot_bounding_box_x_min)) * (plot_width - 1));
        int pixel_y = static_cast<int>(((pos.y - plot_bounding_box_y_min) / 
                      (plot_bounding_box_y_max - plot_bounding_box_y_min)) * (plot_height - 1));

        if (pixel_x >= 0 && pixel_x < plot_width && pixel_y >= 0 && pixel_y < plot_height) {
            int pixel_index = pixel_y * plot_width + pixel_x;
            d_pixels[pixel_index] = 255;
        }
    }
}

std::vector<std::uint8_t> NaiveCudaSimulation::get_pixels(std::uint32_t plot_width, std::uint32_t plot_height, BoundingBox plot_bounding_box, void* d_positions, std::uint32_t num_bodies){
    std::vector<std::uint8_t> pixels(plot_width * plot_height, 0);

    uint8_t* d_pixels_void;
    size_t pixel_data_size = plot_width * plot_height * sizeof(std::uint8_t);

    parprog_cudaMalloc(reinterpret_cast<void**>(&d_pixels_void), pixel_data_size);
    cudaMemset(d_pixels_void, 0, pixel_data_size);


    const int threadsPerBlock = 256;
    const int numBlocks = (plot_width * plot_height + threadsPerBlock - 1) / threadsPerBlock;

    get_pixels_kernel<<<numBlocks, threadsPerBlock>>>(
        num_bodies,                              
        static_cast<double2*>(d_positions),      
        d_pixels_void,                                
        plot_width, plot_height,                
        plot_bounding_box.x_min, plot_bounding_box.x_max,
        plot_bounding_box.y_min, plot_bounding_box.y_max
    );
    

    cudaDeviceSynchronize();

    parprog_cudaMemcpy(pixels.data(), d_pixels_void, pixel_data_size, cudaMemcpyDeviceToHost);

    parprog_cudaFree(d_pixels_void);

    return pixels;

}

__global__
void compress_pixels_kernel(std::uint32_t num_raw_pixels, std::uint8_t* d_raw_pixels, std::uint8_t* d_compressed_pixels){

}

void NaiveCudaSimulation::compress_pixels(std::vector<std::uint8_t>& raw_pixels, std::vector<std::uint8_t>& compressed_pixels){

}

void NaiveCudaSimulation::simulate_epoch(Plotter& plotter, Universe& universe, bool create_intermediate_plots, std::uint32_t plot_intermediate_epochs, void* d_weights, void* d_forces, void* d_velocities, void* d_positions){
    calculate_forces(universe, d_positions, d_weights, d_forces);
    calculate_velocities(universe, d_forces, d_weights, d_velocities);
    calculate_positions(universe, d_velocities, d_positions);

    universe.current_simulation_epoch++;
    if(create_intermediate_plots){
        if(universe.current_simulation_epoch % plot_intermediate_epochs == 0){
            std::vector<std::uint8_t> pixels = get_pixels(plotter.get_plot_width(), plotter.get_plot_height(), plotter.get_plot_bounding_box(), d_positions, universe.num_bodies);
            plotter.add_active_pixels_to_image(pixels);

            // This is a dummy to use compression in plotting, although not beneficial performance-wise
            // ----
            // std::vector<std::uint8_t> compressed_pixels;
            // compressed_pixels.resize(pixels.size()/8);
            // compress_pixels(pixels, compressed_pixels);
            // plotter.add_compressed_pixels_to_image(compressed_pixels);
            // ----

            plotter.write_and_clear();
        }
    }
}

void NaiveCudaSimulation::calculate_forces_kernel_test_adapter(std::uint32_t grid_dim, std::uint32_t block_dim, std::uint32_t num_bodies, void* d_positions, void* d_weights, void* d_forces){
    // adapter function used by automatic tests. DO NOT MODIFY.
    dim3 blockDim(block_dim);
    dim3 gridDim(grid_dim);
    calculate_forces_kernel<<<gridDim, blockDim>>>(num_bodies, (double2*) d_positions, (double*) d_weights, (double2*) d_forces);
}

void NaiveCudaSimulation::calculate_velocities_kernel_test_adapter(std::uint32_t grid_dim, std::uint32_t block_dim, std::uint32_t num_bodies, void* d_forces, void* d_weights, void* d_velocities){
    // adapter function used by automatic tests. DO NOT MODIFY.
    dim3 blockDim(block_dim);
    dim3 gridDim(grid_dim);
    calculate_velocities_kernel<<<gridDim, blockDim>>>(num_bodies, (double2*) d_forces, (double*) d_weights, (double2*) d_velocities);
}

void NaiveCudaSimulation::calculate_positions_kernel_test_adapter(std::uint32_t grid_dim, std::uint32_t block_dim, std::uint32_t num_bodies, void* d_velocities, void* d_positions){
    // adapter function used by automatic tests. DO NOT MODIFY.
    dim3 blockDim(block_dim);
    dim3 gridDim(grid_dim);
    calculate_positions_kernel<<<gridDim, blockDim>>>(num_bodies, (double2*) d_velocities, (double2*) d_positions);
}
