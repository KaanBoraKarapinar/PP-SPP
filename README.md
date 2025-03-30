# PP-SPP
Paralell Programming Coursework Solutions

This project is the result of a three-part practicum in the Parallel Programming module at Technische Universität Darmstadt. The overall goal was to model gravitational interactions between massive bodies in space using numerical simulation, and to iteratively optimize this model across different levels of parallelism.

Over the course of the semester, the implementation evolved from a basic sequential version into a scalable, high-performance system leveraging both CPU and GPU parallelism. Each practicum introduced more advanced techniques while retaining strict constraints on software quality, performance, and correctness.

All components were implemented using modern, standard-compliant C++ and CUDA, and tested on TU Darmstadt’s Lichtenberg high-performance computing system.

## Disclaimer on Repository Scope and Copyright

__This repository includes only those files that were created or modified by our team of three as part of the Parallel Programming practicum at TU Darmstadt (Winter Semester 2024/25).__

__The full software framework, test infrastructure, and utility code provided by the university are not included in this repository to respect the intellectual property and redistribution policies of TU Darmstadt__


## Praktikum 1: Sequential N-Body Simulation
The first phase established the foundation of the simulation. The focus was on implementing core physical behavior and designing a clean, modular architecture.

* Developed a full gravitational simulation from scratch using Newtonian mechanics
* Implemented time-stepped motion updates: force → acceleration → velocity → position
* Built core types such as Vector2d, BoundingBox, and Universe using templates, operator overloading, and strong type safety
* Emphasized const-correctness, exception safety (noexcept), and diagnostic attributes ([[nodiscard]])
* Output visualization through bitmap plotting, enabling visual inspection of body trajectories
* Code was structured for readability, reusability, and unit testability, complying with a fixed build and runtime environment

## Praktikum 2: CPU-Based Parallelization with OpenMP
* Applied OpenMP to parallelize previously sequential simulation steps (forces, velocities, positions, bounding boxes)
* Designed and implemented a Quadtree spatial decomposition structure to group distant bodies efficiently
* Implemented multiple versions of Quadtree construction: serial, task-based, and cutoff-optimized parallel
* Modeled mass aggregation and center-of-mass calculations recursively within the tree
* All functions were bound to strict performance constraints (e.g., < 1-minute per functional test) and required output consistency with predefined formats

## Praktikum 3: CUDA-Based GPU Acceleration
* Implemented memory management and host-device data transfer using wrapped CUDA allocation and copy functions
* Wrote GPU kernels for force calculation, velocity integration, and position updates, assigning one thread per body
* Converted host-side types to CUDA-compatible formats (Vector2d ↔ double2) during memory transfer
* Ensured GPU simulation results matched the CPU-based sequential outputs in both accuracy and structure
