#pragma once

#include "structures/vector2d.h"
#include "structures/universe.h"
#include "quadtreeNode.h"


class Quadtree{
public: 
    Quadtree(Universe& universe, BoundingBox bounding_box, std::int8_t construct_mode);
    ~Quadtree();

    std::vector<QuadtreeNode *> construct(Universe &universe, BoundingBox BB, std::vector<std::int32_t> body_indices);
    void initialize_a_leaf_node(QuadtreeNode *current_quadrant_as_node, int32_t index, Universe &universe);
    std::vector<QuadtreeNode*> construct_task(Universe& universe, BoundingBox BB, std::vector<std::int32_t> body_indices);
    std::vector<QuadtreeNode*> construct_task_with_cutoff(Universe& universe, BoundingBox& BB, std::vector<std::int32_t>& body_indices);

    void calculate_cumulative_masses();
    void calculate_center_of_mass();
    QuadtreeNode* root = nullptr;

    std::vector<BoundingBox> get_bounding_boxes(QuadtreeNode* qtn);

    std::vector<std::int32_t> body_indices_in_a_given_bounding_box( Universe& universe, BoundingBox BB,std::vector<std::int32_t> pre_indices);
    std::vector<std::int32_t> preindice_vector_creator_for_universe(std::size_t size);
};


