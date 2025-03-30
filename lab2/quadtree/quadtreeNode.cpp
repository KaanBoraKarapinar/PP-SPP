#include "quadtreeNode.h"

#include <iostream>


/*
when once calculated, would recalling does not update the value and return the last one.
*/
double QuadtreeNode::calculate_node_cumulative_mass(){
    
    //todo: is it wrong? if calculated once, we doesnt update
    if(cumulative_mass_ready){

        return cumulative_mass;
    }


    int local_cumulative_mass = 0.0;

    if (children.empty()) {
        
        //if leaf do nothing, this functionality has been delegated to quadtree
        return cumulative_mass;
    } else {
        // İç düğümse, tüm çocuk düğümlerin kütlelerini topla
        for (auto child : children) {
            child->calculate_node_cumulative_mass(); // Önce çocuğun kütlesini hesapla
            local_cumulative_mass += child->cumulative_mass;
        }
    }

    cumulative_mass = local_cumulative_mass;

    cumulative_mass_ready = true; // Kütle hesaplandı ve geçerli


    
    return cumulative_mass;

}

QuadtreeNode::QuadtreeNode(BoundingBox arg_bounding_box)   : bounding_box(arg_bounding_box){
//hangi quadrant olduğu verisini nerede saklıyoruz? -> düşünülmemiş, gerekirse eklersin
}

QuadtreeNode::~QuadtreeNode(){
 
 
 for (QuadtreeNode* child : children) {
        if (child) {
            delete child;
        }
    }


//todo: do we need?
//children.clear();


}

Vector2d<double> QuadtreeNode::calculate_node_center_of_mass(){
    
    //todo: is it wrong? if calculated once, we doesnt update
    if(center_of_mass_ready){

        return center_of_mass;

    }

    Vector2d<double> local_vector(0.0,0.0);

    if (children.empty()) {
        
        //if leaf do nothing, this functionality has been delegated to quadtree
        return center_of_mass;
    } else {
        // İç düğümse, tüm çocuk düğümlerin kütlelerini topla
        for (auto child : children) {
            child->calculate_node_center_of_mass();
             Vector2d<double> temp = child->center_of_mass * child->cumulative_mass;
             local_vector = local_vector + temp;

        }

        local_vector = local_vector / calculate_node_cumulative_mass();
        center_of_mass = local_vector;
        center_of_mass_ready = true;
    }
    
    return center_of_mass;

}

