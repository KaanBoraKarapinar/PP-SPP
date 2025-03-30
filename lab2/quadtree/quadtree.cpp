#include "quadtree.h"

#include "quadtreeNode.h"
#include <set>
#include <algorithm>
#include <stdexcept>
#include <omp.h>

Quadtree::Quadtree(Universe& universe, BoundingBox bounding_box, std::int8_t construct_mode){


    //quadtree objesi oluşturduğunda otomatik olarak leaflarında himmelskörper olan bir tree yapısı döndürür
    //adım başı dörte ayırır elimizdeki bounding boxu (bir bounding box 4 quadrat)
    //tahmin: leafların konumu aynı değil, her bir quadratte höchstens 1 olunca (0 KABUL EDİLİYOR) duruyor sistem

    //root node initialization
        //children = calculated ones 
    root = new QuadtreeNode(bounding_box);

    //contains direct pointers to 4 initial quadranten, recursively filled with children of quadranten
    std::vector<QuadtreeNode*> rootNodeChildrenVector;


    //implementation decision: all body indices?  we must calculate them paralelly, and then in recursive steps only check the giben list
            //initially: all elements, todo: can we assume it is 0to last element?
    std::vector<std::int32_t> pre_indices = preindice_vector_creator_for_universe(universe.positions.size());

    if(construct_mode == 0){
       
        rootNodeChildrenVector = construct(universe,bounding_box,body_indices_in_a_given_bounding_box(universe,bounding_box,pre_indices));
    }

    if(construct_mode == 1){
       
        rootNodeChildrenVector = construct_task(universe,bounding_box,body_indices_in_a_given_bounding_box(universe,bounding_box,pre_indices));
    }

    if(construct_mode == 2){
       
       //-> with cutoff requries body_indices as reference, so i create and keep it here

       std::vector<std::int32_t> body_indices_boundingbox = body_indices_in_a_given_bounding_box(universe,bounding_box,pre_indices);
        rootNodeChildrenVector = construct_task_with_cutoff(universe,bounding_box,body_indices_boundingbox);
    }

    root->children = rootNodeChildrenVector;



}


/*
basically adds to vector 0 to index of last element
*/
std::vector<std::int32_t> Quadtree::preindice_vector_creator_for_universe(std::size_t size) {
    std::vector<std::int32_t> indices;

    for (std::int32_t i = 0; i < static_cast<std::int32_t>(size); ++i) {
        indices.push_back(i);
    }

    return indices;
}

/*
This parallel version (PROBABLY) does not return indices in ascending order, it is randomly mixed!
*/
std::vector<std::int32_t> Quadtree::body_indices_in_a_given_bounding_box( Universe& universe, BoundingBox BB,std::vector<std::int32_t> pre_indices) {

    std::vector<std::int32_t> current_indices;

    // OpenMP kullanarak paralel bir for döngüsü
    #pragma omp parallel
    {
        // Paralel ortamda thread-local bir vektör oluştur
        std::vector<std::int32_t> local_indices;

        #pragma omp for nowait
        for (std::size_t i = 0; i < pre_indices.size(); ++i) {
            const auto& position = pre_indices[i];
            if (BB.contains(universe.positions[position])) {
                // Thread-local vektöre ekle
                local_indices.push_back(position);
            }
        }

        // Thread-local sonuçları ana vektöre birleştir
        #pragma omp critical
        {
            current_indices.insert(current_indices.end(), local_indices.begin(), local_indices.end());
        }
    }

    return current_indices;
}


/*

SERIAL VERSION


std::vector<std::int32_t> body_indices_in_a_given_bounding_box( Universe& universe,  BoundingBox BB, std::vector<std::int32_t> pre_indices){

     std::vector<std::int32_t> current_indices;

    
      for (const auto& position : pre_indices) {
        if (BB.contains(universe.positions[position])) {
            current_indices.push_back(position);
        }
    }

    return current_indices;

}*/

Quadtree::~Quadtree(){


/*
goals 
-> delete root
-> trigger destruction of child nodes till end -> delegated that to the quadtreenode
*/
  
   if (root) {
        delete root;
        root = nullptr;
    }

}

void Quadtree::calculate_cumulative_masses(){
    root->calculate_node_cumulative_mass();
}

void Quadtree::calculate_center_of_mass(){
    root->calculate_node_center_of_mass();
}


std::vector<QuadtreeNode*> Quadtree::construct(Universe& universe, BoundingBox BB, std::vector<std::int32_t> body_indices){
    /*
    infos
    -> serial (RECURSIVE!)
    -> BoundingBox quadrant = bounding_box.get_quadrant(quadrant_id); -> CREATES A BOUNDING BOX WITHIN BORDRERS OF AN QUADRANT

    current thoughts:
    -> get amount of

        

    */
    
    std::vector<QuadtreeNode*> vector_of_created_nodes;

        //for all 4 quadrants
    for(int i = 0; i <= 3; ++i){

            //convert current quadrant to a bounding box, create a node
        
        BoundingBox current_quadrant_as_boundingbox = BB.get_quadrant(i);
        QuadtreeNode* current_quadrant_as_node = new QuadtreeNode(current_quadrant_as_boundingbox);

            //get number of contained bodies in quadrant(in new bb)

        std::vector<std::int32_t> body_indices_in_current_quadrant = body_indices_in_a_given_bounding_box(universe,current_quadrant_as_boundingbox,body_indices);

        std::int32_t amountOfBodiesInCurrentQuadrant = body_indices_in_current_quadrant.size();


            //TODO: BIG QUESTION:
                // ARE ALL LEAVES AT SAME LEVEL, DO WE KEEP DIVIDING IT?
                //IMPLEMENTATION HERE -> LEAVES AT DIFFERENT LEVELS

            //if size more than 1, -> INNERE KNOTE
                // we must contiune breaking down
                // assign children 
                //body_identifier = -1, initialized in quadtreenode.h no need to reassign

            if(amountOfBodiesInCurrentQuadrant > 1){
                
                //assign children
                current_quadrant_as_node->children =  construct(universe,current_quadrant_as_boundingbox,body_indices_in_current_quadrant);
                  vector_of_created_nodes.push_back(current_quadrant_as_node);
            }

                //if size = 1, leaf
                    //no children, assign body identifier

            if(amountOfBodiesInCurrentQuadrant == 1){

                initialize_a_leaf_node(current_quadrant_as_node, body_indices_in_current_quadrant[0], universe);

                vector_of_created_nodes.push_back(current_quadrant_as_node);

            }


            //if size = 0, also leaf, currently no specific order
            //if(amountOfBodiesInCurrentQuadrant = 0){} 

          

   }
    
    return vector_of_created_nodes;
}

void Quadtree::initialize_a_leaf_node(QuadtreeNode *current_quadrant_as_node, int32_t index, Universe &universe)
{
    current_quadrant_as_node->body_identifier = index;
    current_quadrant_as_node->cumulative_mass = universe.weights[index];
    current_quadrant_as_node->center_of_mass = universe.positions[index];

    current_quadrant_as_node->cumulative_mass_ready = true;
    current_quadrant_as_node->center_of_mass_ready = true;
}

std::vector<QuadtreeNode*> Quadtree::construct_task(Universe& universe, BoundingBox BB, std::vector<std::int32_t> body_indices) {
    std::vector<QuadtreeNode*> vector_of_created_nodes;

 
    #pragma omp parallel
    {
        // Shared vector and parallel context
        #pragma omp single // Only one thread initializes tasks
        {
            // For all 4 quadrants
            for (int i = 0; i <= 3; ++i) {
                // Create a new task for each quadrant
                #pragma omp task firstprivate(i) shared(vector_of_created_nodes)
                {
                    
                    BoundingBox current_quadrant_as_boundingbox = BB.get_quadrant(i);
                    QuadtreeNode* current_quadrant_as_node = new QuadtreeNode(current_quadrant_as_boundingbox);

                    
                    std::vector<std::int32_t> body_indices_in_current_quadrant =
                        body_indices_in_a_given_bounding_box(universe, current_quadrant_as_boundingbox, body_indices);

                    std::int32_t amountOfBodiesInCurrentQuadrant = body_indices_in_current_quadrant.size();

                    
                    if (amountOfBodiesInCurrentQuadrant > 1) {
                        current_quadrant_as_node->children =
                            construct(universe, current_quadrant_as_boundingbox, body_indices_in_current_quadrant);
                    }

                    
                    if (amountOfBodiesInCurrentQuadrant == 1) {
                       initialize_a_leaf_node(current_quadrant_as_node, body_indices_in_current_quadrant[0], universe);
                    }

                    // Add node to the shared vector (critical section)
                    #pragma omp critical
                    {
                        if(amountOfBodiesInCurrentQuadrant != 0){
                            vector_of_created_nodes.push_back(current_quadrant_as_node);
                        }
                        
                    }
                }
            }
        }
    }

    return vector_of_created_nodes;
}

std::vector<QuadtreeNode*> Quadtree::construct_task_with_cutoff(Universe& universe, BoundingBox& BB, std::vector<std::int32_t>& body_indices){


    /*
    
    how to pick a cutoff criteria?
    -> body indices parameter is passed here by refrence, not by 
    
    */

      std::vector<QuadtreeNode*> vector_of_created_nodes;


    int cutoff= 10;

    
    #pragma omp parallel
    {
        #pragma omp single
        {
            
            for (int i = 0; i <= 3; ++i) {
                #pragma omp task firstprivate(i) shared(vector_of_created_nodes)
                {
                    
                    BoundingBox current_quadrant_as_boundingbox = BB.get_quadrant(i);
                    QuadtreeNode* current_quadrant_as_node = new QuadtreeNode(current_quadrant_as_boundingbox);

                    
                    std::vector<std::int32_t> body_indices_in_current_quadrant =
                        body_indices_in_a_given_bounding_box(universe, current_quadrant_as_boundingbox, body_indices);

                    std::int32_t amountOfBodiesInCurrentQuadrant = body_indices_in_current_quadrant.size();

                    // Cut-off criteria implementation
                    if (amountOfBodiesInCurrentQuadrant > cutoff) {
                        // Alt düğümleri paralel bir şekilde oluştur
                        current_quadrant_as_node->children =
                            construct_task_with_cutoff(universe, current_quadrant_as_boundingbox, body_indices_in_current_quadrant);
                    } else if (amountOfBodiesInCurrentQuadrant > 1) {
                        // Eğer cutoff kriteri sağlanmazsa, seri bir şekilde devam et
                        current_quadrant_as_node->children =
                            construct(universe, current_quadrant_as_boundingbox, body_indices_in_current_quadrant);
                    }

                    // Yaprak düğümse, body_identifier ata
                    if (amountOfBodiesInCurrentQuadrant == 1) {
                       initialize_a_leaf_node(current_quadrant_as_node, body_indices_in_current_quadrant[0], universe);
                    }

                    // Node'u shared vektöre ekle
                    #pragma omp critical
                     if(amountOfBodiesInCurrentQuadrant != 0){
                            vector_of_created_nodes.push_back(current_quadrant_as_node);
                        }
                }
            }
        }
    }

    return vector_of_created_nodes;
}



std::vector<BoundingBox> Quadtree::get_bounding_boxes(QuadtreeNode* qtn){
    // traverse quadtree and collect bounding boxes
    std::vector<BoundingBox> result;
    // collect bounding boxes from children
    for(auto child: qtn->children){
        for(auto bb: get_bounding_boxes(child)){
            result.push_back(bb);
        }
    }
    result.push_back(qtn->bounding_box);
    return result;
}








