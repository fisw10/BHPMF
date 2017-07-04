//
//  Latentnode.h
//
//
//  Created by Farideh Fazayeli on 7/15/13.
//
//

#ifndef ____latentNode__
#define ____latentNode__

#include <iostream>
#include <memory>
#include "utillity.h"
using namespace std;
class NodeLatent{
    
private:
    // A vector of size num_latent_features_ storing latent factor for this node
    std::unique_ptr<double[]> latent_factor_;
    // Store the old value for gradient descent update
    std::unique_ptr<double[]> old_latent_factor_;
    // A vector of size num_latent_features_ storing latent factor for parents
    // and children (add together)
    std::unique_ptr<double[]> prior_;
    // A vecor of size num_latent_features_ storing the gradient for the
    // gradient descent update
    std::unique_ptr<double[]> gradient_;
    
    // Number of parents
    int num_parent_;
    // Number of children
    int num_child_;
    // Index of last children inserted to the tree
    int last_child_;
    // Size of latent factor $k$
    int num_latent_features_;
    // Number of observed value for this node
    int num_observ_;
    // A vector of size num_observe_ containing the observed value for this node
    vector<float> observ_;
    // A vector of index of nodes in other tree which this node have an observed
    // value. It is used for retrieve the corresponding nodes from other tree.
    // For example, if the current node is a plant in Utree, then observ_idx_
    // is index of all traits in which this plant have an observed trait value.
    vector<int> observ_idx_;
    // A vector pointing to parents of the current node
    NodeLatent* parents_;
    // A vector pointing to the children of the current node
    NodeLatent** children_;
    
    NodeLatent() {
        num_parent_ = 0;
        num_child_ = 0;
    }
    
public:
    ~NodeLatent();
    NodeLatent(int child, int feat, int lev);
    
    // Getting num parent and num children and set numPrior
    NodeLatent(int par, int child, int parIdx, NodeLatent ***tree, int feat,
               int lev, int parLevel);
    inline double *GetLatentFactor() {
        return latent_factor_.get();
    }
    
    // Setting the observed value and corresponding node in other tree.
    // For example if current node is a node (e.g., plant) in Utree, it sets the
    // observed value (obs) and index of node (e.g., trait) in Vtree.
    // @obs: the observed value.
    // @ind: index of the pair node in other tree.
    void SetObserv(float obs, int ind);
    
    /**************************************************
     List of functions used in gibbs sampling Inference
     **************************************************/
    // Update the prior information used for gibbs sampling
    void UpdatePrior();
    // Sample latent factors for the current node, given other nodes
    // using the gibbs condition
    // @v_tree: a pointer to other tree. For example if current node is in u_tree
    // v_tree is the pointer to the v_tree
    // @sig_inv: inverse of sigma^2 for observations
    // @sig_u_inv: inverse of sigma_u^2
    // @level: the level at which current node is located
    void GibbsUpdate(NodeLatent*** v_tree, double sig_inv, double sig_u_inv,
                     int level);

    /********************************************
     Functions used in Gradient Descent Inference
     TODO: Add Gradient Descent in HPMF class
     ********************************************/
    // Update the gradient used for gradient descent update
    void UpdateGradient(int freq, double err, double *vv);
    // Update the latent factor used for gradient descent update
    void UpdateLatentFactor(double coef, double momentum, double epsilon);
    // Update the prior information used for gradient descent update
    void UpdatePrior(int mode);
};

#endif /* defined(____latentNode__) */

