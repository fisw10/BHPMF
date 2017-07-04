//
//  HPMF.h
//
//
//  Created by Farideh Fazayeli on 7/8/13.
//
//

#ifndef ____HPMF__
#define ____HPMF__

#include <iostream>
#include <string>
#include <algorithm>    // std::move_backward
#include <fenv.h>
//#include <array>        // std::array
//#include <random>       // std::default_random_engine
//#include <chrono>       // std::chrono::system_clock

#include <R.h>      // R functions
#include <Rmath.h>  // R math
#include <numeric>
//#include <Rinternals.h>

#include "utillity.h"
#include "latentNode.h"

#ifdef _OPENMP
#include "omp.h"
#define HPMF_num_threads 10
#else
#define omp_get_thread_num() 0
#define omp_get_wtime() 0
#endif  /* #ifndef _OPENMP */

#include <vector>



class HPMF{
private:
    /*****************************
     Missing matrix parameters
     ****************************/
    int numPlants; //number of plants
    // number of columns in X matrix, for example number of traits
    const int num_columns_;
    // size of latent factor $k$
    int num_latent_features_;
    
    // latent matrices U represented as a tree
    // number of nodes at each level is equal to number of rows at corresponding
    // matrix at that level.
    // i.e., for a given "level" and one of its "node"
    // u_tree_[level][node] is a vector of size num_latent_features_
    NodeLatent ***u_tree_;
    // latent matrices V with the same structure as u_tree_
    NodeLatent ***v_tree_;
    // A vector containing number of nodes at each level
    vector<double> num_nodes_per_level_;
    
    vector<int> nodePerLevel; //????????
    
    /*****************************
     Hierarchy information:
     num_hierarchy_level_: total number of levels in hierarchy
     for example for Plant Taxonomy: Phylogenetic Group, Famlily, Genus, Species
     num_hierarchy_level_ = 4
     used_hierarchy_level_ = 4 : use all of taxonomy information
                                (Phylogenetic Group + Famlily + Genus + Species)
     used_hierarchy_level_ = 3 : use 3 levels of taxonomy information
                                (Phylogenetic Group + Famlily + Genus )
     used_hierarchy_level_ = 2 : use 2 levels of taxonomy information
                                (Phylogenetic Group + Famlily)
     used_hierarchy_level_ = 1 : use 1 level of taxonomy information
                                (Phylogenetic Group)
     used_hierarchy_level_ = 0 : No taxonomy information
     predict_level_ = 5 : prediction at the Plant level
     predict_level_ = 4 : prediction at the Species level
     predict_level_ = 3 : prediction at the Genus level
     predict_level_ = 2 : prediction at the Family level
     predict_level_ = 1 : prediction at the Phylogenetic Group level
     Note: predict_level_ > used_hierarchy_level_
     *****************************/
    int predict_level_;
    int num_hierarchy_level_;
    int used_hierarchy_level_;
    
    // hiearchy_ is a vector of size sum(num_nodes_per_level_) represented as
    // the long an array link list for all nodes in the tree.
    // Consider hierarchy_ as the following vector
    // ------------------------------------------------------------------------
    // | Level 1 |  Level 2   | Level 3             | Level 4                  |
    // ------------------------------------------------------------------------
    //    ^           |
    //    !-----------!
    // Level 1: hierarchy[0: num_nodes_per_level_[1]] = 0 as nodes at level 1
    // have no parents.
    // Level 2: hierarchy[num_nodes_per_level_[1]:
    //                    num_nodes_per_level_[1]+num_nodes_per_level_[2]]
    // contains the parents for nodes at level 2, for example if parent of a
    // node is 5 it is the 5th element at Level 1.
    // Level 3: similarly contains the parents for nodes at level 3, for example
    // if parent of a node is 3 it is the 3rd element at Level 2.
    // Level 4: similarly contains the parents for nodes at level 4.
    matVecInt hierarchy_;
    // A vector of similar structure and length as hierarchy_.
    // Contains number of parents for each node.
    // TODO: Currently BHPMF only supports single parent, extend it to multiple
    // Parents, require a better structure for hierarchy_.
    vector<int> num_parents_;
    // A vector of similar structure and length as hierarchy_.
    // Contains number of children for each node.
    vector<int> num_children_;
    
    /**************************
     Gibbs sampling Parameters
     ****************************/
    // gaps parameter in gibbs sampling
    int gaps_;
    // burn parameter in gibbs sampling
    int burn_;
    // number of total samples in gibbs sampling
    int num_samples_;
    // effective number of samples in gibbs smapling
    int num_effective_samples_;
    
    /**************************
     IO parameters
     *************************/
    // the fold_id if using cross validation to provide average RMSE
    // fold_id=1 for gap filling
    int fold_id_;
    // input directory containing preprocessed data
    const char* input_dir_;
    // If set TRUE, save the prediction for test data
    int saveFileFlag;
    // if set TRUE, gap fill the whole matrix, write both mean and std
    // for whole dataset in two files
    int outWholeFlag;
    // The file path (including directory and file name) for saving
    // mean of prediction over effective samples for whole matrix
    const char* mean_file_path_;
    // The file path (including directory and file name) for saving
    // std (standard deviation) of prediction in gap filled data, whole matrix
    const char* std_file_path_;
    
    int verbose_;
    
    // A matrix containing the mean of prediction over effective samples
    // for whole matrix. At the end this will be saved into the file in
    // mean_file_path_.
    matVecFloat all_preds_mean_;
    // A matrix containing the std of prediction over effective samples
    // for whole matrix. At the end this will be saved into the file in
    // std_file_path_.
    matVecFloat all_preds_std_;
    
    matVecFloat test_data_;       //?????????????
    vector<float> test_preds_mean_; //???????????
    vector<float> test_preds_std_; //???????????
    int num_row_test_;    //??????????????????
    
    
    // Running gibbs sampling for the matrix at Level "level". Update both
    // U{level} and V{level}. First sample Rows of U in parallel, then using U
    // sample rows of V.
    void GibbsTrain(int level);
    
    // Initialize the latent matrix tree.
    // Where each node is linked to its parent and children.
    // if num_nodes!=0 then all levels have a fixed number of nodes as numNodes
    // For example for initializing V matrix num_nodes=num_columns_;
    void IntializeTree(NodeLatent*** tree, const int num_nodes);
    
    // Set the observed values for both matrices U and V at all levels.
    // It reads observed values from files located at input_dir_ as follows.
    // "Ytrain1.txt" contains the observed values for nodes at level 1.
    // "Ytrain2.txt" contains the observed values for nodes at level 1.
    void SetObservation();
    
    // Initialize Hierarchy informaiton. Read the information for hierarchy_,
    // num_parents_ and num_children_ from files "hierarchy_info.txt",
    // "num_parents.txt", and "num_children.txt", at input_dir, respectively.
    // If (used_hierarchy_level_ == predict_level_-1) &&
    //    (used_hierarchy_level_!=num_hierarchy_level_) updates parents of nodes
    // at lower levels to used_hierarchy_level_ and updates num_children at
    // used_num_hierarchy_level to the num_children at parent of predict_level
    void IntializeHierarchy();

    // Predicts the missing data for a node.
    // @u_ind: index of node at u matrix (row index)
    // @v_ind: index of node at v matrix (column index)
    // @level: at which level the node is located.
    // returns the predicted value as u[u_ind] * v_[v_ind]'
    double Predict(int u_ind, int v_ind, int level);
    
    double findTestErr();
    double findTestErr(FILE * outfile);
    double findTestErr(vector<double> &sum_test_preds, int iter);
    void FindMeanStdAllFields (int iter);
    void FindMeanStdTestData(int iter);
    
public:
    
    ~HPMF();
/*    HPMF(const char* input_dir, const int num_columns, int num_latent_features,
         int num_hierarchy_level, int used_hierarchy_level, int predict_level,
         int fold_id, int saveFileFlag, int outWholeFlag, int gaps, int burn,
         int num_samples, const char* mean_file_path, const char* std_file_path,
         double* num_nodes_per_level);
 */
HPMF(const int num_columns, int num_latent_features, double* num_nodes_per_level,
     int predict_level, int num_hierarchy_level, int used_hierarchy_level,
     int gaps, int burn, int num_samples, int fold_id, const char *input_dir,
     int saveFileFlag, int outWholeFlag, const char* mean_file_path,
     const char* std_file_path, int verbose);
    
    double runBHPMF(double *allRMSE, int &num_columns_);
};


inline int randWrapper(const int n) { return floor(unif_rand()*n); }

#endif /* defined(____HPMF__) */
