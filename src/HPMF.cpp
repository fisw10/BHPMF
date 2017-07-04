//
//  HPMF.cpp
//
//
//  Created by Farideh Fazayeli on 7/8/13.
//
//

#include "HPMF.h"


/*HPMF::HPMF(const char *input_dir, const int num_columns, int num_latent_features,
             int num_hierarchy_level, int used_hierarchy_level,
int predict_level, int fold_id, int saveFileFlag, int outWholeFlag,
int gaps, int burn, int num_samples, const char* mean_file_path,
const char* std_file_path, double* num_nodes_per_level)
*/
HPMF::HPMF(const int num_columns, int num_latent_features, double* num_nodes_per_level,
           int predict_level, int num_hierarchy_level, int used_hierarchy_level,
           int gaps, int burn, int num_samples, int fold_id, const char *input_dir,
           int saveFileFlag, int outWholeFlag, const char* mean_file_path,
           const char* std_file_path, int verbose)
  : num_columns_(num_columns),
    num_latent_features_(num_latent_features),
    predict_level_(predict_level),
    num_hierarchy_level_(num_hierarchy_level),
    used_hierarchy_level_(used_hierarchy_level),
    gaps_(gaps),
    burn_(burn),
    num_samples_(num_samples),
    fold_id_(fold_id),
    input_dir_(input_dir),
    saveFileFlag(saveFileFlag),
    outWholeFlag(outWholeFlag),
    mean_file_path_(mean_file_path),
    std_file_path_(std_file_path),
    verbose_(verbose){
  if(used_hierarchy_level >= predict_level) {
    Rprintf("Error: wrong number of parameters: used_hierarchy_level should be smaller than predict_level\n");
    // error();
  }
  if (verbose_) {
    Rprintf("Initializing BHPMF paramters...\n");
  }
  
  this->num_effective_samples_ = floor((num_samples_ - burn_) / gaps_) + 1;
  
  this->num_nodes_per_level_.assign(num_nodes_per_level,
                                    num_nodes_per_level+num_hierarchy_level+1);
  
  /*****************************
  Initializing the hierarchy
  *****************************/
  if (this->used_hierarchy_level_ > 0) {
    IntializeHierarchy();
  }
  
  /*****************************
  Initializing Latent matrices at all levels of hierarchy, represented as tree
  *****************************/
  nodePerLevel.resize(num_hierarchy_level_+1);
  u_tree_ = new NodeLatent**[used_hierarchy_level_+1];
  IntializeTree(u_tree_, 0);
  
  v_tree_ = new NodeLatent**[used_hierarchy_level_+1];
  IntializeTree(v_tree_, num_columns_);
  SetObservation();
  
  /*****************************
  Initializing all_pred
  *****************************/
  all_preds_mean_.resize(num_nodes_per_level_[predict_level_-1]);
  all_preds_std_.resize(num_nodes_per_level_[predict_level_-1]);
  for (int ii = 0; ii < num_nodes_per_level_[predict_level_-1]; ii++) {
    all_preds_mean_[ii].resize(num_columns_);
    all_preds_std_[ii].resize(num_columns_);
    for (int jj = 0; jj < num_columns_; jj++) {
      all_preds_mean_[ii][jj] = 0;
      all_preds_std_[ii][jj] = 0;
    }
  }
  
  
  /************************************
  Initialize test_pred
  *************************************/
 // string test_var = string(input_dir_) + "/fold" + to_string(static_cast<long double>(fold_id_)) 
 //   + "/Ytest" + to_string(static_cast<long double>(predict_level_)) +".txt";

  string test_var = string(input_dir_) + "/fold" + to_string(fold_id_)
    + "/Ytest" + to_string(predict_level_) +".txt";
	
  dlmreadVec(test_var, test_data_, num_row_test_, 3, verbose_);
  test_preds_std_.resize(num_row_test_);
  test_preds_mean_.resize(num_row_test_);
  
}

/*****************************
Free Heap variables
******************************/
HPMF::~HPMF() {
  for (int ii = used_hierarchy_level_; ii > 0; ii--) {
    for (int jj = 0; jj < num_columns_; jj++) {
      SafeDelete(v_tree_[ii][jj]);
    }
    for(int jj = 0; jj < num_nodes_per_level_[ii]; jj++) {
      SafeDelete(u_tree_[ii][jj]);
    }
    SafeDelete(u_tree_[ii]);
    SafeDelete(v_tree_[ii]);
  }
  SafeDelete(u_tree_);
  SafeDelete(v_tree_);
}

void HPMF::IntializeHierarchy() {
  
  int num_rows = accumulate(num_nodes_per_level_.begin(),
                            num_nodes_per_level_.end(), 0);
  
  // read num_parents info
  string file_path = string(input_dir_) + "/num_parents.txt";
  dlmreadVec(file_path, num_parents_);
  
  // read num_children info
  file_path = string(input_dir_) + "/num_children.txt";
  dlmreadVec(file_path, num_children_);
  
  // read the hierarchy info
  file_path = string(input_dir_) + "/hierarchy_info.txt";
  dlmreadVec(file_path, hierarchy_, num_rows, 1, verbose_); // num_parents_);
  
  // int last_ind = 0;
  if ((used_hierarchy_level_ == predict_level_-1) ||
      (used_hierarchy_level_ == num_hierarchy_level_)) return;
  
  int parent_level_start_ind = accumulate(num_nodes_per_level_.begin(),
                                          num_nodes_per_level_.begin()+used_hierarchy_level_, 0);
  
  int used_num_hierarchy_level_start_ind =
    parent_level_start_ind - num_nodes_per_level_[used_hierarchy_level_-1];
  
  vector<bool> flag_if_first_seen(num_rows, true);
  
  // updates parents of nodes at lower levels to used_hierarchy_level_
  for (int level = used_hierarchy_level_+1; level < predict_level_; level++) {
    int start_ind = parent_level_start_ind+num_nodes_per_level_[level-1];
    int last_ind = start_ind +num_nodes_per_level_[level];
    for (int node_ind = start_ind; node_ind < last_ind; node_ind++) {
      vector<int> new_parents_list;
      for (int parent_ind = 0; parent_ind < num_parents_[node_ind];
      parent_ind++) {
        int parent = hierarchy_[node_ind][parent_ind];
        vector<int> grant_parent =
          hierarchy_[parent_level_start_ind+parent-1];
        //point to parents of your parents
        new_parents_list.insert(new_parents_list.end(),
                                grant_parent.begin(), grant_parent.end());
      }
      hierarchy_[node_ind] = new_parents_list;
      num_parents_[node_ind] = new_parents_list.size();
    }
    parent_level_start_ind = start_ind;
  }
  
  parent_level_start_ind = accumulate(num_nodes_per_level_.begin(),
                                      num_nodes_per_level_.begin()+predict_level_-2, 0);
  
  // updates num_children at used_num_hierarchy_level to the num_children
  // at parent of predict_level
  for (int node_ind = parent_level_start_ind;
       node_ind < parent_level_start_ind+num_nodes_per_level_[predict_level_-2];
       node_ind++) {
    for (int parent_ind = 0; parent_ind < num_parents_[node_ind];
    parent_ind++) {
      int parent = hierarchy_[node_ind][parent_ind];
      int grant_parent_ind =
        used_num_hierarchy_level_start_ind + parent - 1;
      if (flag_if_first_seen[grant_parent_ind]) {
        num_children_[grant_parent_ind] = num_children_[node_ind];
        flag_if_first_seen[grant_parent_ind] = false;
      }
      else {
        num_children_[grant_parent_ind] += num_children_[node_ind];
      }
    }
  }
}

void HPMF::IntializeTree(NodeLatent*** tree, const int num_nodes) {
  int last_level_ind = 0;
  for (int level = 0; level <=  used_hierarchy_level_; level++) {
    int lev;
    if (level == used_hierarchy_level_){
      for (int ii = level; ii < predict_level_-1; ii++) {
        if (!num_nodes)
          last_level_ind += num_nodes_per_level_[ii];
      }
      lev = predict_level_-1;  // to find number of nodes in the level
    } else {
      lev = level;
    }
    
    int parLevel = level - 1;
    
    int num_rows;
    // if num_nodes is not zero, have the same number of nodes per level.
    if (!num_nodes) { // treeU
      num_rows = num_nodes_per_level_[lev];
      nodePerLevel[level] = num_rows;
    } else { // treeV
      num_rows = num_nodes;
    }
    
    tree[level] = (NodeLatent **) R_alloc(num_rows, sizeof(NodeLatent*));
    for (int node = 0; node < num_rows; node++) {
      int num_children;
      int num_parents;
      vector<int> parents;
      vector <int> empty_parent (1,0);
      if (num_nodes) {  // treeV
        num_children = (level == used_hierarchy_level_) ?  0 : 1;
        num_parents = 1;
        parents.push_back(node+1);
      } else { //treeU
        num_children = (level == used_hierarchy_level_)
        ?  0 : num_children_[node+last_level_ind];
        num_parents = (used_hierarchy_level_ == 0)
          ?  0 : num_parents_[node+last_level_ind];
        parents = (used_hierarchy_level_ == 0)
          ?  empty_parent : hierarchy_[node+last_level_ind];
      }
      
      if (!level){  // if it is at the root of tree with no parent
        tree[level][node] =
          new NodeLatent(num_children, num_latent_features_, level);
      } else {
        tree[level][node] = new NodeLatent(num_parents,
                                           num_children,
                                           parents[0],
                                                  tree,
                                                  num_latent_features_,
                                                  level,
                                                  parLevel);
      }
    }
    if (!num_nodes) last_level_ind += num_nodes_per_level_[level];
  }
}

void HPMF::SetObservation() {
  for (int level = 0; level <= used_hierarchy_level_; level++) {
    int lev;
    if (level == used_hierarchy_level_) {
      lev = predict_level_-1;
    } else {
      lev = level;
    }
    
    /************************************
    reading the training data
    check for cross validation or learning
    *************************************/
    int tmp_level = lev+1;
	// string train_var = string(input_dir_) + "/fold" + to_string(static_cast<long double>(fold_id_)) 
    //  + "/Ytrain" + to_string(static_cast<long double>(tmp_level)) +".txt";
	  
    string train_var = string(input_dir_) + "/fold" + to_string(fold_id_)
      + "/Ytrain" + to_string(tmp_level) +".txt";
    
    int num_rows;
    matVecFloat train_data;
    dlmreadVec(train_var, train_data, num_rows, 3, verbose_);
	
    int len_training_data = num_rows;
    
    for (int ind = 0; ind < len_training_data; ind++) {
      int ii = (int) train_data[ind][0]-1; // row index
      int jj = (int) train_data[ind][1]-1; // column index
      float obs = train_data[ind] [2];     // observed value
      
      u_tree_[level][ii]->SetObserv(obs, jj);
      v_tree_[level][jj]->SetObserv(obs, ii);
    }
  }
}

double HPMF::runBHPMF(double *allRMSE, int &noPlants) {
  noPlants =  num_nodes_per_level_[predict_level_-1];
  
  // string test_var = string(input_dir_) + "/fold" + to_string(static_cast<long double>(fold_id_))
  //  + "/Ytest" + to_string(static_cast<long double>(predict_level_)) +".txt";
  
  string test_var = string(input_dir_) + "/fold" + to_string(fold_id_)
    + "/Ytest" + to_string(predict_level_) +".txt";
  
  int len_test_data;
  matVecFloat test_data;
  dlmreadVec(test_var, test_data, len_test_data, 3, verbose_);
  
  vector<double> sum_test_preds(len_test_data, 0);
  int count = 1;
  numPlants = nodePerLevel[predict_level_-1];
  for(int iter = 0; iter < num_samples_; iter++) {
    double t1 = omp_get_wtime();
    if (verbose_) {
      Rprintf("iter:  %i \n", iter);
    }
    if (used_hierarchy_level_ == 0) {
      int level = predict_level_ - 1;
      if (verbose_) {
        Rprintf("level %i * trait \n", level+1);
      }
      numPlants = num_nodes_per_level_[level];
      GibbsTrain(0);
    }
    
    if (verbose_) {
      Rprintf("------TopDown----- \n");
    }
    // for each level from down to top
    for (int indLevel = 1; indLevel <= used_hierarchy_level_ ; indLevel++ ) {
      int level;
      // at the lowest level
      if ( indLevel == used_hierarchy_level_ || used_hierarchy_level_ == 1) {
        level = predict_level_-1;
      } else {
        level = indLevel;
      }
      
      if (verbose_) {
        Rprintf("level %i * trait \n", level+1);
      }
      numPlants = num_nodes_per_level_[level];
      GibbsTrain(indLevel);
    }
    
    if (verbose_) {
      Rprintf("------DownTop----- \n");
    }
    // for each level from down to top
    for (int level = used_hierarchy_level_-1; level >= 0 ; level-- ) {
      if (verbose_) {
        Rprintf("level %i * trait \n", level+1);
      }
      numPlants = num_nodes_per_level_[level];
      GibbsTrain(level);
    }
    
    if ((iter >= burn_) & ((iter-burn_) % gaps_ == 0)) {
      if (outWholeFlag) {
        FindMeanStdAllFields(count);
      }
      FindMeanStdTestData(count);
      
      count++;
    }
    double err = findTestErr();
    if (verbose_) {
      Rprintf("Test Err: %f \n", err);
    }
    allRMSE[iter] = err;
    double t2 = omp_get_wtime();
    if (verbose_) {
      Rprintf("it tooks %f seconds \n",  t2 - t1);
    }
  }
  
  // calculate the average
  if (outWholeFlag) {
    FILE* outMeanFile;
    FILE* outStdFile;
    outMeanFile = fopen(mean_file_path_, "w");
    outStdFile = fopen(std_file_path_, "w");
    for (int ii = 0; ii < num_nodes_per_level_[predict_level_-1]; ii++) {
      for (int jj = 0; jj < num_columns_; jj++) {
        all_preds_std_[ii][jj] /= (count-1);
        fprintf(outMeanFile,"%5.4f\t", all_preds_mean_[ii][jj]);
        fprintf(outStdFile,"%5.4f\t", sqrt(all_preds_std_[ii][jj]));
      }
      fprintf(outMeanFile,"\n");
      fprintf(outStdFile,"\n");
    }
    fclose(outMeanFile);
    fclose(outStdFile);
  }
  double sum_res = 0;
  if (saveFileFlag) {
    FILE* outMeanFileTest;
    FILE* outStdFileTest;
    outMeanFileTest = fopen(mean_file_path_, "w");
    outStdFileTest = fopen(std_file_path_, "w");
    for (int ind = 0; ind < num_row_test_;  ind++) {
      double res = test_preds_mean_[ind] - test_data_[ind][2];
      sum_res += pow(res, 2);
      test_preds_std_[ind] /= (count-1);
      fprintf(outMeanFileTest,"%5.4f\t", test_preds_mean_[ind]);
      fprintf(outStdFileTest,"%5.4f\t", sqrt(test_preds_std_[ind]));
    }
    fprintf(outMeanFileTest,"\n");
    fprintf(outStdFileTest,"\n");
    fclose(outMeanFileTest);
    fclose(outStdFileTest);
  } else {
    for (int ind = 0; ind < num_row_test_;  ind++) {
      double res = test_preds_mean_[ind] - test_data_[ind][2];
      sum_res += pow(res, 2);
      test_preds_std_[ind] /= (count-1);
    }
  }
  return sqrt(sum_res/num_row_test_);
}

double HPMF::Predict(int u_ind, int v_ind, int level) {
  const int inc_one = 1;
  double *uu = u_tree_[level][u_ind]->GetLatentFactor();  // uMat[u_ind]->u; 
  double *vv = v_tree_[level][v_ind]->GetLatentFactor();  // vMat[v_ind]->v;

  // pred = uu .* vv' (inner product)
  double pred = F77_NAME(ddot)(&num_latent_features_,
                         uu, &inc_one, vv, &inc_one);
  return pred;
}

void HPMF::GibbsTrain(int level) {
  
  /************************************
  Initial Parameter
  *************************************/
  double sig_inv = 1/0.01;
  double sig_u_inv = 1/0.01;
  double sig_v_inv = 1/0.01;
  
  // sample each row of U
#ifdef _OPENMP
#pragma omp parallel for num_threads(HPMF_num_threads)
#endif
  for (int node = 0; node < numPlants; node++) {
    u_tree_[level][node]->GibbsUpdate(v_tree_, sig_inv, sig_u_inv, level);
  }
  
  // sample each row of V
#ifdef _OPENMP
#pragma omp parallel for num_threads(HPMF_num_threads)
#endif
  for (int node = 0; node < num_columns_; node++){
    v_tree_[level][node]->GibbsUpdate(u_tree_, sig_inv, sig_v_inv, level);
  }
}

double HPMF::findTestErr() {
  int level = used_hierarchy_level_;
  
  // string test_var = string(input_dir_) + "/fold" + to_string(static_cast<long double>(fold_id_)) 
  //  + "/Ytest" + to_string(static_cast<long double>(predict_level_)) +".txt";
  
  string test_var = string(input_dir_) + "/fold" + to_string(fold_id_)
    + "/Ytest" + to_string(predict_level_) +".txt";
  //    Rprintf("test_var: %s\n", test_var);
  
  int len_test_data, nCol = 3;
  matVecFloat testData;
  dlmreadVec(test_var, testData, len_test_data, nCol, verbose_);
  
  /************************************
  Find Test error
  *************************************/
  double sumRes = 0;
  for (int ind = 0; ind < len_test_data; ind++) {
    double res;
    int ii  = (int) testData[ind][0]-1;
    int jj  = (int) testData[ind][1]-1;
    double pred = Predict(ii, jj, level);
    float obs = testData[ind][2];
    res = pred - obs;
    sumRes += pow(res, 2);
  }
  
  sumRes /= len_test_data;
  double testErr = sqrt(sumRes);
  return testErr;
}

double HPMF::findTestErr(vector<double> &sum_test_preds, int iter) {
  int level = used_hierarchy_level_;
  
  // string test_var = string(input_dir_) + "/fold" + to_string(static_cast<long double>(fold_id_)) 
  //  + "/Ytest" + to_string(static_cast<long double>(predict_level_)) +".txt";
  
  string test_var = string(input_dir_) + "/fold" + to_string(fold_id_)
    + "/Ytest" + to_string(predict_level_) +".txt";
  
  int len_test_data, nCol = 3;
  matVecFloat testData;
  dlmreadVec(test_var, testData, len_test_data, nCol, verbose_);
  
  /************************************
  Find validation error, early stop if ove
  *************************************/
  double sum_res_sq_per_iter = 0;
  
  for (int ind = 0; ind < len_test_data;  ind++) {
    int ii  = (int) testData[ind][0]-1;
    int jj  = (int) testData[ind][1]-1;
    double pred = Predict(ii, jj, level);
    float obs = testData[ind][2];
    double res = pred - obs;
    sum_res_sq_per_iter += pow(res, 2);
    sum_test_preds[ind] += pred;
    
  }
  
  if (iter == num_samples_-1) {
    double sum_res_sq_all = 0;
    for (int ind = 0; ind < len_test_data;  ind++) {
      double avg_pred = sum_test_preds[ind] / num_effective_samples_;
      double res = testData[ind][2] - avg_pred;
      sum_res_sq_all += pow(res, 2);
    }
    return sqrt(sum_res_sq_all/len_test_data);
  }
  
  return sqrt(sum_res_sq_per_iter/len_test_data);
}

double HPMF::findTestErr(FILE *outFile) {
  int level = used_hierarchy_level_;
  
  // string test_var = string(input_dir_) + "/fold" + to_string(static_cast<long double>(fold_id_)) 
  //  + "/Ytest" + to_string(static_cast<long double>(predict_level_)) +".txt";
  
  string test_var = string(input_dir_) + "/fold" + to_string(fold_id_)
    + "/Ytest" + to_string(predict_level_) +".txt";
  
  int len_test_data, nCol = 3;
  matVecFloat testData;
  dlmreadVec(test_var, testData, len_test_data, nCol, verbose_);
  dlmreadVec(test_var, testData, len_test_data, nCol, verbose_);
  
  /************************************
  Find validation error, early stop if over
  *************************************/
  double sumRes = 0;
  
  //  double test_preds[len_test_data];
  
  for (int ind = 0; ind < len_test_data;  ind++) {
    int ii  = (int) testData[ind][0]-1;
    int jj  = (int) testData[ind][1]-1;
    double pred = Predict(ii, jj, level);
    float obs = testData[ind][2];
    // test_preds[ind] = pred;
    
    if(outWholeFlag){
      fprintf(outFile,"%5.4f\t", pred);
    }
    double res = pred - obs;
    sumRes += pow(res, 2);
  }
  
  fprintf(outFile,"\n");
  sumRes /= len_test_data;
  return sqrt(sumRes);
}


void HPMF::FindMeanStdTestData(int iter) {
  int level = used_hierarchy_level_;
  
#ifdef _OPENMP
#pragma omp parallel for num_threads(HPMF_num_threads)
#endif
  for (int ind = 0; ind < num_row_test_;  ind++) {
    int ii  = (int) test_data_[ind][0]-1;
    int jj  = (int) test_data_[ind][1]-1;
    double pred = Predict(ii, jj, level);
    if (iter == 1) {
      test_preds_mean_[ind] = pred;
      test_preds_std_[ind] = 0;
    } else {
      float prev_mean = test_preds_mean_[ind];
      test_preds_mean_[ind] += (pred - prev_mean)/iter;
      test_preds_std_[ind] +=
        (pred - prev_mean) * (pred - test_preds_mean_[ind]);
    }
  }
}

void HPMF::FindMeanStdAllFields(int iter) {
  int level = used_hierarchy_level_;
  int num_nodes = num_nodes_per_level_[predict_level_-1]; 

#ifdef _OPENMP
#pragma omp parallel for num_threads(HPMF_num_threads)
#endif
  for (int ii = 0; ii < num_nodes; ii++) {
    for(int jj = 0; jj < num_columns_; jj++) {
      double pred = Predict(ii, jj, level);
      if (iter == 1) {
        all_preds_mean_[ii][jj] = pred;
        all_preds_std_[ii][jj] = 0;
      } else {
        float prev_mean = all_preds_mean_[ii][jj];
        all_preds_mean_[ii][jj] += (pred - prev_mean)/iter;
        all_preds_std_[ii][jj] +=
          (pred - prev_mean) * (pred - all_preds_mean_[ii][jj]);
      }
    }
  }
}
