//
//  latentNode.cpp
//
//
//  Created by Farideh Fazayeli on 7/15/13.
//
//

#include "latentNode.h"
#include <R.h>      // R functions
#include <Rmath.h>  // R math
#include <Rinternals.h>
#include <R_ext/Linpack.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>

NodeLatent::NodeLatent(int num_child, int num_feat, int level)
    : num_parent_(0),
      num_child_(num_child),
      last_child_(-1),
      num_latent_features_(num_feat),
      num_observ_(0) {
    if (num_child_) {
        //        prior_ = new double[num_latent_features_];
        children_ = new NodeLatent*[num_child_];
    }
    prior_ = std::unique_ptr<double[]>(new double[num_latent_features_]);
    
    //    latent_factor_ = new double[num_latent_features_];
    latent_factor_ = std::unique_ptr<double[]>(new double[num_latent_features_]);
    //fill_n(latent_factor_, num_latent_features_, 0);
    for(int i = 0; i < num_latent_features_; i++) {
        latent_factor_[i] = 2 * unif_rand() -1;
    }
    
    old_latent_factor_ = std::unique_ptr<double[]>(new double[num_latent_features_]);
    gradient_ = std::unique_ptr<double[]>(new double[num_latent_features_]);
    //    old_latent_factor_ = new double[num_latent_features_];
}

NodeLatent::NodeLatent(int num_parent, int num_child, int par_idx,
                       NodeLatent ***tree, int num_feat, int level,
                       int par_level)
    : num_parent_(num_parent),
      num_child_(num_child),
      last_child_(-1),
      num_latent_features_(num_feat),
      num_observ_(0) {
    if (num_child)
        children_ = new NodeLatent*[num_child];
    
    for (int ind = 0; ind < num_parent; ind++){
        parents_ = tree[par_level][par_idx-1];
        parents_->children_[++parents_->last_child_] = this;
    }
    
    if (num_parent_ || num_child_) {
        prior_ = std::unique_ptr<double[]>(new double[num_latent_features_]);
    }
    
    latent_factor_ = std::unique_ptr<double[]>(new double[num_latent_features_]);
    for(int i = 0; i < num_latent_features_; i++) {
        latent_factor_[i] = 2 * unif_rand() -1;
    }
    
    old_latent_factor_ =
          std::unique_ptr<double[]>(new double[num_latent_features_]);
    gradient_ = std::unique_ptr<double[]>(new double[num_latent_features_]);
}



NodeLatent::~NodeLatent() {
    SafeDelete(children_);
}

void NodeLatent::SetObserv(float observ, int ind) {
    observ_.push_back(observ);
    observ_idx_.push_back(ind);
    num_observ_++;
}


void NodeLatent::UpdatePrior(int mode) {
    const int inc_one = 1;
    const double cc = 1;
    
    std::unique_ptr<double[]> prior_parent(new double[num_latent_features_]);
    std::unique_ptr<double[]> prior_child(new double[num_latent_features_]);
    
    fill_n(prior_parent.get(), num_latent_features_, 0);
    fill_n(prior_child.get(), num_latent_features_, 0);
    fill_n(old_latent_factor_.get(), num_latent_features_, 0);
    
    for (int ind = 0; ind < num_parent_; ind++) {
        F77_NAME(daxpy)(&num_latent_features_, &cc, parents_->GetLatentFactor(),
                        &inc_one, prior_parent.get(), &inc_one);
    }
    for (int ind = 0; ind < num_child_; ind++) {
        F77_NAME(daxpy)(&num_latent_features_, &cc,
                        children_[ind]->GetLatentFactor(), &inc_one,
                        prior_child.get(), &inc_one);
    }
    
    if (mode == 1) {
        for (int ii = 0; ii < num_latent_features_; ii++) {
            latent_factor_[ii] = prior_parent[ii]/num_parent_ + rnorm(0, 0.01);
        }
    } else {
        for (int ii = 0; ii < num_latent_features_; ii++) {
            latent_factor_[ii] = prior_child[ii]/num_child_ + rnorm(0, 0.01);
        }
    }
    
    fill_n(prior_.get(), num_latent_features_, 0);
    F77_NAME(daxpy)(&num_latent_features_, &cc, prior_parent.get(), &inc_one,
                    prior_.get(), &inc_one);
    F77_NAME(daxpy)(&num_latent_features_, &cc, prior_child.get(), &inc_one,
                    prior_.get(), &inc_one);
}

void NodeLatent::UpdateGradient(int freq, double err, double *vv) {
    int inc_one = 1;
    if (freq == 0) {
        fill_n(gradient_.get(), num_latent_features_, 0);
    }
    // dE_du[i] = \sum_j^N_i err_ij * V_j
    F77_NAME(daxpy)(&num_latent_features_, &err, vv, &inc_one, gradient_.get(),
                    &inc_one);
}

void NodeLatent::UpdateLatentFactor(double coef, double momentum,
                                    double epsilon) {
    const int inc_one = 1;
    double cc = -1 * coef;
    
    // deV[j] -= coef * ( sum(v_pars) + sum(v_chs) )
    F77_NAME(daxpy)(&num_latent_features_, &cc, prior_.get(), &inc_one,
                    gradient_.get(), &inc_one);
    
    // numPrior = (num_parent_ + numCh)
    // coef = N_j * lambda * (num_parent_ + numCh)_v
    coef *= (num_parent_+num_child_);
    // deV[j] += coef * v[j]
    F77_NAME(daxpy)(&num_latent_features_, &coef, latent_factor_.get(), &inc_one,
                    gradient_.get(), &inc_one);
    
    F77_NAME(dscal)(&num_latent_features_, &momentum, old_latent_factor_.get(),
                    &inc_one);
    // old = old * momentum + epsilon * deV[j]
    F77_NAME(daxpy)(&num_latent_features_, &epsilon, gradient_.get(), &inc_one,
                    old_latent_factor_.get(), &inc_one);
    
    cc = -1;
    // latent_factor_ -= old;
    F77_NAME(daxpy)(&num_latent_features_, &cc, old_latent_factor_.get(), &inc_one,
                    latent_factor_.get(), &inc_one);
}

// Finding iverse of a matrix
// @mat: a doubble matrix for inversion of size dim*dim
// @dim: dimension of matrix
void InverseMat(double* mat, int dim)
{
    std::unique_ptr<int[]> ipiv(new int [dim+1]);
    int mat_size = dim*dim;
    std::unique_ptr<double[]> tmp(new double [mat_size]);
    int info;
    
    F77_NAME(dgetrf)(&dim, &dim, mat, &dim, ipiv.get(), &info);
    F77_NAME(dgetri)(&dim, mat, &dim, ipiv.get(), tmp.get(), &mat_size, &info);
}


void NodeLatent::UpdatePrior() {
    const double one = 1;
    const int inc_one = 1;
    std::unique_ptr<double[]> prior_parent(new double[num_latent_features_]);
    std::unique_ptr<double[]> prior_child(new double[num_latent_features_]);
    
    fill_n(prior_parent.get(), num_latent_features_, 0);
    fill_n(prior_child.get(), num_latent_features_, 0);
    fill_n(old_latent_factor_.get(), num_latent_features_, 0);
    fill_n(prior_.get(), num_latent_features_, 0);
    
    for (int ind = 0; ind < num_parent_; ind++) {
        F77_NAME(daxpy)(&num_latent_features_, &one, parents_->GetLatentFactor(),
                        &inc_one, prior_parent.get(), &inc_one);
    }

    for (int ind = 0; ind < num_child_; ind++) {
        F77_NAME(daxpy)(&num_latent_features_, &one,
                        children_[ind]->GetLatentFactor(),
                        &inc_one, prior_child.get(), &inc_one);
    }
    
    
    F77_NAME(daxpy)(&num_latent_features_, &one, prior_parent.get(), &inc_one,
                    prior_.get(), &inc_one);
    
    F77_NAME(daxpy)(&num_latent_features_, &one, prior_child.get(), &inc_one,
                    prior_.get(), &inc_one);
}

void NodeLatent::GibbsUpdate(NodeLatent*** v_tree, double sig_inv,
                             double sig_u_inv, int level){
    /************************************
     BLAS, Lapack Parameter
     *************************************/
    const char *lower = "L";
    const char *ntrans = "N";
    const char *ytrans = "T";
    // int num_latent_features_L = num_latent_features_;
    int info;
    
    const double zero = 0;
    const double one = 1;
    const int inc_one = 1;

    // Updates prior information for the node
    // Adds parent and children prior together
    UpdatePrior();
    
    /************************************
     Setup variables
     *************************************/
    double numPrior = (double) (num_child_ + num_parent_);
    
    std::unique_ptr<double[]> mu(new double[num_latent_features_]);
    std::unique_ptr<double[]> mu_u(new double[num_latent_features_]);
    std::unique_ptr<double[]> cov(new double[num_latent_features_ * num_latent_features_]);
    std::unique_ptr<double[]> latent_trans(new double[num_observ_ * num_latent_features_]);
    
    fill_n(mu.get(), num_latent_features_, 0);
    
    for (int ind = 0; ind < num_observ_; ind++) {
        double beta = observ_[ind] * sig_inv;
        double *tmp = v_tree[level][observ_idx_[ind]]->GetLatentFactor();
        
        copy(tmp, tmp+num_latent_features_, latent_trans.get()+ind*num_latent_features_);
        
        //mu += v_j * r_{ij} / sig
        F77_NAME(daxpy)(&num_latent_features_, &beta, tmp, &inc_one, mu.get(), &inc_one);
    }
    
    // mu += prior / sigU ; prior = sum(u_ch) + u_pr
    F77_NAME(daxpy)(&num_latent_features_, &sig_u_inv, prior_.get(), &inc_one,
                    mu.get(), &inc_one);
    
    // cov = Identity / sigU
    for(int ii = 0; ii < num_latent_features_*num_latent_features_; ii++)
        cov[ii] = 0.0;
    for(int ii = 0; ii < num_latent_features_; ii++)
        cov[ii*num_latent_features_+ii] = sig_u_inv;
    
    // cov = V_i' * V_i + (|ch| + 1)/sigUInv I
    // V_i' = latentT, V_i sub matrix of V for non zero values of i
    F77_NAME(dgemm)(ntrans, ytrans, &num_latent_features_, &num_latent_features_,
                    &num_observ_, &sig_inv, latent_trans.get(),
                    &num_latent_features_, latent_trans.get(),
                    &num_latent_features_, &numPrior, cov.get(),
                    &num_latent_features_);

    for(int ii = 0; ii < num_latent_features_; ii++) {
        cov[ii*num_latent_features_+ii] =
            cov[ii*num_latent_features_+ii] + 0.00001;
    }
    
    InverseMat(cov.get(), num_latent_features_);
    
    F77_NAME(dgemv)(ntrans, &num_latent_features_, &num_latent_features_, &one,
                    cov.get(), &num_latent_features_, mu.get(), &inc_one, &zero,
                    mu_u.get(), &inc_one);
    
    F77_NAME(dpotrf)(lower, &num_latent_features_, cov.get(),
                     &num_latent_features_, &info);
    
    mvrnorm(latent_factor_.get(), mu_u.get(), cov.get(), num_latent_features_,
            false);
    
}
