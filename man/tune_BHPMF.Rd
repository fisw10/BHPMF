\name{TuneBhpmf}
\alias{TuneBhpmf}
\title{Tune latent factor size in BHPMF.}

\description{
        Tune the latent factor size (k) in BHPMF using cross validation.
}

\usage{
TuneBhpmf(X, hierarchy.info, prediction.level, used.num.hierarchy.levels,
		  	 	num.folds=10, num.samples=1000, burn=100, gaps=2, tmp.dir,
		  	 	verbose=FALSE)
}

\arguments{
  \item{X}{The data matrix (including missing values). Rows are the observations, columns are the features. Missing values must be filled by NA.}
  \item{hierarchy.info}{A matrix containing the hierarchical information in the same order as matrix X. Data should be ordered in descending order (the lowest level being in the last column). The first column is the observation index.}
  \item{prediction.level}{The level at which gaps are filled. For example if set to 4, gaps are filled at the 3rd highest hierarchy level. The default value is at the observation level (filling gaps/missing values at the observation level).}
  \item{used.num.hierarchy.levels}{Number of hierarchy levels used for gap filling. For example, if set to 2, only information from the first and second levels of the hierarchy are used. The default value equals the total number of hierarchy levels.}
    \item{num.folds}{Number of cross validation folds for tuning. The default value is 10 if tuning set to \code{TRUE}.}
  \item{num.samples}{This is aparameter for Gibbs sampling: total number of generated samples at each fold using gibbs sampling. Note: this is not the effective number of samples. The default value is 1000.}
  \item{burn}{This is aparameter for Gibbs sampling: burn in period. Number of initially sampled parameters discarded. The default value is 200. Should be smaller than num.samples.}
  \item{gaps}{This is aparameter for Gibbs sampling: gap between sampled parameters kept. The default value is 2.}
  \item{tmp.dir}{A temporary directory used to save preprocessing files.
        If not provided, a tmp directory in /R/tmp will be created.
    	If provided, each time a function from this package is called, the
  		preprocessing files saved in this directory will be used. Helps to avoid having to re-run preprocessing functions for the same dataset.
		WARNING: This directory should be empty for the first time call on a new dataset!}
  \item{verbose}{If \code{TRUE}, progress of the sampler details are printed to the screen. 
         Otherwise, nothing is printed to the screen. Default is FALSE.}
}

\value{
A list with the following tags:
\item{best.number.latent.features}{The best value for latent vectors length}
\item{min.rmse}{The minimum RMSE for the latent vectors of size best.number.latent.features}
}

\references{
Ecological description:

F. Schrodt, J. Kattge, H. Shan, F. Fazayeli, J. Joswig, A. Banerjee, M. Reichstein, M. Boenisch, S. Diaz, J. Dickie, A. Gillison, A. Karpatne, S. Lavorel, P.W. Leadley, C. Wirth, I. Wright, S.J. Wright, P.B. Reich (2015) BHPMF – a hierarchical Bayesian approach to gap-filling and trait prediction for macroecology and functional biogeography. \emph{Global Ecology and Biogeography} \bold{24}, 1510--1521

Technical descriptions:

F. Fazayeli, A. Banerjee, J. Kattge, F. Schrodt, P.B. Reich (2014) Uncertainty quantified matrix completion using Bayesian Hierarchical Matrix factorization. \emph{International Conference on Machine Learning and Applications (ICMLA)}

H. Shan, J. Kattge, P. B. Reich, A. Banerjee, F. Schrodt, M. Reichstein. (2012) Gap Filling in the Plant Kingdom -- Trait Prediction Using Hierarchical Probabilistic Matrix Factorization. \emph{International Conference on Machine Learning (ICML)}
}

\author{
  Farideh Fazayeli \email{farideh@cs.umn.edu}
}

\seealso{\code{\link{GapFilling}}, \code{\link{CalculateCvRmse}}}

\examples{
\dontrun{
    library(BHPMF)
    data(trait.info)  # Read the matrix X
    data(hierarchy.info) # Read the hierarchy information
	
    # Tune the parameter with default values
    out1 <- TuneBhpmf(trait.info, hierarchy.info)
    num.latent.factors <- out1$best.number.latent.features
	
    # Tune the paramter using partial hierarchy: using only the first 2 level of hierarchy
    out2 <- TuneBhpmf(trait.info, hierarchy.info, used.num.hierarchy.levels=2,
                      num.samples=500, burn=20, gaps=3)
    num.latent.factors <- out2$best.number.latent.features	
}
}
\keyword{model}
