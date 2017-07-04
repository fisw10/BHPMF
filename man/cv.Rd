\name{CalculateCvRmse}
\alias{CalculateCvRmse}
\title{Cross Validation RMSE for BHPMF}

\description{
		Calculate the average RMSE and Standard deviation using cross validation.
}

\usage{
CalculateCvRmse(X, hierarchy.info, num.folds=10, prediction.level,
        used.num.hierarchy.levels, num.samples=1000, burn=100, gaps=2,
        num.latent.feats=15, tuning=FALSE, num.folds.tuning, tmp.dir,
        verbose=FALSE)
}

\arguments{
  \item{X}{The data matrix (including missing values). Rows are the observations, columns are the features. Missing values must be filled by NA.}
  \item{hierarchy.info}{A matrix containing the hierarchical information in the same order as matrix X. Data should be ordered in descending order (the lowest level being in the last column). The first column is the observation index.}
  \item{num.folds}{Number of cross validation folds. The default value is 10.}
  \item{prediction.level}{The level at which gaps are filled. For example if set to 4, gaps are filled at the 3rd highest hierarchy level. The default value is at the observation level (filling gaps/missing values at the observation level).}
  \item{used.num.hierarchy.levels}{Number of hierarchy levels used for gap filling. For example, if set to 2, only information from the first and second levels of the hierarchy are used. The default value equals the total number of hierarchy levels.}
  \item{num.samples}{This is aparameter for Gibbs sampling: total number of generated samples at each fold using gibbs sampling. Note: this is not the effective number of samples. The default value is 1000.}
  \item{burn}{This is aparameter for Gibbs sampling: burn in period. Number of initially sampled parameters discarded. The default value is 200. Should be smaller than num.samples.}
  \item{gaps}{This is aparameter for Gibbs sampling: gap between sampled parameters kept. The default value is 2.}
  \item{num.latent.feats}{Size of latent vectors in BHPMF. If tuning is FALSE, set to value (default value: 10)}
  \item{tuning}{If set to \code{TRUE}, BHPMF is tuned to determine the best value of num.latent.feats. The default value is False.}
  \item{num.folds.tuning}{Number of cross validation folds for tuning. The default value is 10 if tuning set to \code{TRUE}.}
  \item{tmp.dir}{A temporary directory used to save preprocessing files.
        If not provided, a tmp directory in /R/tmp will be created.
    	If provided, each time a function from this package is called, the
  		preprocessing files saved in this directory will be used. Helps to avoid having to re-run preprocessing functions for the same dataset.
		WARNING: This directory should be empty for the first time call on a new dataset!}
  \item{verbose}{If \code{TRUE}, progress of the sampler details are printed to the screen. 
         Otherwise, nothing is printed to the screen. Default is \code{FALSE}}
}

\value{
Return a list with the following tags:
\item{avg.rmse}{Average of cross validation RMSE among all folds.}
\item{std.rmse}{Standard deviation cross validation RMSE among all folds.}
\item{best.number.latent.features}{The best value for latent vectors length if tuning is \code{TRUE}}
\item{min.rmse}{The minimum RMSE for the latent vectors of size best.number.latent.features if tuning is \code{TRUE}}

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

\seealso{\code{\link{GapFilling}}, \code{\link{TuneBhpmf}}}

\examples{
\dontrun{
    # Read the input data
    data(trait.info)  # class matrix
    data(hierarchy.info) # Read the hierarchy information; class data.frame

    # Calculate average RMSE with the default values
    out1 <- CalculateCvRmse(trait.info, hierarchy.info)
    avg.rmse <- out1$AvgRMSE
    std.rmse <- out1$StdRMSE

    # Calculate average RMSE using partial hierarch: using only the first 2 level of hierarchy
    out2 <- CalculateCvRmse(trait.info, hierarchy.info, used.num.hiearchy.levels=2,
    	                    num.samples=500, burn=20, gaps=3, tuning=TRUE)
    avg.rmse <- out2$AvgRMSE
    std.rmse <-	out2$StdRMSE

    out3 <- CalculateCvRmse(trait.info, hierarchy.info,
    	                    num.samples=100, burn=20, gaps=3, tuning=FALSE)
    avg.rmse <- out3$AvgRMSE
    std.rmse <-	out3$StdRMSE
    }
}

\keyword{model}