\name{BHPMF-package}
\alias{BHPMF-package}
\alias{BHPMF-package}
\docType{package}
\title{
An R Package to fill missing values in a data matrix.
}
\description{
An R package to filling gaps in a data matrix with hierarchical information (e.g. taxonomy) using Gibbs sampling (BHPMF).
}
\details{
\tabular{rl}{
Package: \tab BHPMF \cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-11-23\cr
License: \tab GPL-2\cr
}

Inputs to the GapFilling function are a matrix X containing the missing values (e.g., plant trait information), and a matrix containing the hierarchy information. Output is a gap filled matrix (mean of the predicted distribution) in the same format as the matrix X along with another matrix containing the standard deviation informaiton in the same format as X, both are saved into given file paths. The standard deviation will be used for uncertainity quantification. Two parameters "used.num.hierarchy.levels" and "prediction.level" control the usage of hierarchy information and the level to predict missing data. Below, we explain how to use these two parameters effectively.

Figure 1 shows the inputs and output of the GapFilling function. Here the hierarchy has 4 levels where the observation level (e.g. plants) is the 4th level. 
In this example, the GapFilling function will predict the missing values at the observation level i.e., filling the missing values in matrix X by setting the parameter "prediction.level" to 4 (observation (e.g. species) level - default value).
Also, here we use all the hierarchy information (information in levels 1, 2 and 3, e.g. phylogenetic group, genus and family) by setting the parameter "used.num.hierarchy.levels" to 3 (the default value). The usage is shown in the examples.

\if{html}{\figure{Figure1.jpg}{options: width="35\%" alt="Figure: Figure1.jpg"}}
\if{latex}{\figure{Figure1.jpg}{options: width=7cm}}


Figure 2 shows the inputs and output of the GapFilling function where we use only part of the hierarchy information (information in levels 1, 2 (e.g. using only information from the phylogenetic group and family)) by setting the parameter "used.num.hierarchy.levels" to 2.
In this example, the GapFilling function will predict the missing values at the observation level (e.g. species) i.e., filling the missing values in matrix X since we set the parameter "prediction.level" to 4 (default value). The usage is shown in the examples.

\if{html}{\figure{Figure2.jpg}{options: width="35\%" alt="Figure: Figure2.jpg"}}
\if{latex}{\figure{Figure2.jpg}{options: width=7cm}}


Figure 3 shows the inputs and output of the GapFilling function where we use part of the hierarchy information (information in levels 1, 2) by setting the parameter "used.num.hierarchy.levels" to 2.
In this example, the GapFilling function will predict the missing values at level 3 i.e., filling the missing values in matrix at level 3 (smaller dimension, e.g. at the genus level) by setting the parameter "prediction.level" to 3. The usage is shown in the examples.

\if{html}{\figure{Figure3.jpg}{options: width="35\%" alt="Figure: Figure3.jpg"}}
\if{latex}{\figure{Figure3.jpg}{options: width=7cm}}

The paramter "used.num.hierarchy.levels" is useful to generate results similart to Figure 4 in F. Schrodt et al. 2015 and Figure 4 and Table 2 in H. Shan et al. 2012.

}

\author{
  Farideh Fazayeli, Arindam Banerjee, Franziska Schrodt, Jens Kattge,  Peter Reich

  We would like to gratefully acknowledge support from National Science Foundation. The research was supported by NSF grant IIS 1563950 (https://nsf.gov/awardsearch/showAward?AWD_ID=1563950).
}

\references{
Ecological description:

F. Schrodt, J. Kattge, H. Shan, F. Fazayeli, J. Joswig, A. Banerjee, M. Reichstein, M. Boenisch, S. Diaz, J. Dickie, A. Gillison, A. Karpatne, S. Lavorel, P.W. Leadley, C. Wirth, I. Wright, S.J. Wright, P.B. Reich (2015) BHPMF – a hierarchical Bayesian approach to gap-filling and trait prediction for macroecology and functional biogeography. \emph{Global Ecology and Biogeography} \bold{24}, 1510--1521

Technical descriptions:

F. Fazayeli, A. Banerjee, J. Kattge, F. Schrodt, P.B. Reich (2014) Uncertainty quantified matrix completion using Bayesian Hierarchical Matrix factorization. \emph{International Conference on Machine Learning and Applications (ICMLA)}

H. Shan, J. Kattge, P. B. Reich, A. Banerjee, F. Schrodt, M. Reichstein. (2012) Gap Filling in the Plant Kingdom -- Trait Prediction Using Hierarchical Probabilistic Matrix Factorization. \emph{International Conference on Machine Learning (ICML)}
}

\seealso{\code{\link{GapFilling}}, \code{\link{CalculateCvRmse}}, \code{\link{TuneBhpmf}}}

\keyword{ package }

\examples{
\dontrun{
    # Read the input data
    library(BHPMF)
    data(trait.info)  # Read the matrix X
    data(hierarchy.info) # Read the hierarchy information
    
    ###################
    # Usage 1: Figure 1
    ###################
    GapFilling(trait.info, hierarchy.info,
  	          mean.gap.filled.output.path = "/tmp/mean_gap_filled.txt",
    	       	std.gap.filled.output.path="/tmp/std_gap_filled.txt")
    
    ###################
    # Usage 2: Figure 2
    ###################
    GapFilling(trait.info, hierarchy.info,
              prediction.level = 4,
              used.num.hierarchy.levels = 2,
  	          mean.gap.filled.output.path = "/tmp/mean_gap_filled.txt",
    	       	std.gap.filled.output.path="/tmp/std_gap_filled.txt")
    
    ###################
    # Usage 3: Figure 3
    ###################
    GapFilling(trait.info, hierarchy.info,
              prediction.level = 3,
              used.num.hierarchy.levels = 2,
  	          mean.gap.filled.output.path = "/tmp/mean_gap_filled.txt",
    	       	std.gap.filled.output.path="/tmp/std_gap_filled.txt")  
    	       	
    ######################
    # Usage 4: Running PMF
    # not using hierarchy
    #####################
    GapFilling(trait.info, hierarchy.info,
              used.num.hierarchy.levels = 0,
  	          mean.gap.filled.output.path = "/tmp/mean_gap_filled.txt",
    	       	std.gap.filled.output.path="/tmp/std_gap_filled.txt")  
	

    ##########################################
    # Usage 4: Calculate cross validation RMSE
    ##########################################
    # Calculate average RMSE with the default values
    out1 <- CalculateCvRmse(trait.info, hierarchy.info)
    avg.rmse <- out1$avg.rmse
    std.rmse <- out1$std.rmse

    # Calculate average RMSE using partial hierarch: using only the first 2 level of hierarchy
    out2 <- CalculateCvRmse(X, hierarchy.info, used.num.hiearchy.levels=2,
    	                    num.samples=500, burn=20, gaps=3, tuning=TRUE)
    avg.rmse <- out2$avg.rmse
    std.rmse <-	out2$std.rmse
}
}

