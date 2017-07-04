\name{hierarchy.info}
\docType{data}

\alias{hierarchy.info}


\title{Taxonomy hierarchy information for plant traits}

\description{
  
  An example of hierarchy information: plant taxonomy information. The plant kingdom is characterized by a
taxonomic hierarchy. In the taxonomy, individual plants are categorized to species, species to genera, genera to families, and families to phylogenetic groups. Due to the data privacy, we perturbed a subset of taxonomy information from TRY database and not present the phylogenetic groups information.
  This data set represent the hierarchy as a taxa table as follow:
  the 1st column is the plant index in the same order as X data, the 2nd column holds species information,
  the 3rd column holds genus information, and the 4th column holds family information. Figure below represents a schematic of the taxanomy hierarchy structure:
  
\if{html}{\figure{taxa.png}{options: width="80\%" alt="Figure: taxa.png"}}
\if{latex}{\figure{taxa.png}{options: width=12cm}}
  
}

\usage{
data(hierarchy.info)
}
 
\format{
 A data frame containing 20,000 rows and 4 columns.
}

\references{
F. Schrodt, J. Kattge, H. Shan, F. Fazayeli, J. Joswig, A. Banerjee, M. Reichstein, M. Boenisch, S. Diaz, J. Dickie, A. Gillison, A. Karpatne, S. Lavorel, P.W. Leadley, C. Wirth, I. Wright, S.J. Wright, P.B. Reich (2015) BHPMF – a hierarchical Bayesian approach to gap-filling and trait prediction for macroecology and functional biogeography. \emph{Global Ecology and Biogeography} \bold{24}, 1510--1521

  J. Kattge, S. Diaz, S. Lavorel et al. (2011) Try-a global database of plant
traits. \emph{Global Change Biology}, \bold{17}, 2905--2935.
}

\keyword{datasets}
