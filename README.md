## Bayesian Probabilistic Matrix Factorization

This R code provides an algorithm to fill gaps in large hierarchical databases. 
The method was originally developed for plant trait data but is applicable to any hierarchically structured numerical database.
Detailed instructions are given in the Vignette.

Please note, due to compatibility issues with the C compiler in newer versions of R, please use R 3.4.4 when using BHPMF.

For Mac users, please use the following instructions to make sure the C compiler works (thank you to Tom Walker for sorting this out):
1. Install Xcode and command line developer tools.
2. Install openMP: http://thecoatlessprofessor.com/programming/openmp-in-r-on-os-x/
3. Manually reinstall the header files (new since macOS Sierra and Mojave): https://donatstudios.com/MojaveMissingHeaderFiles

