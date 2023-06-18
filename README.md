# aCGH_application_wavelets
This project focuses on the application of wavelet analysis and Markov Chain Monte Carlo (MCMC) methods to analyze an array comparative genomic hybridization (aCGH) data set. The goal is to detect aberrations in DNA copy number by utilizing advanced statistical techniques discussed in my thesis titled "Bayesian estimation of dynamic mixture models by wavelets".

Key Methods:
- Wavelet Regression (WR) Approach: The original aCGH data is transformed into a regression framework, where the regression function represents the dynamic probability to be estimated. This approach leverages wavelet analysis to model and analyze the data effectively.
- Data Augmentation (DA) Approach: This approach utilizes the data augmentation method proposed by Albert and Chibb (1993) to model the allocation data of the two-component mixture model. A discrete wavelet transform (DWT) matrix is employed as the design matrix in the Probit regression. Four different prior distributions are considered for the wavelet coefficients: a diffuse prior, a Gaussian prior, a spike and slab prior with Gaussian slab, and a spike and slab prior with Laplace slab.

By applying these methods, the project aims to identify DNA copy number aberrations and gain insights into the underlying genomic variations.

Overall, this project showcases the effectiveness of wavelet analysis and MCMC methods in analyzing aCGH data. The findings contribute to the field of genomic research, providing valuable insights into DNA copy number alterations and their implications.

