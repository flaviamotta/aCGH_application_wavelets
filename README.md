# aCGH_application_wavelets
Using Wavelets and MCMC methods to analyze an aCGH data set

Array comparative genomic hybridization (aCGH) is a technique that, under certain conditions,
provides intense fluorescence signals which are used for detecting aberrations in DNA copy number.
We use an aCGH data set to illustrate how the MCMC methods discussed in my dissertation ("Bayesian 
estimation of dynamic mixture models by wavelets") can be effective in this kind of application:

WR - Wavelet regression approach: Consists of transforming the original data into a regression,
whose regression function is the dynamic probability to be estimated. 

DA - Data augmentation approach: Uses the data augmentation method by Albert and Chibb (1993)
to model the allocation data of the two-component mixture model. We use a DWT matrix as the
design matrix in the Probit regression. We consider 4 priors to the distribution of the wavelet
coefficients:
1 - Diffuse;
2 - Gaussian;
3 - Spike and Slab (Gaussian slab);
4 - Spike and Slab (Laplace slab).


