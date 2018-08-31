# PhysicallyGPDrosophila

This repository contains the R implementation of two types of physically-inspired GPs based on the reaction-diffusion equation [1]. The main difference between both approaches lies on whether the GP prior is placed: either over mRNA expressions or protein concentrations. Both GP models are tested under different conditions depending on the availability of biological data.

**Authors:** Andrés Felipe López-Lopera (Mines Saint-Étienne) with contributions from Nicolas Durrande (Mines Saint-Étienne, Prowler.io) and Mauricio Álvarez (The University of Sheffield).

### Main dependences (R packages):
1. kergp: Gaussian process models with customised covariance kernels [4].
2. NORMT3: evaluates complex erf, erfc, Faddeeva, and density of sum of Gaussian and Student’s t [5].
3. Pracma: practical numerical math functions [6].

### References:

[1] A. F. López-Lopera, N. Durrande and M. A. Álvarez. 2018. Physically-inspired Gaussian processes for transcriptional regulation in Drosophila melanogaster. ArXiv E-Prints. [[link]](https://arxiv.org/abs/1808.10026)


[2] M. A. Álvarez, D. Luengo, and N. D. Lawrence. 2013. Linear Latent Force Models Using Gaussian Processes. Pattern Analysis and Machine Intelligence, IEEE Transactions on, 35(11):2693–2705.

[3] J. D. Vasquez Jaramillo, M. A. Álvarez and A. A. Orozco (2014). Latent force models for describing transcriptional regulation processes in the embryo development problem for the Drosophila melanogaster. In Engineering in Medicine and Biology Society (EMBC), 36th Annual International Conference of the IEEE, pages 338–341.

[4] Y. Deville, D. Ginsbourger, and O Roustant (2015). kergp: Gaussian process models with customised covariance kernels. (https://cran.r-project.org/web/packages/kergp/index.html). [Online; 23-Dec-2015].

[5] G. Nason (2012). NORMT3: evaluates complex erf, erfc, Faddeeva, and density of sum of Gaussian and Student's t. (https://cran.r-project.org/web/packages/NORMT3/index.html). [Online; 31-Oct-2012].

[6] H. W. Borchers (2012). Pracma: practical numerical math functions. (https://cran.r-project.org/web/packages/pracma/index.html). [Online; 30-Jan-2018].
