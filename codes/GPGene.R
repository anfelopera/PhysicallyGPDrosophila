#### kernel functions for GP-Gene model ####

## Covariance function Kyy
kernY <- function(x1, x2, par) {
  # x1, x2 = matrices with  (time coordinate, space coordinate)
  # par = (variance, time lengthscale, space lengthscale)
  sigma2 <- par[1]; theta_t <- par[2]; theta_x <- par[3]
  
  # precomputing some terms
  dist2_t <- outer(x1[,1]/theta_t, x2[,1]/theta_t,'-')^2
  dist2_x <- outer(x1[,2]/theta_x, x2[,2]/theta_x,'-')^2
  
  # computing the SE kernel
  kern <- sigma2*exp(-dist2_t-dist2_x)
  
  # gradients
  attr(kern, "gradient") <- list(sigma2 = kern/sigma2,
                                 theta_t = 2*kern*dist2_t/theta_t,
                                 theta_x = 2*kern*dist2_x/theta_x)
  return(kern)
}

## Covariance function Kuu
kernU <- function(x1, x2, par) {
  # x1, x2 = matrices with  (time coordinate, space coordinate)
  # par = (sigma2, theta_t, theta_x, decay (lambda), diffusion (D), sensitivity (S))
  sigma2 <- par[1]; theta_t <- par[2]; theta_x <- par[3]
  lambda <- par[4]; D <- par[5]; S <- par[6]
  
  # precomputing some terms
  dist2_t <- outer(x1[,1]/theta_t, x2[,1]/theta_t,'-')^2
  dist2_x <- outer(x1[,2]/theta_x, x2[,2]/theta_x,'-')^2
  
  # computing the SE kernel
  kern0 <- sigma2*exp(-dist2_t-dist2_x)
  # Second derivative of the SE kernel wrt the time 
  kern2t <- -2/theta_t^2 + 4*dist2_t/theta_t^2
  # Second derivative of the SE kernel wrt the space 
  kern2x <- -2/theta_x^2 + 4*dist2_x/theta_x^2
  # Fourth derivative of the SE kernel wrt the space 
  kern4x <- 12/theta_x^4 - 48*dist2_x/theta_x^4 + 16*dist2_x^2/theta_x^4
  # compute Kuu
  kern <- 1/S^2*(D^2*kern4x - 2*D*lambda*kern2x - kern2t + lambda^2) * kern0
  
  # gradients
  Dkern2t_Dtheta_t <- 2/theta_t^3 - 4*dist2_t/theta_t^3
  Dkern2x_Dtheta_x <- 2/theta_x^3 - 4*dist2_x/theta_x^3
  Dkern4x_Dtheta_x <- 12/theta_x^5 + 6*dist2_x/theta_x^5 - 8*dist2_x^2/theta_x^5

  Dkern_Dtheta_t <- -1/S^2*Dkern2t_Dtheta_t * kern0 + dist2_t/theta_t*kern
  Dkern_Dtheta_x <- 1/S^2*(D^2*Dkern4x_Dtheta_x - 2*D*lambda*Dkern2x_Dtheta_x)*kern0 +
    dist2_x/theta_x*kern
  Dkern_Dlambda <- 2/S^2*(-D*kern2x + lambda) * kern0
  Dkern_DD <- 2/S^2*(D*kern4x - lambda*kern2x) * kern0
  attr(kern, "gradient") <- list(sigma2 = kern/sigma2,
                                 theta_t = 2*Dkern_Dtheta_t,
                                 theta_x = 2*Dkern_Dtheta_x,
                                 lambda = Dkern_Dlambda,
                                 D = Dkern_DD, S = -2*kern/S)
  return(kern)
}

## Cross-covariance function Kuy
kernUY <- function(x1, x2, par) {
  # x1, x2 = matrices with  (time coordinate, space coordinate)
  # par = (sigma2, theta_t, theta_x, decay (lambda), diffusion (D), sensitivity (S))
  sigma2 <- par[1]; theta_t <- par[2]; theta_x <- par[3]
  lambda <- par[4]; D <- par[5]; S <- par[6]
  
  # precomputing some terms
  dist_t  <- outer(x1[,1]/theta_t, x2[,1]/theta_t,'-')
  dist2_t <- dist_t^2
  dist2_x <- outer(x1[,2]/theta_x, x2[,2]/theta_x,'-')^2
  
  # computing the SE kernel
  kern0 <- sigma2*exp(-dist2_t-dist2_x)
  # Second derivative of the SE kernel wrt the space
  kern2x <- -2/theta_x^2 + 4*dist2_x/theta_x^2
  # First derivative of the SE kernel wrt the time 
  kern1t <- -2*dist_t/theta_t
  # compute Kuy
  kern <- 1/S*(-D*kern2x + kern1t + lambda) * kern0
  
  # gradients
  Dkern2x_Dtheta_x <- 2/theta_x^3 - 4*dist2_x/theta_x^3
  Dkern1t_Dtheta_t <- 2*dist_t/theta_t^2
  
  Dkern_Dtheta_t <- 1/S*Dkern1t_Dtheta_t * kern0 + dist2_t/theta_t*kern
  Dkern_Dtheta_x <- 1/S*(-D*Dkern2x_Dtheta_x) * kern0 + dist2_x/theta_x*kern
  attr(kern, "gradient") <- list(sigma2 = -kern/sigma2,
                                 theta_t = -2*Dkern_Dtheta_t,
                                 theta_x = -2*Dkern_Dtheta_x,
                                 lambda = -kern0/S,
                                 D = kern2x*kern0/S, S = kern/S)
  return(kern)
}

kernYU <- function(x1, x2, par) {
  # x1, x2 = matrices with  (time coordinate, space coordinate)
  # par = (sigma2, theta_t, theta_x, decay (lambda), diffusion (D), sensitivity (S))
  sigma2 <- par[1]; theta_t <- par[2]; theta_x <- par[3]
  lambda <- par[4]; D <- par[5]; S <- par[6]
  
  # precomputing some terms
  dist_t  <- outer(x1[,1]/theta_t, x2[,1]/theta_t,'-')
  dist2_t <- dist_t^2
  dist2_x <- outer(x1[,2]/theta_x, x2[,2]/theta_x,'-')^2
  
  # computing the SE kernel
  kern0 <- sigma2*exp(-dist2_t-dist2_x)

  # Second derivative of the SE kernel wrt the space
  kern2x <- -2/theta_x^2 + 4*dist2_x/theta_x^2
  # First derivative of the SE kernel wrt the time 
  kern1t <- 2*dist_t/theta_t 
  # compute Kyu
  kern <- 1/S*(-D*kern2x + kern1t + lambda) * kern0
  
  # gradients
  Dkern2x_Dtheta_x <- 2/theta_x^3 - 4*dist2_x/theta_x^3
  Dkern1t_Dtheta_t <- -2*dist_t/theta_t^2
  
  Dkern_Dtheta_t <- 1/S*Dkern1t_Dtheta_t * kern0 + dist2_t/theta_t*kern
  Dkern_Dtheta_x <- 1/S*(-D*Dkern2x_Dtheta_x) * kern0 + dist2_x/theta_x*kern
  attr(kern, "gradient") <- list(sigma2 = kern/sigma2,
                                 theta_t = 2*Dkern_Dtheta_t,
                                 theta_x = 2*Dkern_Dtheta_x,
                                 lambda = kern0/S,
                                 D = -kern2x*kern0/S, S = -kern/S)
  return(kern)
}

## Covariance function of the joint process
jointKern <- function(x1, x2, par) {
  # x1, x2 = matrices with  (time coordinate, space coordinate, label (0 for 'y', 1 for 'u'))
  # par = (sigma2, theta_t, theta_x, decay (lambda), diffusion (D), sensitivity (S))
  sigma2 <- par[1]; theta_t <- par[2]; theta_x <- par[3]
  lambda <- par[4]; D <- par[5]; S <- par[6]
  
  ### compute kernel (per block)
  indY1 = which(x1[,3]==0)
  indU1 = which(x1[,3]==1)
  indY2 = which(x2[,3]==0)
  indU2 = which(x2[,3]==1)
  
  kern <- matrix(0,nrow(x1),nrow(x2))
  Dkern_Dpar <- array(0,c(nrow(x1),nrow(x2),length(par)))
  
  # compute Kyy
  kernyy <- kernY(matrix(x1[indY1,1:2],ncol=2),matrix(x2[indY2,1:2],ncol=2),par[1:3])
  kern[indY1,indY2] <- kernyy
  Dkern_Dpar[indY1,indY2,1:3] <- array(unlist(attr(kernyy,'gradient')))

  # compute Kuu
  kernuu <- kernU(matrix(x1[indU1,1:2],ncol=2),matrix(x2[indU2,1:2],ncol=2),par)
  kern[indU1,indU2] <- kernuu
  Dkern_Dpar[indU1,indU2,] <- array(unlist(attr(kernuu,'gradient')))
  
  # compute cross-covariances
  kernyu <- kernYU(matrix(x1[indY1,1:2],ncol=2),matrix(x2[indU2,1:2],ncol=2),par)
  kernuy <- kernUY(matrix(x1[indU1,1:2],ncol=2),matrix(x2[indY2,1:2],ncol=2),par)
  kern[indY1,indU2] <- kernyu
  kern[indU1,indY2] <- kernuy
  Dkern_Dpar[indY1,indU2,] <- array(unlist(attr(kernyu,'gradient')))
  Dkern_Dpar[indU1,indY2,] <- array(-unlist(attr(kernuy,'gradient')))
  
  # gradients
  attr(kern, "gradient") <- list(sigma2 = Dkern_Dpar[,,1],
                                 theta_t = Dkern_Dpar[,,2],
                                 theta_x = Dkern_Dpar[,,3],
                                 lambda = Dkern_Dpar[,,4],
                                 D = Dkern_Dpar[,,5], S = Dkern_Dpar[,,6])
  return(kern)
}

GPGeneOpts <- function() {
  par <- c(sigma2 = 1, theta_t = 0.2, theta_x = 0.2,
              decay = 0.5, diffusion = 1e-4, sensitivity = 1)
  parLower <- c(sigma2 = 1e-10, theta_t = 0.01, theta_x = 0.01,
               decay = 1e-1, diffusion = 0, sensitivity = 0)
  parUpper <- c(sigma2 = Inf, theta_t = 2, theta_x = 2,
               decay = Inf, diffusion = 1e-3, sensitivity = Inf)
  parNames = c("sigma2", "theta_t", "theta_x", "decay", "diffusion", "sensitivity")
  return(list(par = par, parLower = parLower, parUpper = parUpper, parNames = parNames))
}
