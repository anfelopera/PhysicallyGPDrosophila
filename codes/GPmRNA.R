#### kernel functions for GP-mRNA model ####

## Covariance function Kuu
kernU <- function(x1, x2, par) {
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

## Covariance function Kyy
kernY <- function(x1, x2, par) {
  # x1, x2 = matrices with  (time coordinate, space coordinate)
  # par = (sigma2, theta_t, theta_x, decay (lambda), diffusion (D), sensitivity (S))
  sigma2 <- par[1]; theta_t <- par[2]; theta_x <- par[3]
  lambda <- par[4]; D <- par[5]; S <- par[6]
  
  # fixing and precomputing some terms
  nterms <- 10
  l <- 1
  kern <- 0
  Dkern_Dtheta_x <- Dkern_Dtheta_t <- Dkern_Dlambda <- Dkern_DD <- 0
  
  # computing kernel and its gradient
  for(i in 1:nterms){
    for(j in 1:nterms){
      if((i+j)%%2 == 0){
        # computing kernel
        kernt <- simXsimKernCompute(x1[,1], x2[,1], i, j, par, l)
        kernx <- sheatXsheatKernCompute(x1[,2], x2[,2], i, j, par, l)
        kern <- kern + kernt*kernx
        
        # computing gradient
        Dkernt <- simXsimKernGradient(x1[,1], x2[,1], i, j, par, l)  
        Dkernx_theta_x <- sheatXsheatKernGradient(x1[,2], x2[,2], i, j, par, l)
        Dkern_Dtheta_x <- Dkern_Dtheta_x + kernt*Dkernx_theta_x
        Dkern_Dtheta_t <- Dkern_Dtheta_t + Dkernt$dkern_Dtheta_t*kernx
        Dkern_Dlambda <- Dkern_Dlambda + Dkernt$dkern_Dlambda*kernx
        Dkern_DD <- Dkern_DD + Dkernt$dkern_DD*kernx
      }
    }
  }
  kern <- sigma2*((2*S/l)^2)*kern
  
  # gradients
  Dkern_DS <- 2*kern/S
  Dkern_Dtheta_x <- sigma2*((2*S/l)^2)*Dkern_Dtheta_x
  Dkern_Dtheta_t <- sigma2*((2*S/l)^2)*Dkern_Dtheta_t
  Dkern_Dlambda <- sigma2*((2*S/l)^2)*Dkern_Dlambda
  Dkern_DD <- sigma2*((2*S/l)^2)*Dkern_DD
  attr(kern, "gradient") <- list(sigma2 = kern/sigma2,
                                 theta_t = Dkern_Dtheta_t,
                                 theta_x = Dkern_Dtheta_x,
                                 lambda = Dkern_Dlambda,
                                 D = Dkern_DD, S = Dkern_DS)
  return(kern)
}

## Cross-covariance function Kyu
kernYU <- function(x1, x2, par) {
  # x1, x2 = matrices with  (time coordinate, space coordinate)
  # par = (sigma2, theta_t, theta_x, decay (lambda), diffusion (D), sensitivity (S))
  sigma2 <- par[1]; theta_t <- par[2]; theta_x <- par[3]
  lambda <- par[4]; D <- par[5]; S <- par[6]
  
  # fixing and precomputing some terms
  nterms <- 10
  l <- 1
  kern <- 0
  Dkern_Dtheta_x <- Dkern_Dtheta_t <- Dkern_Dlambda <- Dkern_DD <- 0
  
  # computing kernel and its gradient
  for(i in 1:nterms){
    # computing kernel
    kernt <- simXrbfKernCompute(x1[,1], x2[,1], i, par, l)
    kernx <- sheatXrbfKernCompute(x1[,2], x2[,2], i, par, l)
    kern <- kern + kernt*kernx
    
    # computing gradient
    Dkernt <- simXrbfKernGradient(x1[,1], x2[,1], i, par, l)  
    Dkernx_theta_x <- sheatXrbfKernGradient(x1[,2], x2[,2], i, par, l)
    Dkern_Dtheta_x <- Dkern_Dtheta_x + kernt*Dkernx_theta_x
    Dkern_Dtheta_t <- Dkern_Dtheta_t + Dkernt$dkern_Dtheta_t*kernx
    Dkern_Dlambda <- Dkern_Dlambda + Dkernt$dkern_Dlambda*kernx
    Dkern_DD <- Dkern_DD + Dkernt$dkern_DD*kernx
  }
  kern <- sigma2*(2*S/l)*kern
  
  # gradients
  Dkern_DS <- kern/S
  Dkern_Dtheta_x <- sigma2*(2*S/l)*Dkern_Dtheta_x
  Dkern_Dtheta_t <- sigma2*(2*S/l)*Dkern_Dtheta_t
  Dkern_Dlambda <- sigma2*(2*S/l)*Dkern_Dlambda
  Dkern_DD <- sigma2*(2*S/l)*Dkern_DD
  attr(kern, "gradient") <- list(sigma2 = kern/sigma2,
                                 theta_t = Dkern_Dtheta_t,
                                 theta_x = Dkern_Dtheta_x,
                                 lambda = Dkern_Dlambda,
                                 D = Dkern_DD, S = Dkern_DS)  
  return(kern)
}

## Cross-covariance function Kuy
kernUY <- function(x1, x2, par) {
  # x1, x2 = matrices with  (time coordinate, space coordinate)
  # par = (sigma2, theta_t, theta_x, decay (lambda), diffusion (D), sensitivity (S))
  sigma2 <- par[1]; theta_t <- par[2]; theta_x <- par[3]
  lambda <- par[4]; D <- par[5]; S <- par[6]
  
  # fixing and precomputing some terms
  nterms <- 10
  l <- 1
  kern <- 0
  Dkern_Dtheta_x <- Dkern_Dtheta_t <- Dkern_Dlambda <- Dkern_DD <- 0
  
  # computing kernel and its gradient
  for(i in 1:nterms){
    # computing kernel
    kernt <- simXrbfKernCompute(x2[,1], x1[,1], i, par, l)
    kernx <- sheatXrbfKernCompute(x2[,2], x1[,2], i, par, l)
    kern <- kern + kernt*kernx
    
    # computing gradient
    Dkernt <- simXrbfKernGradient(x2[,1], x1[,1], i, par, l)  
    Dkernx_theta_x <- sheatXrbfKernGradient(x2[,2], x1[,2], i, par, l)
    Dkern_Dtheta_x <- Dkern_Dtheta_x + kernt*Dkernx_theta_x
    Dkern_Dtheta_t <- Dkern_Dtheta_t + Dkernt$dkern_Dtheta_t*kernx
    Dkern_Dlambda <- Dkern_Dlambda + Dkernt$dkern_Dlambda*kernx
    Dkern_DD <- Dkern_DD + Dkernt$dkern_DD*kernx
  }
  kern <- sigma2*(2*S/l)*t(kern)
  
  # gradients
  Dkern_DS <- kern/S
  Dkern_Dtheta_x <- sigma2*(2*S/l)*Dkern_Dtheta_x
  Dkern_Dtheta_t <- sigma2*(2*S/l)*Dkern_Dtheta_t
  Dkern_Dlambda <- sigma2*(2*S/l)*Dkern_Dlambda
  Dkern_DD <- sigma2*(2*S/l)*Dkern_DD
  attr(kern, "gradient") <- list(sigma2 = kern/sigma2,
                                 theta_t = Dkern_Dtheta_t,
                                 theta_x = Dkern_Dtheta_x,
                                 lambda = Dkern_Dlambda,
                                 D = Dkern_DD, S = Dkern_DS) 
  return(kern)
}

## Covariance function of the joint process
jointKern <- function(x1, x2, par) {
  # x1, x2 = matrices with  (time coordinate, space coordinate, label (0 for 'y', 1 for 'f'))
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
  
  # compute Kuu
  kernuu <- kernU(matrix(x1[indU1,1:2],ncol=2),matrix(x2[indU2,1:2],ncol=2),par[1:3])
  kern[indU1,indU2] <- kernuu
  Dkern_Dpar[indU1,indU2,1:3] <- array(unlist(attr(kernuu,'gradient')))
  
  # compute Kyy
  kernyy <- kernY(matrix(x1[indY1,1:2],ncol=2),matrix(x2[indY2,1:2],ncol=2),par)
  kern[indY1,indY2] <- kernyy
  Dkern_Dpar[indY1,indY2,] <- array(unlist(attr(kernyy,'gradient')))
  
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

#### auxiliar functions for the sim kernel ####
Hfun <- function(beta_s, t1, t2, theta_t){
  # t1, t2 = vectors with  (time coordinates)
  # precomputing some terms
  v <- 0.5*theta_t*beta_s
  diff_t <- outer(t2, t1,'-') 
  arg1 <- v + t1/theta_t
  arg2 <- v - diff_t/theta_t
  
  # computing the function H for the sim kernel
  H <- t(matrix(erf(arg1),length(t1),length(t2))) - erf(arg2) 
  return(Re(H))
}

HfunGradient <- function(beta_s, t1, t2, theta_t, omega_m2){
  # t1, t2 = vectors with  (time coordinates)
  # precomputing some terms
  v <- 0.5*theta_t*beta_s
  diff_t <- outer(t2, t1,'-') 
  arg1 <- v + t1/theta_t
  arg2 <- v - diff_t/theta_t
  
  # computing the gradient of the function H for the sim kernel
  dH_Dtheta_t <- t(matrix(exp(-arg1^2)*(0.5*beta_s - t1/theta_t^2),length(t1),length(t2))) -
    exp(-arg2^2)*(0.5*beta_s + diff_t/theta_t^2)
  dH_Dtheta_t <- (2/sqrt(pi))*dH_Dtheta_t
  
  dH_Dlambda <- t(matrix(exp(-arg1^2)*(0.5*theta_t),length(t1),length(t2))) -
    exp(-arg2^2)*(0.5*theta_t)
  dH_Dlambda <- (2/sqrt(pi))*dH_Dlambda
  
  dH_DD <- dH_Dlambda*omega_m2
  
  dH <- list(dH_Dtheta_t=dH_Dtheta_t,dH_Dlambda=dH_Dlambda,dH_DD=dH_DD)
  return(dH)
}

hfun  <- function(t1, t2, beta_s, beta_q, par){
  # t1, t2 = vectors with  (time coordinates)
  # par = (sigma2, theta_t, theta_x, decay (lambda), diffusion (D), sensitivity (S))  
  sigma2 <- par[1]; theta_t <- par[2]; theta_x <- par[3]
  lambda <- par[4]; D <- par[5]; S <- par[6]
  
  # precomputing some terms
  v <- 0.5*theta_t*beta_s
  diff_t <- outer(t2, t1,'-')
  t0 <- 0*t1
  
  # computing the function h for the sim kernel
  hpart1 <- exp(-beta_s*diff_t)*Hfun(beta_s,t1,t2,theta_t)
  hpart2 <- exp(-outer(beta_s*t2,beta_q*t1,'+'))*Hfun(beta_s,t0,t2,theta_t)
  h <- exp(v^2)*(hpart1 - hpart2)/(beta_s + beta_q)
  return(Re(h))
}

hfunGradient  <- function(t1, t2, beta_s, beta_q, omega_n2, omega_m2, par){
  # t1, t2 = vectors with  (time coordinates)
  # par = (sigma2, theta_t, theta_x, decay (lambda), diffusion (D), sensitivity (S))  
  sigma2 <- par[1]; theta_t <- par[2]; theta_x <- par[3]
  lambda <- par[4]; D <- par[5]; S <- par[6]
  
  # precomputing some terms
  v <- 0.5*theta_t*beta_s
  diff_t <- outer(t2, t1,'-')
  t0 <- 0*t1
  
  # computing the gradient of the function h for the sim kernel  
  hpart1 <- exp(-beta_s*diff_t)*Hfun(beta_s,t1,t2,theta_t)
  hpart2 <- exp(-outer(beta_s*t2,beta_q*t1,'+'))*Hfun(beta_s,t0,t2,theta_t)  
  h <- exp(v^2)*(hpart1 - hpart2)/(beta_s + beta_q)
  
  dH_t1t2 <- HfunGradient(beta_s,t1,t2,theta_t,omega_m2)
  dH_t0t2 <- HfunGradient(beta_s,t0,t2,theta_t,omega_m2)
  
  dhpart1_Dtheta_t <- exp(-beta_s*diff_t)*dH_t1t2$dH_Dtheta_t
  dhpart2_Dtheta_t <- exp(-outer(beta_s*t2,beta_q*t1,'+'))*dH_t0t2$dH_Dtheta_t
  dh_Dtheta_t <- 0.5*beta_s^2*theta_t*h + exp(v^2)*(dhpart1_Dtheta_t + dhpart2_Dtheta_t)/(beta_s + beta_q)
  
  dhpart1_Dlambda <- -diff_t*hpart1 + exp(-beta_s*diff_t)*dH_t1t2$dH_Dlambda
  dhpart2_Dlambda <- -outer(t2,t1,'+')*hpart2 + exp(-outer(beta_s*t2,beta_q*t1,'+'))*dH_t0t2$dH_Dlambda
  dh_Dlambda <- (0.5*beta_s*theta_t^2 - 2/(beta_s+beta_q))*h + exp(v^2)*(dhpart1_Dlambda - dhpart2_Dlambda)/(beta_s + beta_q)
  
  dhpart1_DD <- -diff_t*hpart1*omega_m2 + exp(-beta_s*diff_t)*dH_t1t2$dH_DD
  dhpart2_DD <- -outer(omega_m2*t2,omega_n2*t1,'+')*hpart2 + exp(-outer(beta_s*t2,beta_q*t1,'+'))*dH_t0t2$dH_DD
  dh_DD <- (0.5*beta_s*theta_t^2*omega_m2 - (omega_m2+omega_n2)/(beta_s+beta_q))*h + exp(v^2)*(dhpart1_DD - dhpart2_DD)/(beta_s + beta_q)
  
  dh <- list(dh_Dtheta_t=dh_Dtheta_t,dh_Dlambda=dh_Dlambda,dh_DD=dh_DD)
  return(dh)
}

## SIM kernel
simXsimKernCompute <- function(t1, t2, n, m, par, l){
  # t1, t2 = vectors with  (time coordinates)
  # par = (sigma2, theta_t, theta_x, decay (lambda), diffusion (D), sensitivity (S))  
  sigma2 <- par[1]; theta_t <- par[2]; theta_x <- par[3]
  lambda <- par[4]; D <- par[5]; S <- par[6]
  
  # precomputing some terms
  omega_n2 <- (n*pi/l)^2
  omega_m2 <- (m*pi/l)^2
  beta_q <- lambda + D*omega_n2
  beta_s <- lambda + D*omega_m2
  
  # computing the sim kernel
  kern <- t(hfun(t1, t2, beta_s, beta_q, par)) + hfun(t2, t1, beta_q, beta_s, par)
  kern <- 0.5*sqrt(pi)*theta_t*kern
  return(kern)
}

## gradients of the SIM kernel
simXsimKernGradient <- function(t1, t2, n, m, par, l){
  # t1, t2 = vectors with  (time coordinates)
  # par = (sigma2, theta_t, theta_x, decay (lambda), diffusion (D), sensitivity (S))  
  sigma2 <- par[1]; theta_t <- par[2]; theta_x <- par[3]
  lambda <- par[4]; D <- par[5]; S <- par[6]
  
  # precomputing some terms
  omega_n2 <- (n*pi/l)^2
  omega_m2 <- (m*pi/l)^2
  beta_q <- lambda + D*omega_n2
  beta_s <- lambda + D*omega_m2
  kern <- t(hfun(t1, t2, beta_s, beta_q, par)) + hfun(t2, t1, beta_q, beta_s, par)
  dh_t1t2 <- hfunGradient(t1, t2, beta_s, beta_q, omega_n2, omega_m2, par)
  dh_t2t1 <- hfunGradient(t2, t1, beta_q, beta_s, omega_n2, omega_m2, par)
  
  # computing the gradient of the sim kernel
  dkern_Dtheta_t <- t(dh_t1t2$dh_Dtheta_t) + dh_t2t1$dh_Dtheta_t
  dkern_Dtheta_t <- kern + theta_t*dkern_Dtheta_t
  dkern_Dtheta_t <- 0.5*sqrt(pi)*dkern_Dtheta_t
  dkern_Dlambda <- t(dh_t1t2$dh_Dlambda) + dh_t2t1$dh_Dlambda
  dkern_Dlambda <- 0.5*theta_t*sqrt(pi)*dkern_Dlambda
  dkern_DD <- t(dh_t1t2$dh_DD) + dh_t2t1$dh_DD
  dkern_DD <- 0.5*theta_t*sqrt(pi)*dkern_DD
  dkern <- list(dkern_Dtheta_t=dkern_Dtheta_t,dkern_Dlambda=dkern_Dlambda,dkern_DD=dkern_DD)
  return(dkern)
}

## Sheat kernel
sheatXsheatKernCompute <- function(x1, x2, n, m, par, l){
  # x1, x2 = vectors with  (space coordinates)
  # par = (sigma2, theta_t, theta_x, decay (lambda), diffusion (D), sensitivity (S))  
  sigma2 <- par[1]; theta_t <- par[2]; theta_x <- par[3]
  lambda <- par[4]; D <- par[5]; S <- par[6]
  
  # precomputing some terms
  j <- complex(real = 0, imaginary = 1)
  omega_n <- n*pi/l
  gamma_n <- j*omega_n
  z1_n <- 0.5*theta_x*gamma_n
  z2_n <- l/theta_x + z1_n
  wz1_n <- wofz(j*z1_n)
  wz2_n <- wofz(j*z2_n)
  Wox_n <- wz1_n - exp(-(l/theta_x)^2 - gamma_n*l)*wz2_n
  if(n==m){
    omega_m <- omega_n
    c <- Re(Wox_n) - Im(Wox_n)*( ((theta_x*n*pi)^2 + 2*l^2) / (2*pi*n*l^2)) 
    c <- 0.5*theta_x*sqrt(pi)*l*c + 0.5*theta_x^2*( exp(-(l/theta_x)^2)*cos(n*pi) - 1 )
  } else{
    omega_m <- m*pi/l
    gamma_m <- j*omega_m
    z1_m <- 0.5*theta_x*gamma_m
    z2_m <- l/theta_x + z1_m
    wz1_m <- wofz(j*z1_m)
    wz2_m <- wofz(j*z2_m)
    Wox_m <- wz1_m - exp(-(l/theta_x)^2 - gamma_m*l)*wz2_m
    c <- n*Im(Wox_m) - m*Im(Wox_n)
    c <- ( theta_x*l/(sqrt(pi)*(m^2-n^2)) )*c
  }
  # computing the sheat kernel
  kern  <- outer(sin(omega_n*x1),sin(omega_m*x2),'*')*c
  return(kern)
}

## gradients of the sheat kernel
sheatXsheatKernGradient <- function(x1, x2, n, m, par, l){
  # x1, x2 = vectors with  (space coordinates)
  # par = (sigma2, theta_t, theta_x, decay (lambda), diffusion (D), sensitivity (S))  
  sigma2 <- par[1]; theta_t <- par[2]; theta_x <- par[3]
  lambda <- par[4]; D <- par[5]; S <- par[6]
  
  # precomputing some terms
  j <- complex(real = 0, imaginary = 1)
  omega_n <- n*pi/l
  gamma_n <- j*omega_n
  z1_n <- 0.5*theta_x*gamma_n
  z2_n <- l/theta_x + z1_n
  wz1_n <- wofz(j*z1_n)
  wz2_n <- wofz(j*z2_n)
  Wox_n <- wz1_n - exp(-(l/theta_x)^2 - gamma_n*l)*wz2_n
  dWoxn_Dtheta_x <- 0.5*theta_x*gamma_n^2*wz1_n - gamma_n/sqrt(pi) -
    ( 2*l^2/theta_x^3 + 2*( (0.5*theta_x*gamma_n)^2 - (l/theta_x)^2  )/theta_x )*exp(-(l/theta_x)^2 - gamma_n*l)*wz2_n +
    2*(-l/theta_x^2 + 0.5*gamma_n )*exp(-(l/theta_x)^2 - gamma_n*l)/sqrt(pi)
  if(n==m){
    omega_m <- omega_n
    dC_Dtheta_x <- 0.5*sqrt(pi)*l*( Re(Wox_n)-Im(Wox_n)*(0.5*theta_x^2*n*pi/l^2 + 1/(n*pi)) ) +
      0.5*theta_x*sqrt(pi)*l*( Re(dWoxn_Dtheta_x) - Im(dWoxn_Dtheta_x)*(0.5*theta_x^2*n*pi/l^2 + 1/(n*pi)) - Im(Wox_n)*theta_x*n*pi/l^2 ) +
      theta_x*(exp(-(l/theta_x)^2)*cos(n*pi) - 1) + l^2*exp(-(l/theta_x)^2)*cos(n*pi)/theta_x 
  } else{
    omega_m <- m*pi/l
    gamma_m <- j*omega_m
    z1_m <- 0.5*theta_x*gamma_m
    z2_m <- l/theta_x + z1_m
    wz1_m <- wofz(j*z1_m)
    wz2_m <- wofz(j*z2_m)
    Wox_m <- wz1_m - exp(-(l/theta_x)^2 - gamma_m*l)*wz2_m
    dWoxm_Dtheta_x <- 0.5*theta_x*gamma_m^2*wz1_m - gamma_m/sqrt(pi) -
      ( 2*l^2/theta_x^3 + 2*( (0.5*theta_x*gamma_m)^2 - (l/theta_x)^2  )/theta_x )*exp(-(l/theta_x)^2 - gamma_m*l)*wz2_m +
      2*(-l/theta_x^2 + 0.5*gamma_m )*exp(-(l/theta_x)^2 - gamma_m*l)/sqrt(pi)
    dC_Dtheta_x <- n*Im(Wox_m) - m*Im(Wox_n) + theta_x*( n*Im(dWoxm_Dtheta_x)-m*Im(dWoxn_Dtheta_x) )
    dC_Dtheta_x <- ( l/(sqrt(pi)*(m^2-n^2)) )*dC_Dtheta_x
  }
  # computing the gradient of the sheat kernel
  Dkern_Dtheta_x <- outer(sin(omega_n*x1),sin(omega_m*x2),'*')*dC_Dtheta_x
  return(Dkern_Dtheta_x)
}

## cross-kernel between the sim and the SE kernel
simXrbfKernCompute <- function(t1, t2, n, par, l){
  # t1, t2 = vectors with  (time coordinates)
  # par = (sigma2, theta_t, theta_x, decay (lambda), diffusion (D), sensitivity (S))  
  sigma2 <- par[1]; theta_t <- par[2]; theta_x <- par[3]
  lambda <- par[4]; D <- par[5]; S <- par[6]
  
  # precomputing some terms
  omega_n2 <- (n*pi/l)^2
  beta_q <- lambda + D*omega_n2
  z1 <- 0.5*theta_t*beta_q
  diff_t <- outer(t1,t2,'-')
  
  # computing the simXrbf kernel  
  kern <- exp(z1^2-beta_q*diff_t)*Hfun(beta_q,t2,t1,theta_t)
  kern <- 0.5*sqrt(pi)*theta_t*kern
  return(kern)
}

## gradients of the cross-kernel between the sim and the SE kernel
simXrbfKernGradient <- function(t1, t2, n, par, l){
  # t1, t2 = vectors with  (time coordinates)
  # par = (sigma2, theta_t, theta_x, decay (lambda), diffusion (D), sensitivity (S))  
  sigma2 <- par[1]; theta_t <- par[2]; theta_x <- par[3]
  lambda <- par[4]; D <- par[5]; S <- par[6]
  
  # precomputing some terms
  omega_n2 <- (n*pi/l)^2
  beta_q <- lambda + D*omega_n2
  z1 <- 0.5*theta_t*beta_q
  diff_t <- outer(t1,t2,'-')
  kern <- exp(z1^2-beta_q*diff_t)*Hfun(beta_q,t2,t1,theta_t)
  kern <- 0.5*sqrt(pi)*theta_t*kern
  dH_t2t1 <- HfunGradient(beta_q,t2,t1,theta_t,n)
  
  # computing the gradient of the simXrbf kernel
  dkern_Dtheta_t <- kern/theta_t + (z1*beta_q)*kern + 0.5*sqrt(pi)*theta_t*exp(z1^2-beta_q*diff_t)*dH_t2t1$dH_Dtheta_t
  dkern_Dlambda <- (z1*theta_t-diff_t)*exp(z1^2-beta_q*diff_t)*Hfun(beta_q,t2,t1,theta_t) + exp(z1^2-beta_q*diff_t)*dH_t2t1$dH_Dlambda
  dkern_Dlambda <- 0.5*sqrt(pi)*theta_t*dkern_Dlambda
  dkern_DD <- dkern_Dlambda*omega_n2
  dkern <- list(dkern_Dtheta_t=dkern_Dtheta_t,dkern_Dlambda=dkern_Dlambda,dkern_DD=dkern_DD)
  return(dkern)
}

## cross-kernel between the sheat and the rbf kernel
sheatXrbfKernCompute <- function(x1, x2, n, par, l){
  # x1, x2 = vectors with  (space coordinates)
  # par = (sigma2, theta_t, theta_x, decay (lambda), diffusion (D), sensitivity (S))  
  sigma2 <- par[1]; theta_t <- par[2]; theta_x <- par[3]
  lambda <- par[4]; D <- par[5]; S <- par[6]
  
  # precomputing some terms  
  j <- complex(real = 0, imaginary = 1)
  omega_n <- n*pi/l
  gamma_n <- j*omega_n
  z1_n <- 0.5*theta_x*gamma_n + x2/theta_x
  z2_n <- z1_n - l/theta_x
  wz1_n <- wofz(j*z1_n)
  wz2_n <- wofz(j*z2_n)
  Wox_n <- exp(-((x2-l)/theta_x)^2 + gamma_n*l)*wz2_n - exp(-(x2/theta_x)^2)*wz1_n
  c <- 0.5*sqrt(pi)*theta_x*Im(Wox_n) 
  
  # computing the sheatXrbf kernel
  kern <- outer(sin(omega_n*x1),c)
  return(kern)
}

## gradients of the cross-kernel between the sheat and the rbf kernel
sheatXrbfKernGradient <- function(x1, x2, n, par, l){
  # x1, x2 = vectors with  (space coordinates)
  # par = (sigma2, theta_t, theta_x, decay (lambda), diffusion (D), sensitivity (S))  
  sigma2 <- par[1]; theta_t <- par[2]; theta_x <- par[3]
  lambda <- par[4]; D <- par[5]; S <- par[6]
  
  # precomputing some terms  
  j <- complex(real = 0, imaginary = 1)
  omega_n <- n*pi/l
  gamma_n <- j*omega_n
  z1_n <- 0.5*theta_x*gamma_n + x2/theta_x
  z2_n <- z1_n - l/theta_x
  wz1_n <- wofz(j*z1_n)
  wz2_n <- wofz(j*z2_n)
  Wox_n <- exp(-((x2-l)/theta_x)^2 + gamma_n*l)*wz2_n - exp(-(x2/theta_x)^2)*wz1_n
  arg1 <- x2/theta_x
  arg2 <- (x2-l)/theta_x
  dWoxn_Dtheta_x <- (2*arg2^2/theta_x + (2/theta_x)*((0.5*gamma_n*theta_x)^2-arg2^2))*exp(-arg2^2 + gamma_n*l)*wz2_n -
    (2/(theta_x*sqrt(pi)))*(0.5*gamma_n*theta_x - arg2)*exp(-arg2^2 + gamma_n*l) -
    (2*arg1^2/theta_x + (2/theta_x)*((0.5*gamma_n*theta_x)^2-arg1^2))*exp(-arg1^2)*wz1_n +
    (2/(theta_x*sqrt(pi)))*(0.5*gamma_n*theta_x - arg1)*exp(-arg1^2)
  dC_Dtheta_x <- Im(Wox_n) + theta_x*Im(dWoxn_Dtheta_x)
  dC_Dtheta_x <- 0.5*sqrt(pi)*dC_Dtheta_x
  
  # computing the gradient of the sheatXrbf kernel
  Dkern_Dtheta_x <-   outer(sin(omega_n*x1),dC_Dtheta_x)
  return(Dkern_Dtheta_x)
}

GPmRNAOpts <- function() {
  par <- c(sigma2 = 1, theta_t = 0.2, theta_x = 0.2,
           decay = 0.5, diffusion = 1e-4, sensitivity = 1)
  parLower <- c(sigma2 = 1e-10, theta_t = 0.01, theta_x = 0.01,
                decay = 1e-1, diffusion = 0, sensitivity = 0)
  parUpper <- c(sigma2 = Inf, theta_t = 2, theta_x = 2,
                decay = Inf, diffusion = 1e-3, sensitivity = Inf)
  parNames = c("sigma2", "theta_t", "theta_x", "decay", "diffusion", "sensitivity")
  return(list(par = par, parLower = parLower, parUpper = parUpper, parNames = parNames))
}
