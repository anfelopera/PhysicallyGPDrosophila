library(kergp)
library(MASS)
library(plot3D)
library(pracma)
library(viridis)
library(tikzDevice)
library(DiceDesign)

rm(list=ls())
set.seed(2)

viridisPath <- viridis(30)
source("GPDrosMutils.R")
source("GPGene.R")

#############################################################
## Defining the covariance functions for kergp
kernOpts <- GPGeneOpts() # default options

covY <- covMan(kernel = kernY, inputs = c('t','x'),
               hasGrad = TRUE, acceptMatrix = TRUE,
               label = "Kyy",  d = 2,
               par = kernOpts$par[1:3],
               parLower = kernOpts$parLower[1:3],
               parUpper = kernOpts$parUpper[1:3],
               parNames = kernOpts$parNames[1:3]
)

covU <- covMan(kernel = kernU, inputs = c('t','x'),
               hasGrad = TRUE, acceptMatrix = TRUE,
               label = "Kuu", d = 2,
               par = kernOpts$par,
               parLower = kernOpts$parLower,
               parUpper = kernOpts$parUpper,
               parNames = kernOpts$parNames
)

covYU <- covMan(kernel = kernYU, inputs =c('t','x'),
                hasGrad = TRUE, acceptMatrix = TRUE,
                label = "Kyu", d = 2,
                par = kernOpts$par,
                parLower = kernOpts$parLower,
                parUpper = kernOpts$parUpper,
                parNames = kernOpts$parNames
)

covUY<- covMan(kernel = kernUY, inputs =c('t','x'),
               hasGrad = TRUE, acceptMatrix = TRUE,
               label = "Kuy", d = 2,
               par = kernOpts$par,
               parLower = kernOpts$parLower,
               parUpper = kernOpts$parUpper,
               parNames = kernOpts$parNames
)

jointCov <- covMan(kernel = jointKern, inputs =c('t','x','l'),
                   hasGrad = TRUE, acceptMatrix = TRUE,
                   label = "Joint kernel", d = 3,
                   par = kernOpts$par,
                   parLower = kernOpts$parLower,
                   parUpper = kernOpts$parUpper,
                   parNames = kernOpts$parNames
)

#############################################################
## loading data from data-set
nb_exp <- 1
toyName <- paste("demToy", nb_exp, "GPGene", sep = "")
dataName <- paste(toyName,".csv", sep = "")

if (file.exists(dataName)){
  data <- read.csv(dataName)
} else{
  ## Generating synthetic dataset
  l <- 1
  tgrid <- seq(0, 1, 0.025)
  xgrid <- seq(0, 1, 0.025)
  X <- expand.grid(t = tgrid, x = xgrid, l = c(0,1))
  nsamples <- 1

  jointCov@par <- c(sigma2 = 1, theta_t = 0.3, theta_x = 0.3,
                    decay = 0.1, diffusion = 0.01, sensitivity = 1)
  K <- covMat(object = jointCov, X = X, Xnew = X)
  set.seed(7); YUsamples <- mvrnorm(nsamples, rep(0, nrow(K)), K)
  Y <- YUsamples[1:(length(tgrid)*length(xgrid))]
  U <- YUsamples[length(tgrid)*length(xgrid) + 1:(length(tgrid)*length(xgrid))]
  data <- cbind(expand.grid(t=tgrid,x=xgrid), U, Y)
  write.csv(data,dataName, row.names = FALSE)
}

#############################################################
## plotting the results
tgrid <- unique(data[,1])
xgrid <- unique(data[,2])

# filename <- paste(toyName, "output.tex", sep = "")
# tikz(filename, standAlone = TRUE, width=7.0, height=7.0)
# par(cex.axis = 1.5, cex.lab = 2.0, cex.main = 2.5, lwd = 5)
image2D(matrix(data[,4], nrow = length(tgrid)), tgrid, xgrid,
        resfac = 3, col = viridisPath, xlab = "$t$", ylab = "$x$",
        main = "$y(t,x)$", xlim = range(tgrid), ylim = range(xgrid))
# dev.off()
# tools::texi2dvi(filename,pdf=T,clean=TRUE)

# filename <- paste(toyName, "force.tex", sep = "")
# tikz(filename, standAlone = TRUE, width=7.0, height=7.0)
# par(cex.axis = 1.5, cex.lab = 2.0, cex.main = 2.5, lwd = 5)
image2D(matrix(data[,3], nrow = length(tgrid)), tgrid, xgrid,
        resfac = 3, col = viridisPath, xlab = "$t$", ylab = "$x$",
        main = "$u(t,x)$", xlim = range(tgrid), ylim = range(xgrid))
# dev.off()
# tools::texi2dvi(filename,pdf=T,clean=TRUE)

#############################################################
## inference procedure
# training dataset 
expNum <- 4
l <- 1
tgrid <- unique(data[,1])
xgrid <- unique(data[,2])

idxtemp <- lhsDesign(10*expNum, 2, seed = 7)
idxtemp <- maximinSA_LHS(idxtemp$design)
idxtemp <- idxtemp$design
idxtemp[,1] <- round((length(tgrid)-1)*idxtemp[,1]+1)
idxtemp[,2] <- round((length(xgrid)-1)*idxtemp[,2]+1)
idx <- idxtemp[,1] + length(tgrid)*(idxtemp[,2]-1)
dataTrain <- data[sort(idx),]

param <- c(sigma2 = 1, theta_t = 0.3, theta_x = 0.3,
           decay = 0.1, diffusion = 0.01, sensitivity = 1)
xtrain <- dataTrain[,1:2]
utrain <- data.frame(u=dataTrain[,3])
ytrain <- data.frame(y=dataTrain[,4])
xtest <- data[,1:2]

# joint model using data from the output 
T_temp <- expand.grid(t=dataTrain[,1],l=0)
X_temp <- expand.grid(x=dataTrain[,2],l=0)
Xtrain <- data.frame(t=T_temp[,1],x=X_temp[,1],l=T_temp[,2])
Jtrain <- data.frame(j=c(dataTrain[,4]))
jointCov@par <- param
gprJ <- gp(formula = j ~ 1, data = data.frame(Jtrain, Xtrain),
           inputs = c('t','x','l'), cov = jointCov,
           beta = 0, estim = FALSE)
summary(gprJ)

T_temp <- expand.grid(t=data[,1],l=c(0,1))
X_temp <- expand.grid(x=data[,2],l=c(0,1))
Xtest <- data.frame(t=T_temp[,1],x=X_temp[,1],l=T_temp[,2])
predJ <- predict(gprJ,Xtest,type='SK')
Ypred <- predJ$mean[1:(length(tgrid)*length(xgrid))]
Upred <- predJ$mean[length(tgrid)*length(xgrid) + 1:(length(tgrid)*length(xgrid))]
Yvar <- predJ$sd[1:(length(tgrid)*length(xgrid))]^2
Uvar<- predJ$sd[length(tgrid)*length(xgrid) + 1:(length(tgrid)*length(xgrid))]^2

par(mfrow=c(1,1))
# filename <- paste(toyName, "forceDataO", expNum, ".tex", sep = "")
# tikz(filename, standAlone = TRUE, width=7.0, height=7.0)
# par(cex.axis = 1.5, cex.lab = 2.0, cex.main = 2.5, lwd = 5)
image2D(matrix(Upred, nrow = length(tgrid)), tgrid, xgrid,
        resfac = 3, col = viridisPath, xlab = "$t$", ylab = "$x$",
        main = "$u(t,x)$", xlim = range(tgrid), ylim = range(xgrid))
# dev.off()
# tools::texi2dvi(filename,pdf=T,clean=TRUE)
round(errorMeasureRegress(data[idx,3], data[-idx,3], Upred[-idx], Uvar[-idx]), digits = 3)[c(5,7)]

# filename <- paste(toyName, "outputDataO", expNum, ".tex", sep = "")
# tikz(filename, standAlone = TRUE, width=7.0, height=7.0)
# par(cex.axis = 1.5, cex.lab = 2.0, cex.main = 2.5, lwd = 5)
image2D(matrix(Ypred, nrow = length(tgrid)), tgrid, xgrid,
        resfac = 3, col = viridisPath, xlab = "$t$", ylab = "$x$",
        main = "$y(t,x)$", xlim = range(tgrid), ylim = range(xgrid))
points(Xtrain, pch = 20, cex = 2, lwd = 2, col = "white")
points(Xtrain, cex = 1.5, lwd = 2)
# dev.off()
# tools::texi2dvi(filename,pdf=T,clean=TRUE)
round(errorMeasureRegress(data[idx,4], data[-idx,4], Ypred[-idx], Yvar[-idx]), digits = 3)[c(5,7)]

# joint model using data from the force 
T_temp <- expand.grid(t=dataTrain[,1],l=1)
X_temp <- expand.grid(x=dataTrain[,2],l=1)
Xtrain <- data.frame(t=T_temp[,1],x=X_temp[,1],l=T_temp[,2])
Jtrain <- data.frame(j=dataTrain[,3])
jointCov@par <- param
gprJ <- gp(formula = j ~ 1, data = data.frame(Jtrain, Xtrain),
           inputs = c('t','x','l'), cov = jointCov,
           beta = 0, estim = FALSE)
summary(gprJ)

T_temp <- expand.grid(t=data[,1],l=c(0,1))
X_temp <- expand.grid(x=data[,2],l=c(0,1))
Xtest <- data.frame(t=T_temp[,1],x=X_temp[,1],l=T_temp[,2])
predJ <- predict(gprJ,Xtest,type='SK')
Ypred <- predJ$mean[1:(length(tgrid)*length(xgrid))]
Upred <- predJ$mean[length(tgrid)*length(xgrid) + 1:(length(tgrid)*length(xgrid))]
Yvar <- predJ$sd[1:(length(tgrid)*length(xgrid))]^2
Uvar<- predJ$sd[length(tgrid)*length(xgrid) + 1:(length(tgrid)*length(xgrid))]^2

par(mfrow=c(1,1))
# filename <- paste(toyName, "forceDataF", expNum, ".tex", sep = "")
# tikz(filename, standAlone = TRUE, width=7.0, height=7.0)
# par(cex.axis = 1.5, cex.lab = 2.0, cex.main = 2.5, lwd = 5)
image2D(matrix(Upred, nrow = length(tgrid)), tgrid, xgrid,
        resfac = 3, col = viridisPath, xlab = "$t$", ylab = "$x$",
        main = "$u(t,x)$", xlim = range(tgrid), ylim = range(xgrid))
points(Xtrain, pch = 20, cex = 2, lwd = 2, col = "white")
points(Xtrain, cex = 1.5, lwd = 2)
# dev.off()
# tools::texi2dvi(filename,pdf=T,clean=TRUE)
round(errorMeasureRegress(data[idx,3], data[-idx,3], Upred[-idx], Uvar[-idx]), digits = 3)[c(5,7)]

# filename <- paste(toyName, "outputDataF", expNum, ".tex", sep = "")
# tikz(filename, standAlone = TRUE, width=7.0, height=7.0)
# par(cex.axis = 1.5, cex.lab = 2.0, cex.main = 2.5, lwd = 5)
image2D(matrix(Ypred, nrow = length(tgrid)), tgrid, xgrid,
        resfac = 3, col = viridisPath, xlab = "$t$", ylab = "$x$",
        main = "$y(t,x)$", xlim = range(tgrid), ylim = range(xgrid))
# dev.off()
# tools::texi2dvi(filename,pdf=T,clean=TRUE)
round(errorMeasureRegress(data[idx,4], data[-idx,4], Ypred[-idx], Yvar[-idx]), digits = 3)[c(5,7)]

# joint model using data from the output and input 
T_temp <- expand.grid(t=dataTrain[,1],l=c(0,1))
X_temp <- expand.grid(x=dataTrain[,2],l=c(0,1))
Xtrain <- data.frame(t=T_temp[,1],x=X_temp[,1],l=T_temp[,2])
Jtrain <- data.frame(j=c(dataTrain[,4], dataTrain[,3]))
jointCov@par <- param
gprJ <- gp(formula = j ~ 1, data = data.frame(Jtrain, Xtrain),
           inputs = c('t','x','l'), cov = jointCov,
           beta = 0, estim = FALSE)
summary(gprJ)

T_temp <- expand.grid(t=data[,1],l=c(0,1))
X_temp <- expand.grid(x=data[,2],l=c(0,1))
Xtest <- data.frame(t=T_temp[,1],x=X_temp[,1],l=T_temp[,2])
predJ <- predict(gprJ,Xtest,type='SK')
Ypred <- predJ$mean[1:(length(tgrid)*length(xgrid))]
Upred <- predJ$mean[length(tgrid)*length(xgrid) + 1:(length(tgrid)*length(xgrid))]
Yvar <- predJ$sd[1:(length(tgrid)*length(xgrid))]^2
Uvar<- predJ$sd[length(tgrid)*length(xgrid) + 1:(length(tgrid)*length(xgrid))]^2

par(mfrow=c(1,1))
# filename <- paste(toyName, "forceDataB", expNum, ".tex", sep = "")
# tikz(filename, standAlone = TRUE, width=7.0, height=7.0)
# par(cex.axis = 1.5, cex.lab = 2.0, cex.main = 2.5, lwd = 5)
image2D(matrix(Upred, nrow = length(tgrid)), tgrid, xgrid, resfac = 3,
        col = viridisPath, xlab = "$t$", ylab = "$x$",
        main = "$u(t,x)$", xlim = range(tgrid), ylim = range(xgrid))
points(Xtrain, pch = 20, cex = 2, lwd = 2, col = "white")
points(Xtrain, cex = 1.5, lwd = 2)
# dev.off()
# tools::texi2dvi(filename,pdf=T,clean=TRUE)
round(errorMeasureRegress(data[idx,3], data[-idx,3], Upred[-idx], Uvar[-idx]), digits = 3)[c(5,7)]

# filename <- paste(toyName, "outputDataB", expNum, ".tex", sep = "")
# tikz(filename, standAlone = TRUE, width=7.0, height=7.0)
# par(cex.axis = 1.5, cex.lab = 2.0, cex.main = 2.5, lwd = 5)
image2D(matrix(Ypred, nrow = length(tgrid)), tgrid, xgrid,
        resfac = 3, col = viridisPath, xlab = "$t$", ylab = "$x$",
        main = "$y(t,x)$", xlim = range(tgrid), ylim = range(xgrid))
points(Xtrain, pch = 20, cex = 2, lwd = 2, col = "white")
points(Xtrain, cex = 1.5, lwd = 2)
# dev.off()
# tools::texi2dvi(filename,pdf=T,clean=TRUE)
round(errorMeasureRegress(data[idx,4], data[-idx,4], Ypred[-idx], Yvar[-idx]), digits = 3)[c(5,7)]
