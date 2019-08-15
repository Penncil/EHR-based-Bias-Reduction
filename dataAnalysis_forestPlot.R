## This r script is used to draw forest plot.
## Before drawing, two data frames are required to generate from the data.
## data frame one: point_est (columns are different coefficients, dim=(9,6) if there are five covariates.)
 #### row names 
 ## n=100 betahat
 ## n=100 gammabar
 ## n=100 beta_plugin
 ## n=100 betabar
 ## n=200 betahat
 ## n=200 gammabar
 ## n=200 beta_plugin
 ## n=200 betabar
 ## n=400 betahat
 ## n=400 gammabar
 ## n=400 beta_plugin
 ## n=400 betabar
## data frame two: CI (confidence interval) (dim = (9,12) if there are five covariates.)
  #### column names  beta0_lower_bound, beta0_upper_bound, beta1_lower_bound, beta1_upper_bound, etc.
  #### row names 
  ## n=100 betahat
  ## n=100 gammabar
  ## n=100 beta_plugin
  ## n=100 betabar
  ## n=200 betahat
  ## n=200 gammabar
  ## n=200 beta_plugin
  ## n=200 betabar
  ## n=400 betahat
  ## n=400 gammabar
  ## n=400 beta_plugin
  ## n=400 betabar


library(dplyr)
library(matlib)
library(dummies)



##---------------------------------------------------------------------------##
## Load and format BRAVA data
##---------------------------------------------------------------------------##

data <- read.csv("brava_penn_fake.csv")

# pre-process data
Y <- data[,"secBC"]    # true disease status
S <- data[,"SBCE_spec"]    # surrogate 
X <- data[,c("dxyear", "age", "Stage", "ER_PR")]  # predictors

# change the Stage 1, 2 to 0, 1
X$Stage <- X$Stage - 1  # change stage 2 to 1 and stage 1 to 0

# centerize the year
mean_year <- mean(X[,1])
X[,1] <- X[,1]-mean_year


# standardize the age
max <- max(X[,2])
X[,2] <- X[,2]/max

# change three-level categorical variables to two dummy variables 
# outcome 1: dummy1 = 1, dummy2 = 1 
# outcome 2: dummy1 = 0, dummy2 = 1
# outcome 3: dummy1 = 1, dummy2 = 0
x = X
x[,4:6] <- data.frame(dummy(x[,4]))
dummy1 <- rep(0,3152)
dummy2 <- rep(0,3152)

for (i in 1:3152){
  if((x[i,5] == 0)&(x[i,6]==1)){
    dummy2[i] = 1
  }else if ((x[i,5] == 1)&(x[i,6]==0)){
    dummy1[i] = 1
  }else{
    dummy1[i] = 0
    dummy2[i] = 0
  }
}
X_dummy <- X
X_dummy[,4] <- as.data.frame(dummy1)
X_dummy[,5] <- as.data.frame(dummy2) 
colnames(X_dummy)[4] <- "dummy1"
colnames(X_dummy)[5] <- "dummy2"

dxyear <- X_dummy[,1] # predictors in the entire cohort
age <- X_dummy[,2]
stage <- X_dummy[,3]
D1 <- X_dummy[,4]
D2 <- X_dummy[,5]

# true coefficients
glm0 <- glm(Y ~ dxyear+age+factor(stage)+factor(D1)+factor(D2), family = "binomial")
betaTrue <- glm0$coefficients
var_beta = summary(glm0)$coefficients[,2]**2
betareal.ci <- c(t(betaTrue + qnorm(0.975)*sqrt(var_beta)%*%t(c(-1,1))))

##---------------------------------------------------------------------------##
## Function to carry out bias reduction
##  Returns: naive estimate (gammabar), validation only estimate (betahat)
##  adjusted estimate (betabar)
##---------------------------------------------------------------------------##

proposed_var <- function(x_mul,m,val_X,val_Y,val_S,val_Z1,val_Z2,X,Z3,val_Z3,n,betaHat,gammaHat,gammaBar){

  nr <- dim(x_mul)[2]
  rho <- m/nrow(x_mul)
    
  info_matrix1 <- matrix(0 , nrow = nr, ncol = nr)     # information matrix

  for (i in 1:nr){
    for (j in 1:nr){
      if (i == j){
        info_matrix1[i,j] <- sum( val_X[,i]^2 * exp(val_Z1) / (1+exp(val_Z1))^2 ) /m
      }else
        info_matrix1[i,j] <- sum( val_X[,i] * val_X[,j] * exp(val_Z1) / (1+exp(val_Z1))^2) /m
      info_matrix1[j,i] <- info_matrix1[i,j]
    }
  }

  info_matrix2 <- matrix(0 , nrow = nr, ncol = nr)     # information matrix
  for (i in 1:nr){
    for (j in 1:nr){
      if (i == j){
        info_matrix2[i,j] <- sum(X[,i]^2 * exp(Z3) / (1+exp(Z3))^2 ) /n
      }else
        info_matrix2[i,j] <- sum(X[,i] * X[,j] * exp(Z3) / (1+exp(Z3))^2) /n
      info_matrix2[j,i] <- info_matrix2[i,j]
    }
  }

  p1 <- matrix(0, nrow = nr, ncol =nr)
  p2 <- matrix(0, nrow = nr, ncol =nr)
  p3 <- matrix(0, nrow = nr, ncol =nr)

  for (i in 1:m){

    phi1 <- (val_Y[i] - (exp(val_Z1[i])/(1+exp(val_Z1[i])))) %*% as.numeric(val_X[i,])
    phi1 <- inv(info_matrix1) %*% as.vector(phi1)
    p1 <- p1 + phi1 %*% t(phi1)

    phi2 <- (val_S[i] - (exp(val_Z3[i])/(1+exp(val_Z3[i])))) %*% as.numeric(val_X[i,])
    phi2 <- inv(info_matrix2) %*% as.vector(phi2)
    p2 <- p2 + phi2 %*% t(phi2)

    p3 <- p3 + phi1 %*% t(phi2)

  }

  sigma <- p1 / m
  sigma_star <- ((1-rho)) * p2 / m
  omega <- (1-rho) * p3 / m

  betaBar = t(betaHat - omega %*% inv(sigma_star) %*% (gammaHat - gammaBar))
  var_betaBar = (sigma - omega %*% inv(sigma_star) %*% t(omega)) /m
  var_betaBar = diag(var_betaBar)
  
  return(list(betaBar = betaBar, var_betaBar = var_betaBar))
}


plugin <- function(val_X, val_S, val_Y, x_mul, S, betaHat){

  d = dim(x_mul)[2]
  nr <- d

  # calculate plugin estimator and variance
  ##################################################################
  ##### Plugin estimator
  t <- table(val_S,val_Y) # only validation set of y is obvered, so we can only use y in the validation set to do the estimation
  est_pooled_a1 <- t[2,2]/sum(t[,2])   # estimate sensitivity
  est_pooled_a2 <- t[1,1]/sum(t[,1])   # estimate specificity

  ## function for plugin estimator
  nloglik=function(t)
  {
    beta = t[1:d];
    alpha1 = t[length(t)-1]; alpha2 = t[length(t)]
    tem= beta %*% t(x_mul)
    p = (1-alpha2)+(alpha1+alpha2-1)*exp(tem)/(1+exp(tem))
    likelihood= S*log(p)+(1-S)*log(1-p)
    return(-sum(likelihood))
  }

  beta_plugin_temp <- nlminb(c(betaHat,est_pooled_a1,est_pooled_a2),
                             nloglik,lower=c(rep(-5,d),est_pooled_a1,est_pooled_a2),
                             upper=c(rep(5,d),est_pooled_a1,est_pooled_a2))

  beta_plugin <-beta_plugin_temp$par[1:d]

  #### variance
  Zp <- beta_plugin %*% t(x_mul)

  info_matrix <- matrix(0 , nrow = nr, ncol = nr)     # information matrix
  for (i in 1:nr){
    for (j in 1:nr){
      if (i == j){
        info_matrix[i,j] <- sum( x_mul[,i]^2 * exp(Zp) / (1+exp(Zp))^2 ) /n
      }else
        info_matrix[i,j] <- sum( x_mul[,i] * x_mul[,j] * exp(Zp) / (1+exp(Zp))^2) /n
      info_matrix[j,i] <- info_matrix[i,j]
    }
  }

  var_beta_plugin_matrix <- inv(info_matrix)/n
  var_beta_plugin <- diag(var_beta_plugin_matrix)

  return(list(beta_plugin = beta_plugin, var_beta_plugin=var_beta_plugin))
}


par_est <- function(m){
  n = dim(data)[1]

  sample <- F
  while (!sample){
    id <- sample(n, m) 
    
    val_Y <- Y[id]  # validation set of Y
    val_S <- S[id]  # validation set of S
    val_X <- X_dummy[id,] # validation set of X
    
    v_dxyear <- dxyear[id] # predictors in the validation set
    v_age <- age[id]
    v_stage <- stage[id]
    v_D1 <- D1[id]
    v_D2 <- D2[id]
  
    
    #check for empty cells
    if (min(c(c(table(val_Y,v_D1)),c(table(val_Y,v_D2)),c(table(val_S,v_D2)),c(table(val_S,v_D1))))>0) sample <- T 
  }
  
  # feed the glm
  ###############################################################

  # with only validation Y and X
  glm1 <- glm(val_Y~v_dxyear+v_age+factor(v_stage)+factor(v_D1)+factor(v_D2), family = "binomial")
  betaHat <- glm1$coefficients
  var_beta_hat = summary(glm1)$coefficients[,2]**2
  
  # validation S and X
  glm2 <- glm(val_S~v_dxyear+v_age+factor(v_stage)+factor(v_D1)+factor(v_D2), family = "binomial")
  gammaHat <- glm2$coefficients
  var_gamma_hat = summary(glm2)$coefficients[,2]**2
  
  # all S and X
  glm3 <- glm(S ~ dxyear+age+factor(stage)+factor(D1)+factor(D2), family = "binomial")
  gammaBar <- glm3$coefficients
  var_gamma_bar = summary(glm3)$coefficients[,2]**2
  
  
  # calcualte Z and get validation set of Z
  ##################################################################
  x_mul <- cbind(1,X_dummy)
  Z1 <- betaHat %*% t(x_mul)
  Z2 <- gammaHat %*% t(x_mul)
  Z3 <- gammaBar %*% t(x_mul)
  val_Z1 <- Z1[id]
  val_Z2 <- Z2[id]
  val_Z3 <- Z3[id]
  val_X <- x_mul[id,]

  
  # calculate proposed variance
  ##################################################################
  adj <- proposed_var(x_mul,m,val_X,val_Y,val_S,val_Z1,val_Z2,x_mul,Z3,val_Z3,n,betaHat,gammaHat,gammaBar)
  
  var_betaBar <- adj$var_betaBar
  betaBar <- adj$betaBar
  
  plugin <- plugin(val_X, val_S, val_Y, x_mul, S, betaHat)
  
  beta_plugin <- plugin$beta_plugin
  var_beta_plugin <- plugin$var_beta_plugin
  
  est <- rbind(betaHat, gammaBar, beta_plugin, betaBar)
  
  betahat.ci <- c(t(betaHat + qnorm(0.975)*sqrt(var_beta_hat)%*%t(c(-1,1))))
  gammabar.ci <- c(t(gammaBar + qnorm(0.975)*sqrt(var_gamma_bar)%*%t(c(-1,1))))
  beta_plugin.ci <- c(t(beta_plugin + qnorm(0.975)*sqrt(var_beta_plugin)%*%t(c(-1,1))))
  betabar.ci <- c(t(c(betaBar) + qnorm(0.975)*sqrt(var_betaBar)%*%t(c(-1,1))))
  
  ci <- rbind(betahat.ci, gammabar.ci, beta_plugin.ci, betabar.ci)
  var <- rbind(var_beta_hat, var_gamma_bar, var_beta_plugin, var_betaBar)
  
  return(list(est = est, CI = ci, VAR=var))
}


##---------------------------------------------------------------------------##
## Apply bias reduction method to get point estimates and standard errors
##  for validation sample sizes: 200, 400, 800
##---------------------------------------------------------------------------##

set.seed(4321)
est1 <- par_est(200)
set.seed(4321)
est2 <- par_est(400)
set.seed(4321)
est4 <- par_est(800)

point_est <- rbind(est1$est, est2$est, est4$est)
CI <- rbind(est1$CI, est2$CI, est4$CI)


### forest plot
forest.plot <- function(i){
  ## xlim depends on the value of CI and point_est
  plot(0,0,type="n", xlab="", ylab="", yaxt="n", xaxt="n", xaxs="i", yaxs="i",ylim=c(1,7),xlim=c(min(CI[,(2*i-1)])-0.5,max(CI[,(2*i)])+0.5), main="")
  axis(1)
  mtext("95% Confidence Interval", outer=FALSE, line=0.4,cex=0.9)
  
  
  # segments for betahat
  segments(CI[c(1,5,9),2*i-1], c(2.2,4.2,6.2),CI[c(1,5,9),2*i], c(2.2,4.2,6.2), lwd=2)
  segments(point_est[c(1,5,9),i], c(2.2,4.2,6.2)-0.1,point_est[c(1,5,9),i], c(2.2,4.2,6.2)+0.1, lwd=2)
  segments(CI[c(1,5,9),2*i-1], c(2.2,4.2,6.2)-0.05,CI[c(1,5,9),2*i-1], c(2.2,4.2,6.2)+0.05, lwd=2)
  segments(CI[c(1,5,9),2*i], c(2.2,4.2,6.2)-0.05,CI[c(1,5,9),2*i], c(2.2,4.2,6.2)+0.05, lwd=2)
  
  # segments for gammabar
  segments(CI[c(2,6,10),2*i-1], c(1.9,3.9,5.9),CI[c(2,6,10),2*i], c(1.9,3.9,5.9), lwd=1, lty = 1)
  segments(point_est[c(2,6,10),i], c(1.9,3.9,5.9)-0.1,point_est[c(2,6,10),i], c(1.9,3.9,5.9)+0.1, lwd=1)
  segments(CI[c(2,6,10),2*i-1], c(1.9,3.9,5.9)-0.05,CI[c(2,6,10),2*i-1], c(1.9,3.9,5.9)+0.05, lwd=1)
  segments(CI[c(2,6,10),2*i], c(1.9,3.9,5.9)-0.05,CI[c(2,6,10),2*i], c(1.9,3.9,5.9)+0.05, lwd=1)
  
  # segments for beta_plugin
  segments(CI[c(3,7,11),2*i-1], c(1.6,3.6,5.6),CI[c(3,7,11),2*i], c(1.6,3.6,5.6), lwd=1, lty = 2)
  segments(point_est[c(3,7,11),i], c(1.6,3.6,5.6)-0.1,point_est[c(3,7,11),i], c(1.6,3.6,5.6)+0.1, lwd=1)
  segments(CI[c(3,7,11),2*i-1], c(1.6,3.6,5.6)-0.05,CI[c(3,7,11),2*i-1], c(1.6,3.6,5.6)+0.05, lwd=1)
  segments(CI[c(3,7,11),2*i], c(1.6,3.6,5.6)-0.05,CI[c(3,7,11),2*i], c(1.6,3.6,5.6)+0.05, lwd=1)
  
  # segments for betabar
  segments(CI[c(4,8,12),2*i-1], c(1.3,3.3,5.3),CI[c(4,8,12),2*i], c(1.3,3.3,5.3), lwd=1,col='red', lty = 3)
  segments(point_est[c(4,8,12),i], c(1.3,3.3,5.3)-0.1,point_est[c(4,8,12),i], c(1.3,3.3,5.3)+0.1, lwd=1,col="red")
  segments(CI[c(4,8,12),2*i-1], c(1.3,3.3,5.3)-0.05,CI[c(4,8,12),2*i-1], c(1.3,3.3,5.3)+0.05, lwd=1,col="red")
  segments(CI[c(4,8,12),2*i], c(1.3,3.3,5.3)-0.05,CI[c(4,8,12),2*i], c(1.3,3.3,5.3)+0.05, lwd=1,col="red")
  
  # real value of the ith data (confidence interval)
  abline(v=betaTrue[i],lty=2,col="blue")
  abline(v=betareal.ci[2*i-1],lty=2,col="green")
  abline(v=betareal.ci[2*i],lty=2,col="green")

  
  if (i != 2){
    u = par("usr")
    segments(u[1], c(2,4,6), u[1]-0.18, c(2,4,6), xpd=TRUE)
    text(u[1]-0.18,2,'n=200',xpd=T,adj=1,cex=1)
    text(u[1]-0.18,4,'n=400',xpd=T,adj=1,cex=1)
    text(u[1]-0.18,6,'n=800',xpd=T,adj=1,cex=1)
    
    legend("topleft",legend=c( expression(paste(beta[paste("V")])),
                               expression(paste(gamma[paste("F")])),
                               expression(paste(beta[paste("P")])),
                               expression(paste(beta[paste("A")]))),
           lty=c(1,1,2,3),lwd=c(2,1,1,1),cex=0.8,col = c(1,1,1,2),bty="n")
  }else{
    u = par("usr")
    segments(u[1], c(2,4,6), u[1]-0.03, c(2,4,6), xpd=TRUE)
    text(u[1]-0.03,2,'n=200',xpd=T,adj=1,cex=1)
    text(u[1]-0.03,4,'n=400',xpd=T,adj=1,cex=1)
    text(u[1]-0.03,6,'n=800',xpd=T,adj=1,cex=1)
    
    legend("topleft",legend=c( expression(paste(beta[paste("V")])),
                               expression(paste(gamma[paste("F")])),
                               expression(paste(beta[paste("P")])),
                               expression(paste(beta[paste("A")]))),
           lty=c(1,1,2,3),lwd=c(2,1,1,1),cex=0.8,col = c(1,1,1,2),bty="n")
  }
  

}

par(mfrow=c(3,2))
forest.plot(1) 
forest.plot(2)
forest.plot(3)
forest.plot(4)
forest.plot(5)
forest.plot(6)
# dev.off()


