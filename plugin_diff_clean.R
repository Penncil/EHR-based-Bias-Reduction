library(matlib)
library(ggplot2)

simulation_plugin_diff <- function(m,n,niters,beta0,beta1){
  
  alpha1 = 0.9          # sensitivity pr(s=1|y=1) = 0.9 for non-exposure group
  alpha2 = 0.95           # specificity pr(s=0|y=0) = 0.95  for non-exposure group
  
  myest=array(dim=c(niters,6))
  b1 =array(dim=c(niters,2))    # array of beta_hat
  b2 =array(dim=c(niters,2))    # array of beta_bar
  r1 = array(dim=c(niters,2))   # array of gamma_bar
  b_real =array(dim=c(niters,2)) # array of real beta
  b3 =array(dim=c(niters,1))
  
  set.seed(1)
  
  for (iter in 1:niters){
    

    
    # ## differential setting 1
    # x = rnorm(n)   # some continuous variables
    x = rbinom(n,1,0.4)
    z = beta0 + beta1*x
    pr = 1/(1+exp(-z))     # pass through an inv-logit function
    y = rbinom(n,1,pr)     # bernoulli response variable


    # # sensiticity 0.9 -> 0.95
    # # specificity 0.95 -> 0.9
    # a1 <- 0.95
    # a2 <- 0.9
    # 
    # # sensiticity 0.9 -> 0.85
    # # specificity 0.95 -> 0.99
    a1 <- 0.85
    a2 <- 0.99
    # 
    pr_s = vector(mode = "numeric",length = n)
    pr_s[x==1] = a1*(y[x==1]==1) + (1-a2)*(y[x==1]==0)
    pr_s[x==0] = alpha1*(y[x==0]==1) + (1-alpha2)*(y[x==0]==0)
    s = rbinom(n,1,pr_s)
    
    m = m                 # size of validation set
    k = m/n
    
    v <- sample(1:n, m, replace = FALSE)
    
    t <- table(s[v],y[v]) # only validation set of y is obvered, so we can only use y in the validation set to do the estimation
    est_pooled_a1 <- t[2,2]/sum(t[,2])   # estimate sensitivity
    est_pooled_a2 <- t[1,1]/sum(t[,1])   # estimate specificity

    
    nloglik=function(t)
    {
      beta0 = t[1]; beta1 = t[2];
      alpha1 = t[3]; alpha2 = t[4]
      tem=beta0+beta1*x
      p = (1-alpha2)+(alpha1+alpha2-1)*exp(tem)/(1+exp(tem))
      likelihood= s*log(p)+(1-s)*log(1-p)
      return(-sum(likelihood))
    }
    

    val_x = x[v]           # validation set of x
    val_y = y[v]           # validation set of y
    val_s = s[v]           # validation set of s
    df_yval <- data.frame(val_y,val_x)      # dataframe of validation x and y
    df_sval <- data.frame(val_s,val_x)      # dataframe of validation x and s
    df_s <- data.frame(s, x)           # dataframe of full x and s
    df <- data.frame(y,x)
    
    
    # feed it to glm and get mle of estimators:
    glm0 <- glm(y~x, data = df, family = "binomial")
    beta_real <- glm0$coefficients
    
    glm1 <- glm(val_y~val_x, data=df_yval, family="binomial")
    beta_hat <- glm1$coefficients
    
    glm2 <- glm(val_s~val_x,data=df_sval,family="binomial")
    gamma_hat <- glm2$coefficients
    
    glm3 <- glm(s~x,data=df_s,family="binomial")
    gamma_bar <- glm3$coefficients
   
    beta_plugin_temp <- nlminb(c(beta0,beta1,est_pooled_a1,est_pooled_a2),
                               nloglik,lower=c(-5,-5,est_pooled_a1,est_pooled_a2),
                               upper=c(5,5,est_pooled_a1,est_pooled_a2))
    
    beta_plugin <-beta_plugin_temp$par[2]

    # variance of beta_plugin
    z_p <- beta_plugin_temp$par[1]+beta_plugin_temp$par[2]*x
    
    val_x <- cbind(1, val_x)
    x <- cbind(1,x)
    
    info_matrix <- matrix(0, nrow = 2, ncol = 2)     # information matrix 
    info_matrix[1,1] <- sum(exp(z_p) / (1+exp(z_p))^2) /n
    info_matrix[1,2] <- sum(x[,2] * exp(z_p) / (1+exp(z_p))^2)/n
    info_matrix[2,1] <- sum(x[,2] * exp(z_p) / (1+exp(z_p))^2)/n
    info_matrix[2,2] <- sum(x[,2]^2 * exp(z_p) / (1+exp(z_p))^2) /n
    
    var_beta_plugin_matrix <- inv(info_matrix)/n
    var_beta_plugin <- var_beta_plugin_matrix[2,2]
    
    val_z1 <- beta_hat[1] + beta_hat[2]*val_x
    z2 <- gamma_hat[1] + gamma_hat[2]*x
    z3 <- gamma_bar[1] + gamma_bar[2]*x
    val_z2 <- z2[v]
    val_z3 <- z3[v]
  
    
    info_matrix1 <- matrix(0, nrow = 2, ncol = 2)     # information matrix 
    info_matrix1[1,1] <- sum(exp(val_z1) / (1+exp(val_z1))^2) /m
    info_matrix1[1,2] <- sum(val_x[,2] * exp(val_z1) / (1+exp(val_z1))^2) /m
    info_matrix1[2,1] <- sum(val_x[,2] * exp(val_z1) / (1+exp(val_z1))^2) /m
    info_matrix1[2,2] <- sum(val_x[,2]^2 * exp(val_z1) / (1+exp(val_z1))^2) /m
    
    info_matrix2 <- matrix(0, nrow = 2, ncol = 2)     # information matrix
    info_matrix2[1,1] <- sum(exp(z3) / (1+exp(z3))^2) /n
    info_matrix2[1,2] <- sum(x[,2] * exp(z3) / (1+exp(z3))^2) /n
    info_matrix2[2,1] <- sum(x[,2] * exp(z3) / (1+exp(z3))^2) /n
    info_matrix2[2,2] <- sum(x[,2]^2 * exp(z3) / (1+exp(z3))^2) /n

    rho = m/n   # proportion: size of validation set / size of full set
    
    
    p1 <- matrix(0, nrow = 2, ncol =2)
    p2 <- matrix(0, nrow = 2, ncol =2)
    p3 <- matrix(0, nrow = 2, ncol =2)
    
    for (i in 1:m){
      
      phi1 <- (val_y[i] - (exp(val_z1[i])/(1+exp(val_z1[i])))) %*% val_x[i,]
      phi1 <- inv(info_matrix1) %*% as.vector(phi1)
      p1 <- p1 + phi1 %*% t(phi1)
      
      phi2 <- (val_s[i] - (exp(val_z3[i])/(1+exp(val_z3[i])))) %*% val_x[i,]
      phi2 <- inv(info_matrix2) %*% as.vector(phi2)
      p2 <- p2 + phi2 %*% t(phi2)
      
      p3 <- p3 + phi1 %*% t(phi2)
      
    }
    
    sigma <- p1 / m
    sigma_star <- ((1-rho)) * p2 / m
    omega <- (1-rho) * p3 / m
    

    beta_bar = t(beta_hat -  omega %*% inv(sigma_star) %*% (gamma_hat - gamma_bar))
    
    # varainces
    var_beta = summary(glm0)$coefficients[2,2]**2
    var_beta_hat = summary(glm1)$coefficients[2,2]**2
    var_gamma_hat = summary(glm2)$coefficients[2,2]**2
    var_gamma_bar = summary(glm3)$coefficients[2,2]**2
    
    var_beta_bar = (sigma - omega * ((sigma_star)^(-1)) * omega) / m
    var_beta_bar = var_beta_bar[2,2]
    
    
    
    diff <- var_beta_bar - var_beta_hat  # difference between after bias reduction and before bias reduction
    
    myest[iter, ] = c(var_beta,var_beta_hat, var_gamma_bar,var_beta_plugin,var_beta_bar, diff)
    b1[iter,] = beta_hat
    b2[iter,] = beta_bar
    r1[iter,] = gamma_bar
    b_real[iter,] = beta_real
    b3[iter,] = beta_plugin
  
  }
  
  est <- data.frame(b1[,2],r1[,2],b3,b2[,2])
  names(est) = c("Model (1)", "Model (2)", "Model (3)","Model (4)")

  p <- ggplot(stack(est), aes(x=ind, y=values,fill=ind)) +
    geom_boxplot(notch=TRUE) +
    labs(title="SCENARIO 3",x=" ", y = "Value") +
    geom_hline(yintercept = 1,linetype="dotted") +
    theme_classic(base_size=30) +
    ylim(0.5,2)
  print(p)

  # deduction percentage
  reduction = (var(b1[,2]) - var(b2[,2])) / var(b1[,2])
  
  # bias
  bias <- mean(b_real[,2]) - beta1
  bias_b1 <- mean(b1[,2]) - beta1
  bias_r1 <- mean(r1[,2]) - beta1
  bias_b2 <- mean(b2[,2]) - beta1
  bias_b3 <- mean(b3) - beta1
  # bias_list <- list(bias,bias_b1, bias_r1, bias_b2)
  
  
  # empirical variance
  var_b_real <- var(b_real[,2])
  var_b1 <- var(b1[,2])
  var_r1 <- var(r1[,2])
  var_b2 <- var(b2[,2])
  var_b3 <- var(b3)
  # var_list <- list(var_b_real,var_b1,var_r1,var_b2)
  
  
  result <- c(bias_b1,bias_r1,bias_b3,bias_b2,
              sqrt(var_b1),sqrt(var_r1),sqrt(var_b3),sqrt(var_b2))
  
  return(list(result=result,reduction=reduction))
}

## result frame for differential misclassificaiton scenario 1
result_frame1 <- data.frame(matrix(NA,nrow = 18,ncol=8))

## result frame for differential misclassificaiton scenario 2
# result_frame2 <- data.frame(matrix(NA,nrow = 18,ncol=8))

## result frame for differential misclassificaiton scenario 3
# result_frame3 <- data.frame(matrix(NA,nrow = 18,ncol=8))
i <- 1
for (m in c(100,200,400)){
  for (beta0 in c(-0.5,-1,-1.5)){
    for (beta1 in c(1.0,0.25)){
      print(paste(i,m,beta0,beta1))
      temp <- simulation_plugin_diff(m=m, n=10000, niters=100,beta0 = beta0,beta1 = beta1)
      result1 <- temp$result
      reduction <- temp$reduction
      print(reduction*100)
      result_frame2[i,] <- result1
      i <- i + 1
    }
  }
}

save.image("diff_sce1.Rdata")

save.image("diff_sce2.Rdata")

save.image("diff_sce3.Rdata")
