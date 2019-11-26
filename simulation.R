##simulation A1
DataGeneratorA1 <- function(N,p1,p2,a0,a1,a2,a3,r,b0,b1,b2,b3,b4){
  
  M1 <- rbinom(N,1,p1)
  
  M2 <- rbinom(N,1,p2)
  
  U <- rnorm(N,0,1)
  
  e1 <- runif(N,0,1)
  
  e2 <- runif(N,0,1)
  G <- a0+a1*M1+a2*M2+r*M1*M2+a3*U+e1
  
  Y <- b0+b1*G+b2*M1+b3*M2+b4*U+e2
  datasimulation <- data.frame(M1,U,G,Y,M2)
  return(datasimulation)
}

##simulation A2
library(MASS)
DataGenerator <- function(N,p1,p2,a0,a1,a2,a3,r,b0,b1,b2,b3,b4,c,d){
  Sigma <- matrix(c(1,d,c,d,1,c,c,c,1),3,3)
  h <- mvrnorm(N, rep(0,3), Sigma)
  M1 <- h[,1]
  M1 <- ifelse(M1<=quantile(h[,1],p1),0,1)
  M2 <- h[,2]
  M2 <- ifelse(M2<=quantile(h[,1],p2),0,1)
  U <-h[,3] 
  e1 <- runif(N,0,1)
  e2 <- runif(N,0,1)
  G <- a0+a1*M1+a2*M2+a3*U+r*M1*M2+e1
  Y <- b0+b1*G+b2*M1+b3*M2+b4*U+e2
  
  datasimulation <- data.frame(M1,M2,G,Y,U)
  
  return(datasimulation)
}

##simulation B1
DataGenerator <- function(N,p1,p2,a0,a1,a2,a3,r,b0,b1,b2,b3,b4){
  
  M1 <- rbinom(N,1,p1)
  M2 <- rbinom(N,1,p2)
  U <- runif(N,0,1)
  e <- runif(N,0,1)
  
  pG <- a0+a1*M1+a2*M2+r*M1*M2+a3*U
  G <-rbinom(N,1,pG) 
  
  Y <- b0+b1*G+b2*M1+b3*M2+b4*U+e
  datasimulation <- data.frame(M1,U,G,Y,M2)
  return(datasimulation)
}

##simulation B2
library(MASS)
DataGenerator <- function(N,p1,p2,a0,a1,a2,a3,r,b0,b1,b2,b3,b4,c,d){
  Sigma <- matrix(c(1,d,c,d,1,c,c,c,1),3,3)
  h <- mvrnorm(N, rep(0,3), Sigma)
  M1 <- h[,1]
  M1 <- ifelse(M1<=quantile(h[,1],p1),0,1)
  M2 <- h[,2]
  M2 <- ifelse(M2<=quantile(h[,1],p2),0,1)
  UU <-h[,3] 
  U <- pnorm(UU)
  e <- runif(N,0,1)
  pG <- a0+a1*M1+a2*M2+r*M1*M2+a3*U
  G <- rbinom(N,1,pG)
  Y <- b0+b1*G+b2*M1+b3*M2+b4*U+e
  
  datasimulation <- data.frame(M1,M2,G,Y,U)
  
  return(datasimulation)
}


main <- function(NN,N,p1,p2,a0,a1,a2,a3,r,b0,b1,b2,b3,b4){  
  
  result <- NULL
  
  for(i in 1:NN){
    fdata <- DataGenerator(N,p1,p2,a0,a1,a2,a3,r,b0,b1,b2,b3,b4) 
    
    beta1 <- estimate(fdata,N)
    mbias1 <- beta1-b1
    sd <- sqrt((0.01+1/12)/(t(fdata$M1)%*%fdata$G*(1/t(fdata$M1)%*%fdata$M1)*(t(fdata$G)%*%fdata$M1*N)))
    
    mod2 <- lm(fdata$Y~fdata$G)
    beta2 <- mod2$coef[2]
    mbias2 <- mod2$coef[2]-b1
    
    
    mod3 <- lm(fdata$Y~fdata$G+fdata$M1+fdata$M2)
    beta3 <- mod3$coef[2]
    mbias3 <- mod3$coef[2]-b1
    
    result_once <- c(b1,beta1,mbias1,beta2,mbias2,beta3,mbias3,sd)
    result <- rbind(result,result_once)
  }
  
  end_result1 <- apply(result,2,mean)
  end_result2 <- data.frame(end_result1)
  se1 <- sd(result[,2]) 
  se2 <- sd(result[,4])
  se3 <- sd(result[,6])
  mse1 <- (result[,3])^2+se1^2
  mse2 <- (result[,5])^2+se2^2
  mse3 <-(result[,7])^2+se3^2
  end_result3 <- rbind(end_result2,se1,se2,se3,mse1,mse2,mse3) 
  end_result <- t(end_result3)
  return(end_result)
}

##
bias1_N <- NULL
bias2_N <- NULL
bias3_N <- NULL
se1_N <- NULL
se2_N <- NULL
se3_N <- NULL
mse1_N <- NULL
mse2_N <- NULL
mse3_N <- NULL
N <- seq(1000,10000,1000)
for(i in 1:length(N)){
  result <- main(NN=2000,N[i],p1=0.5,p2=0.5,a0=0,a1=0.6,a2=0.4,a3=0.2,r=0.6,b0=0,b1=1,b2=0.4,b3=0.3,b4=0.1)
  bias1_N[i] <- result[,3]
  se1_N[i] <- result[,8]
  bias2_N[i] <- result[,5]
  se2_N[i] <- result[,9]
  bias3_N[i] <- result[,7]
  se3_N[i] <- result[,10]
  mse1_N[i] <-result[,11] 
  mse2_N[i] <-result[,12] 
  mse3_N[i] <-result[,13] 
}
a0 <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
a1<- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
a3 <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
r <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
b0<- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
b1<- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
b2<- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
b4<- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
c <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
d <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)

bias1_a0<-rep(0,10)
bias2_a0<-rep(0,10)
bias3_a0<-rep(0,10)
se1_a0<-rep(0,10)
se2_a0<-rep(0,10)
se3_a0<-rep(0,10)
mse1_a0<-rep(0,10)
mse2_a0<-rep(0,10)
mse3_a0<-rep(0,10)
for(i in 1:10){
  result <- main(NN=2000,N=2000,p1=0.5,p2=0.5,a0[i],a1=0.6,a2=0.4,a3=0.2,r=0.6,b0=0,b1=1,b2=0.4,b3=0.3,b4=0.1)
  bias1_a0[i] <- result[,3]
  se1_a0[i] <- result[,8]
  bias2_a0[i] <- result[,5]
  se2_a0[i] <- result[,9]
  bias3_a0[i] <- result[,7]
  se3_a0[i] <- result[,10]
  mse1_a0[i] <-result[,11] 
  mse2_a0[i] <-result[,12] 
  mse3_a0[i] <-result[,13] 
}

bias1_a1<-rep(0,10)
bias2_a1<-rep(0,10)
bias3_a1<-rep(0,10)
se1_a1<-rep(0,10)
se2_a1<-rep(0,10)
se3_a1<-rep(0,10)
mse1_a1<-rep(0,10)
mse2_a1<-rep(0,10)
mse3_a1<-rep(0,10)
for(i in 1:10){
  result <- main(NN=2000,N=2000,p1=0.5,p2=0.5,a0=0,a1[i],a2=0.4,a3=0.2,r=0.6,b0=0,b1=1,b2=0.4,b3=0.3,b4=0.1)
  bias1_a1[i] <- result[,3]
  se1_a1[i] <- result[,8]
  bias2_a1[i] <- result[,5]
  se2_a1[i] <- result[,9]
  bias3_a1[i] <- result[,7]
  se3_a1[i] <- result[,10]
  mse1_a1[i] <-result[,11] 
  mse2_a1[i] <-result[,12] 
  mse3_a1[i] <-result[,13] 
}

bias1_a3<-rep(0,10)
bias2_a3<-rep(0,10)
bias3_a3<-rep(0,10)
se1_a3<-rep(0,10)
se2_a3<-rep(0,10)
se3_a3<-rep(0,10)
mse1_a3<-rep(0,10)
mse2_a3<-rep(0,10)
mse3_a3<-rep(0,10)
for(i in 1:10){
  result <- main(NN=2000,N=2000,p1=0.5,p2=0.5,a0=0,a1=0.6,a2=0.4,a3[i],r=0.6,b0=0,b1=1,b2=0.4,b3=0.3,b4=0.1)
  bias1_a3[i] <- result[,3]
  se1_a3[i] <- result[,8]
  bias2_a3[i] <- result[,5]
  se2_a3[i] <- result[,9]
  bias3_a3[i] <- result[,7]
  se3_a3[i] <- result[,10]
  mse1_a3[i] <-result[,11] 
  mse2_a3[i] <-result[,12] 
  mse3_a3[i] <-result[,13] 
}
bias1_b0<-rep(0,10)
bias2_b0<-rep(0,10)
bias3_b0<-rep(0,10)
se1_b0<-rep(0,10)
se2_b0<-rep(0,10)
se3_b0<-rep(0,10)
mse1_b0<-rep(0,10)
mse2_b0<-rep(0,10)
mse3_b0<-rep(0,10)
for(i in 1:10){
  result <- main(NN=2000,N=2000,p1=0.5,p2=0.5,a0=0,a1=0.6,a2=0.4,a3=0.2,r=0.6,b0[i],b1=1,b2=0.4,b3=0.3,b4=0.1)
  bias1_b0[i] <- result[,3]
  se1_b0[i] <- result[,8]
  bias2_b0[i] <- result[,5]
  se2_b0[i] <- result[,9]
  bias3_b0[i] <- result[,7]
  se3_b0[i] <- result[,10]
  mse1_b0[i] <-result[,11] 
  mse2_b0[i] <-result[,12] 
  mse3_b0[i] <-result[,13] 
}

bias1_b1<-rep(0,10)
bias2_b1<-rep(0,10)
bias3_b1<-rep(0,10)
se1_b1<-rep(0,10)
se2_b1<-rep(0,10)
se3_b1<-rep(0,10)
mse1_b1<-rep(0,10)
mse2_b1<-rep(0,10)
mse3_b1<-rep(0,10)
for(i in 1:10){
  result <- main(NN=2000,N=2000,p1=0.5,p2=0.5,a0=0,a1=0.6,a2=0.4,a3=0.2,r=0.6,b0=0,b1[i],b2=0.4,b3=0.3,b4=0.1)
  bias1_b1[i] <- result[,3]
  se1_b1[i] <- result[,8]
  bias2_b1[i] <- result[,5]
  se2_b1[i] <- result[,9]
  bias3_b1[i] <- result[,7]
  se3_b1[i] <- result[,10]
  mse1_b1[i] <-result[,11] 
  mse2_b1[i] <-result[,12] 
  mse3_b1[i] <-result[,13] 
}

bias1_b2<-rep(0,10)
bias2_b2<-rep(0,10)
bias3_b2<-rep(0,10)
se1_b2<-rep(0,10)
se2_b2<-rep(0,10)
se3_b2<-rep(0,10)
mse1_b2<-rep(0,10)
mse2_b2<-rep(0,10)
mse3_b2<-rep(0,10)
for(i in 1:10){
  result <- main(NN=2000,N=2000,p1=0.5,p2=0.5,a0=0,a1=0.6,a2=0.4,a3=0.2,r=0.6,b0=0,b1=1,b2[i],b3=0.3,b4=0.1)
  bias1_b2[i] <- result[,3]
  se1_b2[i] <- result[,8]
  bias2_b2[i] <- result[,5]
  se2_b2[i] <- result[,9]
  bias3_b2[i] <- result[,7]
  se3_b2[i] <- result[,10]
  mse1_b2[i] <-result[,11] 
  mse2_b2[i] <-result[,12] 
  mse3_b2[i] <-result[,13] 
}

bias1_b4<-rep(0,10)
bias2_b4<-rep(0,10)
bias3_b4<-rep(0,10)
se1_b4<-rep(0,10)
se2_b4<-rep(0,10)
se3_b4<-rep(0,10)
mse1_b4<-rep(0,10)
mse2_b4<-rep(0,10)
mse3_b4<-rep(0,10)
for(i in 1:10){
  result <- main(NN=2000,N=2000,p1=0.5,p2=0.5,a0=0,a1=0.6,a2=0.4,a3=0.2,r=0.6,b0=0,b1=1,b2=0.4,b3=0.3,b4[i])
  bias1_b4[i] <- result[,3]
  se1_b4[i] <- result[,8]
  bias2_b4[i] <- result[,5]
  se2_b4[i] <- result[,9]
  bias3_b4[i] <- result[,7]
  se3_b4[i] <- result[,10]
  mse1_b4[i] <-result[,11] 
  mse2_b4[i] <-result[,12] 
  mse3_b4[i] <-result[,13] 
}

r <- seq(0.1,1,0.1)
r1=r2=r3=r4=r5=r
bias1_r<-rep(0,10)
bias2_r<-rep(0,10)
bias3_r<-rep(0,10)
se1_r<-rep(0,10)
se2_r<-rep(0,10)
se3_r<-rep(0,10)
mse1_r<-rep(0,10)
mse2_r<-rep(0,10)
mse3_r<-rep(0,10)
for(i in 1:10){
  result <- main(NN=2000,N=2000,p1=0.5,p2=0.2,p3=0.5,p4=0.1,a0=0,a1=0.6,a2=0.4,a3=0.2,a4=0.1,a5=0.2,r[i],r1[i],r2[i],r3[i],r4[i],r5[i],b0=0,b1=1,b2=0.4,b3=0.3,b4=0.1,b5=0.3,b6=0.1)
  bias1_r[i] <- result[,3]
  se1_r[i] <- result[,8]
  bias2_r[i] <- result[,5]
  se2_r[i] <- result[,9]
  bias3_r[i] <- result[,7]
  se3_r[i] <- result[,10]
  mse1_r[i] <-result[,11] 
  mse2_r[i] <-result[,12] 
  mse3_r[i] <-result[,13] 
}

bias1_c<-rep(0,10)
bias2_c<-rep(0,10)
bias3_c<-rep(0,10)
se1_c<-rep(0,10)
se2_c<-rep(0,10)
se3_c<-rep(0,10)
mse1_c<-rep(0,10)
mse2_c<-rep(0,10)
mse3_c<-rep(0,10)
for(i in 1:10){
  result <- main(NN=2000,N=2000,p1=0.5,p2=0.5,a0=0,a1=0.6,a2=0.4,a3=0.2,r=0.6,b0=0,b1=1,b2=0.4,b3=0.3,b4=0.1,c[i])
  bias1_c[i] <- result[,3]
  se1_c[i] <- result[,8]
  bias2_c[i] <- result[,5]
  se2_c[i] <- result[,9]
  bias3_c[i] <- result[,7]
  se3_c[i] <- result[,10]
  mse1_c[i] <-result[,11] 
  mse2_c[i] <-result[,12] 
  mse3_c[i] <-result[,13] 
}

bias1_d<-rep(0,10)
bias2_d<-rep(0,10)
bias3_d<-rep(0,10)
se1_d<-rep(0,10)
se2_d<-rep(0,10)
se3_d<-rep(0,10)
mse1_d<-rep(0,10)
mse2_d<-rep(0,10)
mse3_d<-rep(0,10)
for(i in 1:10){
  result <- main(NN=2000,N=2000,p1=0.5,p2=0.2,a0=0,a1=0.3,a2=0.4,a3=0.6,r=0.8,b0=0,b1=0.5,b2=0.1,b3=0.1,b4=0.1,c=0.4,d[i])
  bias1_d[i] <- result[,3]
  se1_d[i] <- result[,8]
  bias2_d[i] <- result[,5]
  se2_d[i] <- result[,9]
  bias3_d[i] <- result[,7]
  se3_d[i] <- result[,10]
  mse1_d[i] <-result[,11] 
  mse2_d[i] <-result[,12] 
  mse3_d[i] <-result[,13] 
}