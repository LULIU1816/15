TBCestimate<-function(fdata,N){
  
  M1 <- fdata$M1
  M2 <- fdata$M2
  
  m00 <- mean(fdata$G[M1==0 & M2==0])
  m01 <- mean(fdata$G[M1==0 & M2==1])
  m10 <- mean(fdata$G[M1==1 & M2==0])
  m11 <- mean(fdata$G[M1==1 & M2==1])
  
  n00 <- mean(fdata$Y[M1==0 & M2==0])
  n01 <- mean(fdata$Y[M1==0 & M2==1])
  n10 <- mean(fdata$Y[M1==1 & M2==0])
  n11 <- mean(fdata$Y[M1==1 & M2==1])
  
  M_a <- matrix(c(1,1,1,1,0,0,1,1,0,1,0,1,m00,m01,m10,m11),4,4)
  M_b <- matrix(c(n00,n01,n10,n11),4,1)
  
  M_beta <- solve(M_a,M_b)  
  beta <-M_beta[4,1]
  
  return(c(beta))
  
}