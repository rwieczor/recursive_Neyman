noptcond <- function(dh ,mh ,Mh ,n)
# function from article: Siegfried Gabler, Matthias Ganninger, Ralf Münnich, 
# Optimal allocation of the sample size to strata under box constraints, 
# Metrika (2012) 75:151–161.
{
  H <- length(dh) # Number of strata
  m <- sum(mh) # Total of lower bounds
  M <- sum(Mh) # Total of upper bounds
  U1 <- order(dh/mh ,decreasing=FALSE) # Ordered set equivalent to U(L_1)
  U3 <- order(dh/Mh ,decreasing=TRUE) # Ordered set equivalent to U(L_2)
  nopt <- n*dh/sum(dh) # Naiive allocation
  Hk <- 1:H # Index 1 to H
  
  R <- 0
  while(R<=H)
  {
    i <- 0
    while(i<=R)
    {
      s1 <- as.integer (0)
      s3 <- as.integer(H+1)
      if(i<R){s1 <- as.integer(U1 [1:(R-i)])}
      if(i >0){ s3 <- as.integer(U3[1:i])}
      if(!any(s1%in%s3))
      {
        noptU1 <- numeric(length(s1))
        noptU3 <- numeric(length(s3))
        if(i<R){ noptU1 <- mh[s1]}
        if(i >0){ noptU3 <- Mh[s3]}
        # Omit solutions which violate the requirement
        # sum(noptU1 ,noptU3 )<=n
        if(sum(noptU1 ,noptU3 )<=n)
        {
          if(i==0 & R==0){ s2 <- Hk}
          if(i==0 & R>0) {s2 <- Hk[-s1]}
          if(i>0 & R==i){s2 <- Hk[-s3]}
          if(i>0 & R>i) {s2 <- Hk[-c(s1 ,s3)]}
          noptU2 <- (n-sum(noptU1)-sum(noptU3 ))*dh[s2]/
            sum(dh[s2])
          # If all conditions are met , stop and return solution
          if(sum(noptU2 <mh[s2 ])==0 & sum(noptU2 >Mh[s2 ])==0)
          {
            if(i==0 & R==0){ nopt1 <- cbind(s2 ,noptU2 )}
            if(i==0 & R >0){ nopt1 <- cbind(c(s1 ,s2),c(noptU1 ,noptU2 ))}
            if(i>0 & R==i){ nopt1 <- cbind(c(s2 ,s3),c(noptU2 ,noptU3 ))}
            if(i>0 & R>i){ nopt1 <- cbind(c(s1 ,s2 ,s3),c(noptU1 ,noptU2 ,noptU3 ))}
            
            return(nopt1[order(nopt1 [ ,1]) ,2])
          }
        }
      }
      i <- i+1
      
    }
    R <- R+1
  }
}




# example

#Sh=c(1,2,3,4,5)
#Nh=c(10,20,20,40,10)
#mh=c(1,2,2,4,1)
#Mh=c(2.5,5,5,10,2.5)
#dh<-Sh*Nh
#n<-20

#alok<-noptcond(dh ,mh ,Mh ,n)
#print(alok)



  


