#'
#' @title Optimal univariate allocation under upper constraints for stratified sampling
#' @description Change of monotonicity algorithm from chapter 2.3.3 of Wójciak Master's diploma thesis (2019)
#'
#' @param n -  target sample size for allocation
#' @param Nh - population sizes in strata
#' @param Sh - standard deviations for given variable in strata
#' @param Mh - upper constraints for sample sizes in strata
#'
#' @return list with: vector of optimal allocation sizes and value of variance for obtained allocation
#'
#' @references Wójciak W. (2019), Optimal allocation in stratified sampling schemes,
#'   Master's diploma thesis, Warsaw University of Technology
#'
#' @export


coma <- function(n, Nh, Sh, Mh=NULL)
{
  
  if (is.null(Mh)) Mh <- Nh
  
  dh <- Sh*Nh
  H <- length(Nh)
  jk <- order(dh/Mh, decreasing = TRUE)
  a <- n
  b <- sum(Nh*Sh)
  
  for (k in 1:(H-1)) {
    a1 <- a - Mh[jk[k]]
    b1 <- b - dh[jk[k]]
    
    if ((a1/b1) >= (a/b)) {
      a <- a1
      b <- b1
      if (k==(H-1)) kstar <- H
    }
    else { kstar <- k; break }
    
  }
  #cat('kstar = ',kstar,"\n")
  # cat((a/b)*Sh[jk[kstar]],"\n")
  #  
  Jstar <- (1:H)[jk[kstar:H]]
  Jstarc <- setdiff(1:H,Jstar)
  
  nh <- numeric(H)
  nh[Jstarc] <- Mh[Jstarc]
  nh[Jstar] <- (n-sum(Mh[Jstarc]))*dh[Jstar]/sum(dh[Jstar])

  
  v <- sum(Nh * (Nh - nh) * Sh^2 / nh)
  
  return(list(nh = nh, v = v))
  
}

