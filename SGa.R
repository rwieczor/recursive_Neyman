#'
#' @title Optimal univariate allocation under upper constraints for stratified sampling
#' @description Algorithm from paper: Stenger, Gabler (2005), 
#' implementation based on Algorithm from chapter 2.3.1 
#' of Wójciak Master's diploma thesis (2019)
#' 
#'
#' @param n -  target sample size for allocation
#' @param Nh - population sizes in strata
#' @param Sh - standard deviations for given variable in strata
#' @param Mh - upper constraints for sample sizes in strata
#'
#' @return list with: vector of optimal allocation sizes and value of variance for obtained allocation
#'
#' @references  Horst Stenger, Siegfried Gabler (2005), 
#' Combining random sampling and census strategies - Justification of inclusion probabilities equal to 1,
#' Metrika 61, 137-156,
#'  Wójciak W. (2019), Optimal allocation in stratified sampling schemes,
#'   Master's diploma thesis, Warsaw University of Technology
#'   
#'
#' @export


SGa <- function(n, Nh, Sh, Mh=NULL)
{
  
  if (is.null(Mh)) Mh <- Nh
  
  dh <- Sh*Nh
  H <- length(Nh)
  jk <- order(dh/Mh, decreasing = TRUE)
  a <- n
  b <- sum(dh)
  
  for (k in 1:H) {
    if (a*dh[jk[k]]/(b*Mh[jk[k]]) >= 1) {
      a <- a - Mh[jk[k]] 
      b <- b - dh[jk[k]]
    }
    else { kstar <- k; break }
    
  }
  #cat('kstar = ',kstar,"\n")
  
  Jstar <- (1:H)[jk[kstar:H]]
  Jstarc <- setdiff(1:H,Jstar)
  
  nh <- numeric(H)
  nh[Jstarc] <- Mh[Jstarc]
  nh[Jstar] <- (a/b)*dh[Jstar]
  
  v <- sum(Nh * (Nh - nh) * Sh^2 / nh)
  
  return(list(nh = nh, v = v))
  
}

