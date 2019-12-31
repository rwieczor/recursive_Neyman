#'
#' @title Optimal univariate allocation under upper constraints for stratified sampling
#' @description Change of monotonicity algorithm from chapter 2.3.3 of Wojciak Master's diploma thesis (2019)
#'
#' @param n -  total sample size
#' @param N - population sizes in strata
#' @param S - standard deviations of a given variable in strata
#' @param M - upper constraints for sample sizes in strata
#' @param sort_ - logical, when TRUE strata are sorted by N\*S/M in non-increasing order, otherwise, max(N\*S/M) is executed at every iteration.
#'
#' @return  vector of optimal allocation sizes 
#'
#' @references Wojciak W. (2019), Optimal allocation in stratified sampling schemes,
#'   Master's diploma thesis, Warsaw University of Technology.
#'
#' @export
#' @examples
#' 
#' comaR(n = 190, N = c(3000, 4000, 5000, 2000), S = rep(1, 4), M = c(100, 90, 70, 80)) # 40.71429 54.28571 67.85714 27.14286
#' comaR(n = 270, N = c(3000, 4000, 5000, 2000), S = rep(1, 4), M = c(100, 90, 70, 80)) # 66.66667 88.88889 70.00000 44.44444
#' comaR(n = 300, N = c(3000, 4000, 5000, 2000), S = rep(1, 4), M = c(100, 90, 70, 80)) # 84 90 70 56
#' comaR(n = 330, N = c(3000, 4000, 5000, 2000), S = rep(1, 4), M = c(100, 90, 70, 80)) # 100  90  70  70
#' 
#'  
comaR <- function(n, N, S, M = N, sort_ = TRUE) {
  
  H <- length(N)
  d <- S * N
  
  d_sum <- sum(d)
  ksi <- n / d_sum
  
  if(sort_) {
    
    J <- order(d/M, decreasing = TRUE) # set of strata indices sorted by d/M in non-increasing order
    
    for (k in 1:(H-1)) {
      
      n1 <- n - M[J[k]]
      d_sum1 <- d_sum - d[J[k]]
      ksi1 <- n1 / d_sum1
      
      if (ksi <= ksi1) { # if ksi_k <= \ksi_{k+1}, continue 
        n <- n1
        d_sum <- d_sum1
        ksi <- ksi1
      } else
        break
    }
    
    Jstar <- J[ifelse(k == H-1 && ksi <= ksi1, H, k):H]  # strata ids for which n* < M
    
  } else {
    
    JdM <- cbind(d/M, 1:H) # 2nd column - stratum id (it starts from 1)
    
    for (k in 1:(H-1)) {
      
      indmax <- which.max(JdM[, 1]) # index of maximum d/M
      
      n1 <- n - M[JdM[indmax, 2]] # JdM[imax, 2] stratum id of maximum d/M
      d_sum1 <- d_sum - d[JdM[indmax, 2]]
      ksi1 <- n1 / d_sum1
      
      if (ksi <= ksi1) { # if ksi_k <= \ksi_{k+1}, continue 
        n <- n1
        d_sum <- d_sum1
        ksi <- ksi1
        JdM <- JdM[-indmax, , drop = FALSE]
      }
      else
        break
    }
    
    Jstar <- JdM[, 2] # strata ids for which n* < M
    
  }

  nopt <- M  
  nopt[Jstar] <- ksi * d[Jstar]
  return(nopt)
  
  #v <- sum(N[Jstar] * (N[Jstar] - nopt) * S[Jstar]^2 / nopt)
  #return(list(nopt = nopt, v = v))
  
}