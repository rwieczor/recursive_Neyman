#' @title Optimal sample allocation under upper constraints in stratified random sampling scheme
#' 
#' @description 
#' Maximize: \deqn{ \sum_{h \in J} \frac{N_h^2 S_h^2}{n_h} - \sum_{h \in J} N_h S_h^2 }
#' Subject to: \deqn{ \sum_{h \in J} = n } and \deqn{ n_h \leq u_h \forall h \in J }
#'
#' @param n - total sample size
#' @param N - strata sizes
#' @param S - standard deviations of a given variable in strata
#' @param M - upper constraints for sample sizes in strata
#'
#' @return vector of optimal sample allocations in strata
#'
#' @references Wojciak W. (2019), Optimal allocation in stratified sampling schemes,
#'   Master's diploma thesis, Warsaw University of Technology.
#'
#' @export
#' @examples
#'
#' SGaR(n = 190, N = c(3000, 4000, 5000, 2000), S = rep(1, 4), M = c(100, 90, 70, 80)) # 40.71429 54.28571 67.85714 27.14286
#' SGaR(n = 270, N = c(3000, 4000, 5000, 2000), S = rep(1, 4), M = c(100, 90, 70, 80)) # 66.66667 88.88889 70.00000 44.44444
#' SGaR(n = 300, N = c(3000, 4000, 5000, 2000), S = rep(1, 4), M = c(100, 90, 70, 80)) # 84 90 70 56
#' SGaR(n = 330, N = c(3000, 4000, 5000, 2000), S = rep(1, 4), M = c(100, 90, 70, 80)) # 100  90  70  70
#' SGaR(n = 340, N = c(3000, 4000, 5000, 2000), S = rep(1, 4), M = c(100, 90, 70, 80)) # 100  90  70  80
#' SGaR(n = 350, N = c(3000, 4000, 5000, 2000), S = rep(1, 4), M = c(100, 90, 70, 80)) # 100  90  70  80 
#'
#'
SGa <- function(n, N, S, M = N, onebyone = TRUE) {
  
  if( n >= sum(M) )
    return(M)
  
  H <- length(N)
  d <- N * S
  dsum <- sum(d)
  dM <- d / M
  J <- order(dM, decreasing = TRUE) # set of strata indices sorted by d/M in non-increasing order
  # it preserves the relative order of elements with equivalent values
  
  if(onebyone) {
    
    for (k in 1:H) {
      h <- J[k]
      ksi <- n / dsum
      if( ksi * dM[h] < 1) 
        break
      else {
        n <- n - M[h]
        dsum <- dsum - d[h]
      }
    }

  } else {
    
    k <- 1
    i <- 1
    repeat {
      ksi <- n / dsum 
      if ( ksi * dM[J[i]] >= 1) {
        i <- i + 1
      } else if (k != i) {
        h <- J[k:(i-1)] 
        n <- n - sum(M[h])
        dsum <- dsum - sum(d[h])
        k <- i
      } else 
        break
    }

  }
      
  if(k == 1) # ksi * d < M for all strata
    nopt <-  ksi * d
  else {
    nopt <- M
    J <- J[k:H]
    nopt[J] <- ksi * d[J]
  }
  
  return(nopt)
  
}

# Change of monotonicity algorithm from chapter 2.3.3 of Wojciak Master's diploma thesis (2019)
coma <- function(n, N, S, M = N) {
  
  if( n >= sum(M) )
    return(M)
  
  H <- length(N)
  d <- N * S
  dsum <- sum(d)
  J <- order(d / M, decreasing = TRUE) # set of strata indices sorted by d/M in non-increasing order
  # it preserves the relative order of elements with equivalent values
  
  ksi <- n / dsum
  
  for (k in 1:(H-1)) {
    
    h <- J[k]
    n <- n - M[h]
    dsum <- dsum - d[h]
    ksi_1 <- n / dsum

    if (ksi > ksi_1) # if ksi_k > ksi_{k+1} - change of monotonicity found
      break
    else
      ksi <- ksi_1
    
  }
  
  if(k == 1) # ksi * d < M for all strata
    nopt <-  ksi * d
  else {
    nopt <- M
    J <- J[ifelse(k == H-1 && ksi <= ksi_1, H, k):H] 
    nopt[J] <- ksi * d[J]
  }
  
  return(nopt)
  
}

# Recursive Neyman
rNa <- function(n, N, S, M = N) {
  
  if( n >= sum(M) )
    return(M)
  
  H <- length(N)
  J <- 1:H # strata ids
  d <- N * S
  dsum <- sum(d)
  
  nopt_ <- ( n / dsum ) * d # Neyman allocation
  i <- which(nopt_ >= M)
  
  repeat {
    
    if(length(i) == 0)
      break
    else {
      h <- J[i] # h - strata ids (unique) for which nopt_ >= M 
      J <- J[-i] # remove strata h from set J
      n <- n - sum(M[h])
      dsum <- dsum - sum(d[h])
      nopt_ <- ( n / dsum ) * d[J] # Neyman allocation
      i <- which(nopt_ >= M[J])
      
    }
    
  }
  
  if(length(J) == H) # nopt_ < M for all strata
    nopt <- nopt_
  else {
    nopt <- M
    nopt[J] <- nopt_
  }
  
  return(nopt)
  
}
