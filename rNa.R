
#' @title Optimal univariate allocation under upper constraints for stratified sampling
#' @description Classical recursive Neyman algorithm for optimal allocation in 
#' stratified sampling with upper constraints
#'
#' @param n -  target sample size for allocation
#' @param Nh - population sizes in strata
#' @param Sh - standard deviations for given variable in strata
#' @param Mh - upper constraints for sample sizes in strata
#'
#' @return  vector of optimal allocation sizes,
#'  and number of iterations 
#'
#' @references Sarndal, Swensson and Wretman (1992), Model Assisted 
#' Survey Sampling
#'
#' @export


rNa_old <- function(n, Nh, Sh, Mh=NULL)
  # Mh - ograniczenia gorne
{
  if (is.null(Mh)) Mh <- Nh

  H <- length(Nh)
  
  M <- sum(Mh)
  dh <- Nh*Sh
  nh <- n * dh / sum(dh)
  
  iter <- 0
  len_ic <- H
  
  while (1) {
    iter <- iter + 1
    ic <- which(nh < Mh)
    if (length(ic)==len_ic) break
    len_ic <- length(ic)
    r <- (n - M + sum(Mh[ic]))/ sum(dh[ic])
    nh[ic] <- dh[ic] * r
  }
  nh <- pmin(nh,Mh)
    
  #v <- sum(Nh * (Nh - nh) * Sh^2 / nh)
  
  return(list(nh = nh, iter=iter))
  
}




