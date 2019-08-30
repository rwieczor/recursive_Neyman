
#' @title Integer-valued optimal univariate allocation under constraints for stratified sampling
#' @description Better algorithm from paper Friedrich et al. (2015) for integer-valued optimal
#'   allocation in stratified sampling
#'
#' @param v0 - upper limit for value of variance which must be attained for computed optimal allocation
#' @param Nh - population sizes in strata
#' @param Sh - standard deviations for given variable in strata
#' @param mh - lower constraints for sample sizes in strata
#' @param Mh - upper constraints for sample sizes in strata
#'
#' @return list with: vector of optimal allocation sizes and value of variance for obtained allocation
#'
#' @references Ulf Friedrich, Ralf Münnich, Sven de Vries, Matthias Wagner,
#' Fast integer-valued algorithms for optimal allocations under constraints in stratified sampling,
#' Computational Statistics and Data Analysis 92 (2015), 1-12
#'
#' @export

CapacityScaling2 <- function(v0, Nh, Sh,
                             mh = rep(1, length(Nh)),
                             Mh = Nh)
{
  nh <- mh
  
  if (any(nh > Mh)) {
    stop("There are no feasible solutions")
  }
  
  dh <- Nh *Sh
  ah <- dh/sum(dh)
  n <- sum((Sh/ah)*Sh*Nh*Nh)/(v0+sum(Nh*Sh*Sh))
  n<-round(n)

  vmin <- sum(Nh * (Nh - nh) * Sh^2 / Mh)
  if (vmin>v0) {
    cat("Target variance is too high for lower constraints !","\n")
    return(list(nh=Mh,v=vmin))   
  }
  
  H<-length(Nh)
  s<-ceiling(n/(2*H))
  #cat("n s",n,s,"\n")
  
  while (!is.na(s) && s>1) {
    r <- 0L
    v <- sum(Nh * (Nh - nh) * Sh^2 / nh)
    
    while (v > v0 && sum(nh) < sum(Nh)) {
      r <- r + 1
      Vh <- Nh * Sh / sqrt(nh * (nh+1)) * (nh+1 <= Mh)
      h <- which.max(Vh)
      if (nh[h]+1<=Mh[h]) {
        if (nh[h]+s<=Mh[h])  nh[h]<-nh[h]+s
        else nh[h]<-nh[h]+1
      }
      v <- sum(Nh * (Nh - nh) * Sh^2 / nh)
    }
    
    nh <- pmax(nh-s,mh)
    s <- ceiling(s/2)
  }
  
  nh <- SimpleGreedy2(v0,Nh,Sh,mh,Mh,nh)$nh
  v <- sum(Nh * (Nh - nh) * Sh^2 / nh)
  
  return(list(nh = nh, v = v))
}

