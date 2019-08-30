
# R code with numerical experiments from paper:
# Wesolowski J., Wieczorkowski R., Wójciak W. (2019), Recursive Neyman-type approach 
# is optimal under upper bounds on sample strata sizes.

library(dplyr)
library(microbenchmark)
library(ggplot2)
library(stratification)


# R codes for used algorithms
source("CapacityScaling.R")
source("SimpleGreedy.R")
source("CapacityScaling2.R")
source("SimpleGreedy2.R")
source("noptcond.R")
source("rNa.R")
source("SGa.R")
source("coma.R")



# random rounding
ran_round<-function(x)
{
  return( floor(x)+(runif(length(x))<(x-floor(x))) )
}


# random rounding with multinomial distribution
ran_round_mult<-function(x)
{
  n <- round(sum(x))
  x1 <- floor(x)
  xf <- x-x1
  if (all(xf==0)) return(x)
  else {
    p <- xf/sum(xf)
    n1 <- n-sum(x1)
    return( x1 + rmultinom(1,n1,p) )
  }
}



# rounding based on article:
# Rama Count, Massoud Heidari, Optimal rounding under integer constraints
# December 2014, arxiv
round_oric <- function(x)
{
  n <- round(sum(x))
  m <- floor(x)
  y <- x - m
  Ix <- sum(y)
  
  if (Ix==0) return(x)
  else {
    iy <- order(-y)
    u <- unique(y[iy])
    z <- integer(length(x))
    for (i in 1:length(u)) z[iy] <- z[iy] + (y[iy]==u[i])*i
    iy2 <- order(-y,z,-m)
    #m[iy][iy2][1:Ix] <- ceiling(x[iy][iy2][1:Ix])
    m[iy2][1:Ix] <- (m[iy2][1:Ix])+1
    return(m)
  }

}




# generation of artificial population

source("gen_population.R")
pop <- gen_population(Nrep=200)
Nh<-pop$Nh
Sh<-pop$Sh
NROW(Nh)
# plot(Nh*Sh)

dh <- Sh*Nh
mh<-rep(0,length(Nh)) # lower bounds
Mh <- Nh # upper bounds
(s1<-sum(mh))
(s2<-sum(Mh))


# improving population to have allocation with values greater than 0
# using integer allocation algorithm
(n <- 0.1*sum(Nh))
(alc<-CapacityScaling(n, Nh, Sh, mh = mh, Mh = Mh))
alc<-alc$nh
sum(alc)

length(alc[alc>2])/length(alc)
ix <- which(alc>2)
Nh <- Nh[ix]
Sh <- Sh[ix]
dh <- dh[ix]
mh <- mh[ix]
Mh <- Mh[ix]
(s1<-sum(mh))
(s2<-sum(Mh))




N <- sum(Nh) # population size

# additional function used in simulations;
# fraction of strata not meeting the condition: nh<=Mh
# for given sampling fraction  'x'
fover<-function(x) {
  N <- sum(Nh)
  return( 100* sum(round(x*N)*dh/sum(dh) > Mh)/length(Nh) )
  
}
fover <- Vectorize(fover)

# number of strata instead of percentages
hover<-function(x) {
  N <- sum(Nh)
  return( sum(round(x*N)*dh/sum(dh) > Mh) )
  
}
hover <- Vectorize(hover)

curve(fover(x),0.1,0.9,xlab="Sampling fraction",
      ylab=substitute(paste("Percent of strata with Neyman allocation    ",n[h]>N[h])))

curve(hover(x),0.1,0.9,xlab="Sampling fraction",
      ylab=substitute(paste("Number of strata with Neyman allocation    ",n[h]>N[h])))





# variance comparison for rounding

# variance estimation for given allocation
varal <- function(Nh,Sh,nh)
{
  return( sum(Nh * (Nh - nh) * Sh*Sh / nh) )
}



options(digits=6)

mh<-rep(0,length(Nh))
Mh <- Nh # simple variant
(s1<-sum(mh))
(s2<-sum(Mh))


tab <- NULL
for (f in seq(0.1,0.5,0.1)) {
  print(f)
  N<-sum(Nh)
  n<-round(f*N)
  
  if (s1<n && n<s2) {
    
    alc<-CapacityScaling(n, Nh, Sh, mh = mh, Mh = Mh)
    V0<-alc$v
    V0/varal(Nh,Sh,alc$nh)
    
    ix <- which(alc$nh==1)
    nix <- sum(ix)
    
    SGal <- SGa(n,Nh,Sh, Mh)
    nh_SGa <- SGal$nh
    v_SGa <- SGal$v
    nh1_SGa <- round_oric(SGal$nh)
    
    v1_SGa <- varal(Nh,Sh,round_oric(nh_SGa))
    
    comal <- coma(n,Nh,Sh, Mh)
    nh_coma <- comal$nh
    v_coma <- comal$v
    v1_coma <- varal(Nh,Sh,round_oric(nh_coma))
    
    
    rNal <- rNa(n,Nh,Sh, Mh)
    nh_rNa <- rNal$nh
    v_rNa <- rNal$v
    v1_rNa <- varal(Nh,Sh,round_oric(nh_rNa))
    
    
    nh_noptcond <- noptcond(Nh*Sh , mh , Mh , n)
    v_noptcond <- varal(Nh,Sh,nh_noptcond)
    v1_noptcond <- varal(Nh,Sh,round_oric(nh_noptcond))
    
    
    
    tabi <- data.frame(N=N,f=f,n=n,
                       rv_SGa=v_SGa/V0,rv_coma=v_coma/V0,
                       rv_rNa=v_rNa/V0,rv_noptcond=v_noptcond/V0,
                       rv1_SGa=v1_SGa/V0,rv1_coma=v1_coma/V0,
                       rv1_rNa=v1_rNa/V0,rv1_noptcond=v1_noptcond/V0
    )
    
    tab<-bind_rows(tab,tabi)
    
  }
}

#saveRDS(tab,"tab.rds")

View(tab)

library(knitr)

# information about number of strata and population size
(info_strata_pop <- paste(length(Nh),"strata,  N=", sprintf("%3.1e",N)))


tab1 <- tab[,1:7]
colnames(tab1)<-c("N","fraction","n","V(SGa)/V0","V(coma)/V0","V(rNa)/V0","V(noptcond)/V0")
tab_tex1 <- kable(tab1,format="latex",digits=6,align="r",
      caption=paste(info_strata_pop,"\n No rounding: ratio of variances for different allocations"))


tab2 <- tab[,c(1:3,8:11)]
colnames(tab2)<-c("N","fraction","n","V(SGa)/V0","V(coma)/V0","V(rNa)/V0","V(noptcond)/V0")
tab_tex2 <- kable(tab2,format="latex",digits=6,align="r",
      caption=paste(info_strata_pop,"\n ORIC rounding: ratio of variances for different allocations"))


cat(tab_tex1,file=paste0("tab1_compvar_",info_strata_pop,".tex"))
cat(tab_tex2,file=paste0("tab2_compvar_",info_strata_pop,".tex"))








# creating data with times for selected algorithms and different fractions
tab<-NULL

for (f in seq(0.01,0.5,0.05)) {
print(f)
N<-sum(Nh)
n<-round(f*N)

if (s1<n && n<s2) {

alc<-CapacityScaling(n, Nh, Sh, mh = mh, Mh = Mh)
#V0<-alc$v
  
recursive_Neyman0=round(rNa(n,Nh,Sh, Mh)$nh)
sequential_al=round(coma(n,Nh,Sh, Mh)$nh)
print(all(recursive_Neyman0==sequential_al))

  
options(digits=3)

ex<-microbenchmark(times=100,unit="ms",
                   #noptcond=ranroundm(noptcond(dh , mh , Mh , n),n),
                   rNa=(rNa(n,Nh,Sh, Mh)$nh),
                   coma=(coma(n,Nh,Sh, Mh)$nh),
                   SGa=(SGa(n,Nh,Sh, Mh)$nh)
)
summary(ex)
autoplot(ex)

                                            

exi<-group_by(ex,expr) %>% 
  summarise(Median_time=median(time)/1e6,
            Mean_time=mean(time)/1e6) # from nanoseconds to miliseconds
exi<-mutate(exi,N=N,f=f,H=length(Nh), hover=hover(f))

# adding information about number of iterations for rNa algorithm
exi$niter <- rNa(n,Nh,Sh, Mh)$iter



tab<-bind_rows(tab,exi)

}
}

saveRDS(tab,"tab2.rds")


tab1 <- readRDS("tab1.rds")
tab2 <- readRDS("tab2.rds")
tab <- bind_rows(tab1,tab2)



# creation of plots 

library(ggplot2)
library(ggpubr)
library(ggrepel)

options(digits=2)

#tab<-readRDS("tab.rds")

tab<-mutate(tab,H=as.factor(H),algorithm=expr , 
            flab=paste0(as.character(round(hover)),"(",niter,")")  )

as.vector(table(tab$N))

xN<-count(tab,N)$N

levels(tab$H)<-paste(levels(tab$H),"strata,  N=", sprintf("%3.1e",xN))
round(table(tab$N))




p<-
  ggplot(data=tab,aes(x=f,y=Median_time, color=algorithm)) +
  geom_point(size=2) +
  geom_line(data=tab,aes(x=f,y=Median_time)) +
  geom_text_repel(data=filter(tab,algorithm=="rNa"),aes(x=f,y=Median_time,label=flab)) + 
  facet_wrap(~H,scale="free") +
  #scale_x_continuous(breaks = seq(0.1,0.9,0.1)) +
  labs(x="sample fraction", y="Time [miliseconds]" ,color="Algorithms: ", 
       title="Time comparison of selected methods",
       subtitle = "using microbenchmark package from R") + 
  theme_bw(base_size=12) + theme(legend.position = "right")
  #theme(legend.position = "right",legend.text=element_text(size=rel(1.2)))
  #coord_flip()

p

##ggsave("fig_times.png",p,device="png", dpi=600, width = 8, height = 8/1.618)
ggsave("fig_times.pdf",p,device="pdf", dpi=600, width = 8, height = 8/1.618)







