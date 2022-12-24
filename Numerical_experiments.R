
# R code with numerical experiments from paper:
# Wesolowski J., Wieczorkowski R., Wójciak W., Optimality of the recursive Neyman allocation,
# Journal of Survey Statistics and Methodology (2022) 10, 1263–1275.

library(dplyr)
library(microbenchmark)
library(ggplot2)
library(stratification)


# R codes for used algorithms
source("CapacityScaling.R")
source("noptcond.R")
source("SimpleGreedy.R")
source("rNa.R")
source("SGa.R")
source("com_.R")

source("nopt_u.R") # new W. Wójciak codes


# random rounding
ran_round<-function(x)
{
  return( floor(x)+(runif(length(x))<(x-floor(x))) )
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




# generation of artificial populations

# seed for Nrep=100
#set.seed(2234)

# seed for Nrep=200
set.seed(876)

source("gen_population.R")
pop <- gen_population(Nrep=200)
Nh<-pop$Nh
Sh<-pop$Sh
NROW(Nh)
# plot(Nh*Sh)


# # population for fig.3
# Nh <- rep(1000,20)
# Sh <- 10^(1:20)
# Sh <- sample(Sh,20) 
# mh <- rep(100,20)
# Mh <- Nh
# dh <- Sh*Nh


# # ad. table 1.
# Nh <- rep(1000,20)
# Sh <- c(
# 0.33,2.65,0.15,0.66,0.15,15.45,1.49,1.74,0.30,0.93,2.37,0.36,0.14,0.37,
# 4.25,0.39,10.21,0.10,0.23,0.51
# )


dh <- Sh*Nh
mh<-rep(0,length(Nh)) # lower bounds
Mh <- Nh # upper bounds
(s1<-sum(mh))
(s2<-sum(Mh))
n <- 8000


# improving population to have allocation with values greater than 0
# using integer allocation algorithm
(n <- round(0.1*sum(Nh)))
(alc<-CapacityScaling(n, Nh, Sh, mh = mh, Mh = Mh))
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
NROW(Nh)



# variance comparison for rounding

# variance estimation for given allocation
varal <- function(Nh,Sh,nh)
{
  return( sum(Nh * (Nh - nh) * Sh*Sh / nh) )
}



tab <- NULL
for (f in seq(0.1,0.5,0.1)) {
  print(f)
  N<-sum(Nh)
  n<-round(f*N)
  
  if (s1<n && n<s2) {
    
    alc<-CapacityScaling(n, Nh, Sh, mh = mh, Mh = Mh)
    V0<-varal(Nh,Sh,alc)
    
    ix <- which(alc==1)
    nix <- sum(ix)
    
    nh_SGa <- SGa(n,Nh,Sh, Mh)
    v_SGa <- varal(Nh,Sh,nh_SGa)
    nh1_SGa <- round_oric(nh_SGa)
    
    v1_SGa <- varal(Nh,Sh,nh1_SGa)
    
    #comal <- coma(n,Nh,Sh, Mh, 0)
    #nh_coma <- comal$nh
    #v_coma <- comal$v
    nh_coma <- coma(n,Nh,Sh, Mh)
    v_coma <- varal(Nh,Sh,nh_coma)
    v1_coma <- varal(Nh,Sh,round_oric(nh_coma))
    
    
    nh_rNa <- rNa(n,Nh,Sh, Mh)
    v_rNa <-  varal(Nh,Sh,nh_rNa)
    v1_rNa <- varal(Nh,Sh,round_oric(nh_rNa))
    
    
    # nh_noptcond <- noptcond(Nh*Sh , mh , Mh , n)
    # v_noptcond <- varal(Nh,Sh,nh_noptcond)
    # v1_noptcond <- varal(Nh,Sh,round_oric(nh_noptcond))
    
    
    
    tabi <- data.frame(N=N,J=length(Nh),f=f,n=n,
                       rv_SGa=v_SGa/V0,rv_coma=v_coma/V0,
                       rv_rNa=v_rNa/V0,
                       #rv_noptcond=v_noptcond/V0,
                       rv1_SGa=v1_SGa/V0,
                       rv1_coma=v1_coma/V0,
                       rv1_rNa=v1_rNa/V0
                       #rv1_noptcond=v1_noptcond/V0
    )
    
    tab<-bind_rows(tab,tabi)
    
  }
}

#saveRDS(tab,"tab1.rds")

View(tab)

#tab <- readRDS("tab1.rds")


library(knitr)

# information about number of strata and population size
(info_strata_pop <- paste(length(Nh),"strata,  N=", sprintf("%3.1e",N)))


tab1 <- tab[,1:7]
colnames(tab1)<-c("N","J","fraction","n","V(SGa)/V0","V(coma)/V0","V(rNa)/V0")
tab_tex1 <- kable(tab1,format="latex",digits=6,align="r",
      caption=paste(info_strata_pop,"\n No rounding: ratio of variances for different allocations"))


tab2 <- tab[,c(1:4,8:10)]
colnames(tab2)<-c("N","J","fraction","n","V(SGa)/V0","V(coma)/V0","V(rNa)/V0")
tab_tex2 <- kable(tab2,format="latex",digits=6,align="r",
      caption=paste(info_strata_pop,"\n ORIC rounding: ratio of variances for different allocations"))


# ew. zapis do tablic w TEX
#cat(tab_tex1,file=paste0("tab1_compvar_",info_strata_pop,".tex"))
#cat(tab_tex2,file=paste0("tab2_compvar_",info_strata_pop,".tex"))




# time comparison
# creating data with times for selected algorithms and different fractions
tab<-NULL

for (f in seq(0.01,0.5,0.05)) {
print(f)
N<-sum(Nh)
n<-round(f*N)

if (s1<n && n<s2) {

alc<-CapacityScaling(n, Nh, Sh, mh = mh, Mh = Mh)
#V0<-alc$v

h_over <- sum(rNa_old(n,Nh,Sh, Mh)$nh>=Mh)  # number of take-all strata
recursive_Neyman0=round(rNa_old(n,Nh,Sh, Mh)$nh)
recursive_Neyman=round(rNa(n,Nh,Sh, Mh))

print(all(recursive_Neyman0==recursive_Neyman))

#sequential_al=round(coma(n,Nh,Sh, Mh)$nh)
sequential_al=round(coma(n,Nh,Sh, Mh))
print(all(recursive_Neyman0==sequential_al))

  
options(digits=3)

ex<-microbenchmark(times=100,unit="ms",
                   CapScal=CapacityScaling(n, Nh, Sh, mh = mh, Mh = Mh),
                   noptcond=noptcond(dh , mh , Mh , n),
                   #rNa=(rNa(n,Nh,Sh, Mh)$nh),
                   #coma=comaR(n,Nh,Sh, Mh, 1),
                   rNa=rNa(n,Nh,Sh, Mh),
                   coma=coma(n,Nh,Sh, Mh),
                   
                   #coma_R0=comaR(n,Nh,Sh, Mh, 0),
                   #coma_1=coma(n,Nh,Sh, Mh, 1),
                   #coma_0=coma(n,Nh,Sh, Mh, 0),
                   SGa=SGa(n,Nh,Sh, Mh, onebyone = TRUE)
)
summary(ex)
autoplot(ex)

                                            

exi<-group_by(ex,expr) %>% 
  summarise(Median_time=median(time)/1e6,
            Mean_time=mean(time)/1e6) # from nanoseconds to miliseconds
exi<-mutate(exi,N=N,f=f,H=length(Nh), hover=h_over)

# adding information about number of iterations for rNa algorithm
exi$niter <- rNa_old(n,Nh,Sh, Mh)$iter



tab<-bind_rows(tab,exi)

}
}

#saveRDS(tab,"tab2.rds")


tab1 <- readRDS("tab1.rds")
tab2 <- readRDS("tab2.rds")
#tab <- readRDS("tab3.rds")
#tab4 <- readRDS("tab4.rds")
tab <- bind_rows(tab1,tab2)



# creation of plots 

library(ggplot2)
#library(ggpubr)
library(ggrepel)

options(digits=6)


# version of plots based on WW code
source("to_facet_grid_layout.R")

tab1 <- readRDS("tab1.rds")
tab2 <- readRDS("tab2.rds")
#tab <- readRDS("tab3.rds")
#tab4 <- readRDS("tab4.rds")
tab <- bind_rows(tab1,tab2)
tab <- filter(tab,expr!='CapScal',expr!='noptcond')

tab<-mutate(tab,algorithm=expr)
tab <- mutate(tab,algorithm=relevel(algorithm,"rNa"))
count(tab,algorithm)

df <- as.data.frame(tab)
df <- tab

# preapre population column, and algorithms facet
population <- paste0(df$H, " strata, N = ", df$N) # population info (no of strata, pop size)
population <- factor(population, levels = unique(population)[order(df$H)]) # levels order influences plotting order

df <- data.frame(facet = "algorithms",
                 plyr::rename(df[, c("expr", "Median_time", "f", "niter", "hover")], c('expr' = 'series', 'Median_time' = 'value')),
                 population = population,
                 takeall_pct = (df$hover/df$H )*100) # Take-strata in %

# prepare takeall facet
df_takeall <- df[, c("f", "population", "takeall_pct", "niter", "hover")]
df_takeall <- df_takeall[!duplicated(df_takeall), ]
df_takeall <- data.frame(facet = "takeall",
                         series = "takeall_pct",
                         plyr::rename(df_takeall, c('takeall_pct' = 'value')))

# rbind algorithms and takeall facets
df_plot <- plyr::rbind.fill(df[, colnames(df_takeall)], df_takeall)


## plot it
p <- 
  ggplot(data = df_plot, mapping = aes(x = f, y = value)) +
  geom_line(data = subset(df_plot, facet == "algorithms"), mapping = aes(linetype = series)) +
  geom_point(data = subset(df_plot, facet == "algorithms"), mapping = aes(shape = series), size = 2) +
  #geom_text(data = subset(df_plot, facet == "algorithms" & series == "rNa"), aes(label = niter), vjust = -.6) +
  ggrepel::geom_text_repel(data = subset(df_plot, facet == "algorithms" & series == "rNa"), aes(label = niter)) +
  geom_bar(data = subset(df_plot, facet == "takeall"), mapping = aes(y = value), stat = "identity") +
  geom_text(data = subset(df_plot, facet == "takeall"), mapping = aes(label = hover), size = 3.0, vjust = 1.2, color="white") +
  #geom_blank(data = subset(df_plot, facet == "takeall"), mapping = aes(y = 100)) +
  theme(panel.background = element_rect(fill = NA, colour = "black"),
        panel.grid = element_line(colour = "grey"),
        #strip.text.y = element_text(size = 8),
        strip.background.y = element_blank()) +
  #plot.title = element_text(hjust = 0.5)) +
  labs(shape = "Algorithm", linetype = "Algorithm", x = "Sample fraction", 
       title = "Time comparison of selected algorithms") +
  facet_wrap(~ facet + population, nrow = 2, scales = "free_y",  
             labeller = labeller(facet = c('algorithms' = 'Time [milisec]', 'takeall' = 'Take-all \nstrata [%]'))) 
  #scale_y_continuous(breaks = scales::pretty_breaks(5), expand = expand_scale(mult = c(0.05, 0.06), add = c(0, 13)))  
  
p

gt <- to_facet_grid_layout(p,panel_resize_factor = 4)
grid::grid.draw(gt)




