

# Why low coverage (77%) even with eta = 1?
# Persists with simpler simulation setup


# 
# ######### FOR CLUSTER USE #########
# 
# # because Sherlock 2.0 restores previous workspace
# rm( list = ls() )
# 
# # load command line arguments
# args = commandArgs(trailingOnly = TRUE)
# jobname = args[1]
# scen = args[2]  # this will be a letter
# 
# # get scen parameters created by genSbatch
# setwd("/home/groups/manishad/SAPB")
# scen.params = read.csv( "scen_params.csv" )
# p = scen.params[ scen.params$scen.name == scen, ]
# 
# print(p)
# 
# 
# # simulation reps to run within this job
# # this need to match n.reps.in.doParallel in the genSbatch script
# sim.reps = 200
# 
# 
# # EDITED FOR C++ ISSUE WITH PACKAGE INSTALLATION
# library(crayon, lib.loc = "/home/groups/manishad/Rpackages/") # for dplyr
# library(dplyr, lib.loc = "/home/groups/manishad/Rpackages/")
# library(foreach, lib.loc = "/home/groups/manishad/Rpackages/")
# library(doParallel, lib.loc = "/home/groups/manishad/Rpackages/")
# library(boot, lib.loc = "/home/groups/manishad/Rpackages/")
# library(metafor, lib.loc = "/home/groups/manishad/Rpackages/")
# library(data.table, lib.loc = "/home/groups/manishad/Rpackages/")
# library(robumeta, lib.loc = "/home/groups/manishad/Rpackages/")
# library(metafor, lib.loc = "/home/groups/manishad/Rpackages/")
# library(fansi, lib.loc = "/home/groups/manishad/Rpackages/")  # for dplyr
# library(utf8, lib.loc = "/home/groups/manishad/Rpackages/")  # for dplyr
# library(cli, lib.loc = "/home/groups/manishad/Rpackages/")  # for dplyr
# library(nleqslv, lib.loc = "/home/groups/manishad/Rpackages/")
# 
# 
# # for use in ml load R
# # install.packages( c("doParallel", "foreach", "mvtnorm", "StepwiseTest", "matrixcalc"), lib = "/home/groups/manishad/Rpackages/" )
# 
# source("helper_sim_study_SAPB.R")
# source("weightfunct from package.R")
# 
# # set the number of cores
# registerDoParallel(cores=16)
# ######### END OF CLUSTER PART #########


######### FOR LOCAL USE #########


# REMEMBER TO INCREASE BOOT REPS IF NEEDED! :)
# note: total studies is k * per.cluster
k = c(20)  # number of clusters prior to selection
per.cluster = c(5)  # studies per cluster
mu = c(0.2)  # RE distribution mean
V = c(1)  # RE heterogeneity
V.gam = c(0)  # variance of random intercepts (can't be > V because that's total heterogeneity!)
sei.min = 1  # runif lower bound for study SEs
sei.max = c(1.5)  # runif upper bound for study SEs
# sei.min = .1  # runif lower bound for study SEs   ~ REDUCE SEs
# sei.max = c(.3)  # runif upper bound for study SEs
eta = c(1)  # selection prob
q = c(2)
boot.reps = c(0)
bt.meta.model = c("rma.uni")
bt.type = c("wtd.vanilla")
orig.meta.model = rev( c("robumeta.lazy") )
true.dist = "exp" # ~~~ CHANGED
SE.corr = FALSE
select.SE = FALSE


# matrix of scenario parameters
scen.params = expand.grid(k,
                          per.cluster,
                          mu,
                          V,
                          V.gam,
                          sei.min,
                          sei.max,
                          eta,
                          q,
                          boot.reps,
                          orig.meta.model,
                          bt.meta.model,
                          bt.type,
                          true.dist,
                          SE.corr,
                          select.SE )

names(scen.params) = c("k",
                       "per.cluster",
                       "mu",
                       "V",
                       "V.gam",
                       "sei.min",
                       "sei.max",
                       "eta",
                       "q",
                       "boot.reps",
                       "orig.meta.model",
                       "bt.meta.model",
                       "bt.type",
                       "true.dist",
                       "SE.corr",
                       "select.SE")


# remove scenarios with V = 0 but V.gam > 0
scen.params = scen.params[ !(scen.params$V == 0 & scen.params$V.gam > 0), ]

# remove scenarios with nonsense combinations of methods and parameters
stupid = rep( FALSE, nrow(scen.params) )
stupid[ scen.params$V > 0 & scen.params$orig.meta.model == "fixed" ] = TRUE
scen.params = scen.params[ stupid == FALSE, ]

# meta-analysis choices: fixed, wtd.score, robumeta, robumeta.lazy
# (last one uses a lazier initial guess for t2 to avoid calling wtd.score method)
# # when V.gam is 0, always use Vevea
# scen.params$orig.meta.model = c("robumeta")
# scen.params$orig.meta.model[ scen.params$V.gam == 0 &
#                                scen.params$V > 0 ] = "wtd.score"
# scen.params$orig.meta.model[ scen.params$V == 0 ] = "fixed"

# compute true P
scen.params$trueP = 1 - pnorm( ( scen.params$q - scen.params$mu ) / sqrt(scen.params$V) )


# scenario names
start.at = 1
scen.params$scen.name = paste( "scen.",
                               start.at : ( start.at + nrow(scen.params) - 1 ),
                               sep="" )
( n.scen = length(scen.params[,1]) )



# # first make parameters csv file in genSbatch
# library(here)
# setwd( here:::here("Code/Simulation study") )
# ( scen.params = read.csv("scen_params.csv") )

# WAS FOR LOCAL USE
library(foreach)
library(doParallel)
library(dplyr)
library(boot)
library(robumeta)
library(metafor)
library(nleqslv)

sim.reps = 500

library(here)
setwd(here:::here("Code"))
source("helper_sim_study_SAPB.R")
source("weightfunct from package.R")

# set the number of cores
registerDoParallel(cores=8)

scen = "scen.1"

######### END OF LOCAL PART #########


########################### THIS SCRIPT COMPLETELY RUNS 1 SIMULATION  ###########################


# j is the number of simulation iterations to run sequentially
# so for j=10, we are generating 10 observed datasets,
# each with reps original datasets
#  and 500 bootstrap iterates for each


# for ( scen in scen.params$scen.name ) {  # LOCAL ONLY
#  cat("\n\nscen ", scen)  # LOCAL ONLY
rm(agg)
rep.time = system.time({
  
  rs = foreach( i = 1:sim.reps,
                .combine=rbind,
                .errorhandling = "remove"  # this shouldn't happen
  ) %dopar% {
    
    # extract simulation params for this scenario (row)
    # exclude the column with the scenario name itself (col) 
    p = scen.params[ scen.params$scen.name == scen, names(scen.params) != "scen.name"]
    
    boot.reps = p$boot.reps
    
    
    # ##### Simulate Dataset #####
    # # make sure there's at least 1 nonsignificant study
    # n.nonsig = 0
    # while( n.nonsig == 0 ) {
    #   d = sim_data2(p)
    #   n.nonsig = sum(d$pval > 0.05 | d$yi < 0)
    #   SE.corr.emp = d$SE.corr.emp[1]
    #   SE.mean.emp = d$SE.mean.emp[1]
    #   
    #   P.select.SE.emp = d$P.select.SE.emp[1]
    #   P.publish.emp = d$P.publish.emp[1]
    # }
    
    d = sim_data2(p)
    #d = sim_data_DEBUG(p)
    
    # dim(d)
    # prop.table( table(d$weight) )
    

    # unlike the above, this one uses a lazier initial guess for tau^2
    if ( p$orig.meta.model == "robumeta.lazy" ) {
      
      # initialize t2
      re = rma.uni( yi = d$yi,
                      vi = d$vi )
      t2.guess = re$tau2

      
      meta.naive = robu( yi ~ 1,
                         data = d,
                         studynum = cluster,
                         
                         #studynum = 1:nrow(d),
                         userweights = 1 / (d$vi + t2.guess),
                         var.eff.size = vi,
                         small = TRUE )
      
      # # ~~~ NEW TRY: no userweights, no clustering
      # meta.naive = robu( yi ~ 1, 
      #                    data = d, 
      #                    studynum = cluster,
      #                    var.eff.size = vi,
      #                    small = TRUE )
      
      muhat.naive = meta.naive$b.r
      mu.se.naive = meta.naive$reg_table[["SE"]]  # used for Phat later
      t2.naive = NA
      mu.lo.naive = meta.naive$reg_table$CI.L
      mu.hi.naive = meta.naive$reg_table$CI.U
      t2.lo.naive = NA
      t2.hi.naive = NA
      
      Phat.naive = list( lo = NA, hi = NA )
    }
    

    ##### Write Results #####
    # fill in new row of summary dataframe with bias, coverage, and CI width for DM and bootstrap
    # dataframe with 3 rows, one for each method
    
    rows =     data.frame( 
      # BELOW IS FOR MULTIPLE METHODS
      # method of calculating CI: exponentiate logit or not?
      Method = c( "Naive" ),
      
      # stats for mean estimate
      # note that both boot CI methods use the same point estimate
      MuEst = c( muhat.naive ),  
      
      MuCover = c( covers( p$mu, mu.lo.naive, mu.hi.naive ) ),
      
      MuLo = c( mu.lo.naive ),
      
      MuHi = c( mu.hi.naive )
      
      # # observed sample size (after selection)
      # k.obs = rep( nrow(d), 2 ),
      # 
      # # number of nonsignificant studies before selection
      # n.nonsig = rep( sum(d$pval > 0.05), 2 ),
      # 
      # # FE mean in NS studies
      # nonsig.fe.mean = rep( nonsig.fe.mean,
      #                       2 ),
      # 
      # # empirical correlation between SE and true effects
      # SE.corr.emp = rep( d$SE.corr.emp[1], 
      #                    2),
      # 
      # # empirical correlation between SE and true effects
      # SE.mean.emp = rep( d$SE.mean.emp[1], 
      #                    2),
      # 
      # # E[Fi]
      # P.select.SE.emp = rep( d$P.select.SE.emp[1], 
      #                        2),
      # 
      # # E[Di]
      # P.publish.emp = rep( d$P.publish.emp[1], 
      #                      2)
    )
    
    # add in scenario parameters
    rows$scen.name = scen
    rows = as.data.frame( merge(scen.params, rows) )
    rows
    
  }  ### end foreach loop
  
} )[3]  # end timer


#}  # ends loop over scens


sum(is.na(rs$MuCover))
mean(rs$MuCover, na.rm = TRUE)





nrow(rs)
print(head(rs))

print(table(rs$scen.name))


# add rep time and CI width
rs$doParallel.min = rep.time/60
rs$MuCIWidth = rs$MuHi - rs$MuLo
mean(rs$MuCIWidth, na.rm=TRUE)

# LOCAL ONLY
my.vars = c("MuEst", "MuCover",
            "MuLo", "MuHi")

# mean performance among non-NA reps
options(scipen=999)
agg = rs %>% group_by(Method) %>% summarise_at( vars(my.vars), mean, na.rm = TRUE )
View(agg)

rs %>% group_by(Method) %>% summarise( median(MuEst) )
round( summary(rs$MuEst, na.rm = TRUE), 3 )

hist(rs$MuEst[ !is.na(rs$MuEst)])

# # how often did bootstrap fail?
# prop.miss = rs %>% group_by(Method) %>% summarise( PropMiss = mean( is.na(MuCover) ) )
# agg$PropMiss = prop.miss$PropMiss


agg$scen.name = scen
( agg = merge( scen.params, agg ) )


########################### REWRITE FROM SCRATCH ###########################


# THIS ONE WORKS with regular call to robu (i.e., no userweights)! NOMINAL COVERAGE. ???
# With t2.guess and userweights: still works
# With sim_data2: still works! 


# REMEMBER TO INCREASE BOOT REPS IF NEEDED! :)
# note: total studies is k * per.cluster
k = c(20*5)  # number of clusters prior to selection
per.cluster = c(1)  # studies per cluster
mu = c(0.2)  # RE distribution mean
V = c(1)  # RE heterogeneity
V.gam = c(0)  # variance of random intercepts (can't be > V because that's total heterogeneity!)
sei.min = 1  # runif lower bound for study SEs
sei.max = c(1.5)  # runif upper bound for study SEs
eta = c(1)  # selection prob
q = c(1)
true.dist = "exp" 
SE.corr = FALSE
select.SE = FALSE


# matrix of scenario parameters
scen.params = expand.grid(k,
                          per.cluster,
                          mu,
                          V,
                          V.gam,
                          sei.min,
                          sei.max,
                          eta,
                          q,
                          true.dist,
                          SE.corr,
                          select.SE )

names(scen.params) = c("k",
                       "per.cluster",
                       "mu",
                       "V",
                       "V.gam",
                       "sei.min",
                       "sei.max",
                       "eta",
                       "q",
                       "true.dist",
                       "SE.corr",
                       "select.SE")



rm(rs)
rs = foreach( i = 1:sim.reps,
              .combine=rbind,
              .errorhandling = "remove"  # this shouldn't happen
) %dopar% {
  
  # # simulate my own exponential data - WORKS
  # #N = p$k * p$per.cluster
  # N = 20*5
  # 
  # mui = rexp( n = N,
  #             rate = 1^(-1/2) )
  # # shift to have desired mean
  # mui = (mui - 1^(-1/2)) + 0.2
  # 
  # # individual study SEs
  # sei = runif( n = N,
  #              min = 1,
  #              max = 1.5 )
  # vi = sei^2
  # 
  # # individual study point estimates
  # yi = rnorm( n = N, mean = mui, sd = sei )
  # 
  # dat = data.frame(yi, 
  #                  vi)
  
  # simulate with existing fn
  dat = sim_data2( p = data.frame(k = 20*5,
                                  per.cluster = 1,
                                  mu = 0.2,
                                  V = 1,
                                  V.gam = 0,
                                  sei.min = 1,
                                  sei.max = 1.5,
                                  eta = 1,
                                  q = 1,
                                  true.dist = "exp",
                                  SE.corr = FALSE,
                                  select.SE = FALSE ) )
  
  
  # initialize t2
  re = rma.uni( yi = dat$yi,
                vi = dat$vi )
  t2.guess = re$tau2
  
  meta.naive = robu( yi ~ 1,
                     data = dat,
                     #studynum = cluster,
                     
                     studynum = 1:nrow(dat),
                     userweights = 1 / (dat$vi + t2.guess),
                     var.eff.size = vi,
                     small = TRUE )
  
  return( data.frame(MuEst = meta.naive$b.r,
                     MuLo = meta.naive$reg_table$CI.L,
                     MuHi = meta.naive$reg_table$CI.U,
                     MuCover = (meta.naive$reg_table$CI.L < 0.2) & (meta.naive$reg_table$CI.U > 0.2)
                     ) )
  
}  ### end foreach loop


dim(rs)

any(is.na(rs$MuCover))
mean(rs$MuCover)

any(is.na(rs$MuEst))
mean(rs$MuEst)






























