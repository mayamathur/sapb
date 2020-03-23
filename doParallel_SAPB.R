
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

# rm(list=ls())
#
# # REMEMBER TO INCREASE BOOT REPS IF NEEDED! :)
# # note: total studies is k * per.cluster
# k = c(20, 40, 80, 200)  # number of clusters prior to selection
# per.cluster = c(5)  # studies per cluster
# mu = c(0.6, .8)  # RE distribution mean
# V = c(1, 0)  # RE heterogeneity
# V.gam = c(0, 0.5)  # variance of random intercepts (can't be > V because that's total heterogeneity!)
# sei.min = 1  # runif lower bound for study SEs
# sei.max = c(1.5)  # runif upper bound for study SEs
# eta = c(100, 50, 20, 10, 1)  # selection prob
# q = c(2, 1)
# #q = c(2.2, 1.2)
# boot.reps = c(0)
# bt.meta.model = c("rma.uni")
# bt.type = c("wtd.vanilla")
# orig.meta.model = rev( c("robumeta.lazy") )
# true.dist = "exp"
# SE.corr = TRUE
# select.SE = TRUE

# REMEMBER TO INCREASE BOOT REPS IF NEEDED! :)
# note: total studies is k * per.cluster
k = c(20)  # number of clusters prior to selection
per.cluster = c(5)  # studies per cluster
mu = c(0.2)  # RE distribution mean
V = c(1)  # RE heterogeneity
V.gam = c(0)  # variance of random intercepts (can't be > V because that's total heterogeneity!)
sei.min = 1  # runif lower bound for study SEs
sei.max = c(1.5)  # runif upper bound for study SEs
eta = c(1)  # selection prob


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
                     
               
                          orig.meta.model,
       
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
                
           
                       "orig.meta.model",
                    
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
rm(rs)
rep.time = system.time({

  rs = foreach( i = 1:sim.reps,
                .combine=rbind,
                .errorhandling = "remove"  # this shouldn't happen
                ) %dopar% {
                  
    # extract simulation params for this scenario (row)
    # exclude the column with the scenario name itself (col) 
    p = scen.params[ scen.params$scen.name == scen, names(scen.params) != "scen.name"]
    
    
    ##### Simulate Dataset #####
    # make sure there's at least 1 nonsignificant study
    n.nonsig = 0
    while( n.nonsig == 0 ) {
      d = sim_data2(p)
      n.nonsig = sum(d$pval > 0.05 | d$yi < 0)
      SE.corr.emp = d$SE.corr.emp[1]
      SE.mean.emp = d$SE.mean.emp[1]
      
      P.select.SE.emp = d$P.select.SE.emp[1]
      P.publish.emp = d$P.publish.emp[1]
    }

    # dim(d)
    # prop.table( table(d$weight) )
    
 
    
    ##### Corrected Meta-Analyses #####
    
    # fixed-effects model
    if ( p$orig.meta.model == "fixed" ) {
      
      meta.naive = correct_est_fe( yi = d$yi,
                                   vi = d$vi,
                                   eta = p$eta )
      muhat.naive = meta.naive$est.adj
      mu.se.naive = NA
      t2.naive = NA
      mu.lo.naive = meta.naive$lo.adj
      mu.hi.naive = meta.naive$hi.adj
      t2.se.naive = NA
      t2.lo.naive = NA
      t2.hi.naive = NA
      
      # for this specification, can also get analytic Phat
      Phat.naive = data.frame( Est = NA,
                               lo = NA,
                               hi = NA
                                )
    }

    
    # Vevea-Woods RE model
    if ( p$orig.meta.model == "vevea" ) {
      meta.naive = correct_est_re( effect = d$yi,
                                    v = d$vi,
                                    weights = c(1, 1/p$eta) )
      muhat.naive = meta.naive$est[ meta.naive$param == "b" ]
      mu.se.naive = meta.naive$se[ meta.naive$param == "b" ]  # used for Phat later
      t2.naive = meta.naive$est[ meta.naive$param == "t2" ]
      mu.lo.naive =  meta.naive$lo[ meta.naive$param == "b" ]
      mu.hi.naive = meta.naive$hi[ meta.naive$param == "b" ]
      t2.se.naive = meta.naive$se[ meta.naive$param == "t2" ]
      t2.lo.naive = meta.naive$lo[ meta.naive$param == "t2" ]
      t2.hi.naive = meta.naive$hi[ meta.naive$param == "t2" ]
      
      # for this specification, can also get analytic Phat
      Phat.naive = suppressWarnings( prop_stronger2( q = p$q,
                                                  M = muhat.naive,
                                                  t2 = t2.naive,
                                                  se.M = mu.se.naive,
                                                  se.t2 = t2.se.naive,
                                                  tail = "above",
                                                  boot = "never",
                                                  always.analytical = TRUE ) )
    }
    
    # my version of Vevea-Woods model
    if ( p$orig.meta.model == "wtd.score" ) {
      
      
      meta.naive = correct_est_re_score( yi = d$yi,
                                          vi = d$vi,
                                          weights = d$weight,
                                          start.ests = c(0,0) )
      
      muhat.naive = meta.naive$est[ meta.naive$param == "b" ]
      mu.se.naive = meta.naive$se[ meta.naive$param == "b" ]  # used for Phat later
      t2.naive = meta.naive$est[ meta.naive$param == "t2" ]
      mu.lo.naive =  meta.naive$lo[ meta.naive$param == "b" ]
      mu.hi.naive = meta.naive$hi[ meta.naive$param == "b" ]
      t2.se.naive = meta.naive$se[ meta.naive$param == "t2" ]
      t2.lo.naive = meta.naive$lo[ meta.naive$param == "t2" ]
      t2.hi.naive = meta.naive$hi[ meta.naive$param == "t2" ]
      
      if ( muhat.naive > 500 ) break
      
      # for this specification, can also get analytic Phat
      Phat.naive = suppressWarnings( prop_stronger2( q = p$q,
                                                     M = muhat.naive,
                                                     t2 = t2.naive,
                                                     se.M = mu.se.naive,
                                                     se.t2 = t2.se.naive,
                                                     tail = "above",
                                                     boot = "never",
                                                     always.analytical = TRUE ) )
    }
    
    
    if ( p$orig.meta.model == "robumeta" ) {
      # if we did bootstrapping to initialize t2hat, use that
      if ( boot.reps > 0 ) t2.guess = t2hat.bt
      
      # if we haven't bootstrapped, fit the wtd score RE model to guess tau^2
      if ( boot.reps == 0 ) {
        re = correct_est_re_score( yi = d$yi,
                                  vi = d$vi,
                                  weights = d$weight,
                                  start.ests = c(0,0) )
        t2.guess = re$est[ re$param == "t2" ]
      }
          
      meta.naive = robu( yi ~ 1, 
                         data = d, 
                         studynum = cluster,
                         userweights = d$weight / (d$vi + t2.guess),
                         #userweights = weight,
                         var.eff.size = vi,
                         small = TRUE )
      
      muhat.naive = meta.naive$b.r
      mu.se.naive = meta.naive$reg_table[["SE"]]  # used for Phat later
      t2.naive = NA
      mu.lo.naive = meta.naive$reg_table$CI.L
      mu.hi.naive = meta.naive$reg_table$CI.U
      t2.lo.naive = NA
      t2.hi.naive = NA
      
      Phat.naive = list( lo = NA, hi = NA )
    }

    # unlike the above, this one uses a lazier initial guess for tau^2
    if ( p$orig.meta.model == "robumeta.lazy" ) {
   
      # guess t2 by fitting the lazy, biased rma.uni model
      re = rma.uni( yi = d$yi,
                    vi = d$vi )
      t2.guess = re$tau2
      
      meta.naive = robu( yi ~ 1, 
                         data = d, 
                         studynum = cluster,
                         userweights = d$weight / (d$vi + t2.guess),
                         #userweights = weight,
                         var.eff.size = vi,
                         small = TRUE )
      
      muhat.naive = meta.naive$b.r
      mu.se.naive = meta.naive$reg_table[["SE"]]  # used for Phat later
      t2.naive = NA
      mu.lo.naive = meta.naive$reg_table$CI.L
      mu.hi.naive = meta.naive$reg_table$CI.U
      t2.lo.naive = NA
      t2.hi.naive = NA
      
      Phat.naive = list( lo = NA, hi = NA )
    }
    
    

    ##### Estimate FE mean in nonsignificant/negative studies #####
    
    NS = ( d$pval > 0.05 | d$yi < 0 )
    
    nonsig.fe.mean = rma.uni( yi = d$yi[ NS ],
                              vi = d$vi[ NS ],
                              method = "FE" )$b

    ##### Write Results #####
    # fill in new row of summary dataframe with bias, coverage, and CI width for DM and bootstrap
    # dataframe with 3 rows, one for each method

    rows =     data.frame( 
                           
                           Method = c( "Naive" ),
                           
                           # stats for mean estimate
                           # note that both boot CI methods use the same point estimate
                           MuEst = c( muhat.naive ),  
                           
                           MuCover = c( covers( p$mu, mu.lo.naive, mu.hi.naive ) ),

                           MuLo = c( mu.lo.naive ),
                           
                           MuHi = c( mu.hi.naive ),
                           
                           # stats for t2 estimate
                           T2Est = c( t2.naive ),
                           
                           T2Cover = c( covers( p$V, t2.lo.naive, t2.hi.naive ) ),
                           
                           T2Lo = c( t2.lo.naive ),
                           
                           T2Hi = c( t2.hi.naive ),
                          
                           
                           # observed sample size (after selection)
                           k.obs = nrow(d),
                           
                           # number of nonsignificant studies before selection
                           n.nonsig = sum(d$pval > 0.05),
                           
                           # FE mean in NS studies
                           nonsig.fe.mean = nonsig.fe.mean,
                           
                           # empirical correlation between SE and true effects
                           SE.corr.emp = d$SE.corr.emp[1],
                           
                           # empirical correlation between SE and true effects
                           SE.mean.emp = d$SE.mean.emp[1],
                           
                           # E[Fi]
                           P.select.SE.emp = d$P.select.SE.emp[1],
                           
                           # E[Di]
                           P.publish.emp = d$P.publish.emp[1],
                  
                           # sanity check for normality of true effects
                           true.effect.shapiro.pval = shapiro.test(d$mui)$p.value
    )
    
    # add in scenario parameters
    rows$scen.name = scen
    rows = as.data.frame( merge(scen.params, rows) )
    rows

    }  ### end foreach loop
  
} )[3]  # end timer


#}  # ends loop over scens


nrow(rs)
print(head(rs))

print(table(rs$scen.name))


# add rep time and CI width
rs$doParallel.min = rep.time/60
rs$MuCIWidth = rs$MuHi - rs$MuLo
mean(rs$MuCIWidth, na.rm=TRUE)


# LOCAL ONLY
mean(rs$MuCover, na.rm=TRUE)

my.vars = c("MuEst", "MuCover",
            "MuLo", "MuHi",
            "MuCIWidth",
            "T2Est", "T2Cover",
            "T2Lo", "T2Hi",
            "PEst", "PCover",
            "k.obs", "n.nonsig", "SE.mean.emp" )

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


########################### WRITE LONG RESULTS  ###########################

setwd("/home/groups/manishad/SAPB/results/long")
write.csv( rs, paste( "long_results", jobname, ".csv", sep="_" ) )




