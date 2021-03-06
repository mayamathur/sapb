
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
q = c(2)
#q = c(2.2, 1.2)
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
rm(rs)
rep.time = system.time({

  rs = foreach( i = 1:sim.reps,
                .combine=rbind,
                .errorhandling = "remove"  # this shouldn't happen
                ) %dopar% {
                  
    # extract simulation params for this scenario (row)
    # exclude the column with the scenario name itself (col) 
    p = scen.params[ scen.params$scen.name == scen, names(scen.params) != "scen.name"]
    
    boot.reps = p$boot.reps
    
    
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
    
  
    
    ##### Get T2 and Phat Via IPCW Bootstrap #####
    # purpose is to get Phat and tau^2 since robumeta can't help with that
    if ( boot.reps > 0 ) {
      
      if ( p$bt.type == "wtd.cluster" ) {
        
        t2b = rep(NA, boot.reps)
        Mb = rep(NA, boot.reps)
        
        for ( i in 1:boot.reps ) {
          # stage 1: resample clusters, leaving observations intact
          #  no weights
          #https://stats.stackexchange.com/questions/46821/bootstrapping-hierarchical-multilevel-data-resampling-clusters
          cluster.ids = data.frame(cluster = sample(d$cluster, replace = TRUE))
          b = d %>% inner_join(cluster.ids, by = 'cluster')
          
          # stage 2: resample observations with individual weights, ignoring clusters
          b = sample_n( b, size = nrow(b), replace = TRUE, weight = b$weight )
          
          # meta-analyze the bootstraps
          if ( p$bt.meta.model == "rma.uni" ) {
            
            #tryCatch( {
            mb = rma.uni( yi = b$yi,
                          vi = b$vi,
                          knha = TRUE,
                          method = "REML",
                          control=list(maxiter=500, stepadj=0.5) )
            
            Mb[i] = mb$b
            t2b[i] = mb$tau2
            # }, error = function(err) {
            #   stop("asdfasdf")
            #   Mb[i] = "rma.uni error"
            #   t2b[i] = "rma.uni error"
            # } )
            
          }
          
          if ( p$bt.meta.model == "robumeta" ) {
            # clustered version
            # version that's "naive" wrt pub bias, but accounts for clustering
            mb = robu( yi ~ 1, 
                       data = b, 
                       var.eff.size = vi,
                       small = TRUE )
            Mb[i] = mb$b.r
            t2b[i] = mb$mod_info$tau.sq  # this works because no userweights
          }
          
        }  # end loop over boot.reps iterates
        
        # mean estimates (using IPCW)
        t2hat.bt = mean( t2b )
        muhat.bt = mean( Mb ) # no longer necessary because we're IPCW-weighting the naive models
        t2.se.bt = sd(t2b)
        
        # in case there is some weird problem with 
        #  computing the boot CI
        tryCatch( {
          percCIs = list( quantile( t2b, c(.025, 0.975) ), 
                          quantile( Mb, c(.025, 0.975) ) )
        }, error = function(err) {
          percCIs <<- list( c(NA, NA), c(NA, NA) )
        } )
        
      }  # end wtd cluster bootstrap
      
      
      if ( p$bt.type == "wtd.vanilla" ) {
        # NOT clustered
        boot.res = boot( data = d, 
                         parallel = "multicore",
                         R = boot.reps, 
                         weights = d$weight,
                         statistic = function(original, indices) {
                           
                           b = original[indices,]
                           
                           if ( p$bt.meta.model == "rma.uni" ) {
                             # meta-analyze the bootstraps
                             tryCatch( {
                               mb = rma.uni( yi = b$yi,
                                             vi = b$vi,
                                             knha = TRUE,
                                             method = "REML",
                                             control=list(maxiter=500, stepadj=0.5) )
                               
                               Mb = mb$b
                               t2b = mb$tau2
                             }, error = function(err) {
                               Mb = NA
                               t2b = NA
                             } )
                           }
                           
                           if ( p$bt.meta.model == "robumeta" ) {
                             # clustered version
                             # version that's "naive" wrt pub bias, but accounts for clustering
                             mb = robu( yi ~ 1, 
                                        data = b, 
                                        var.eff.size = vi,
                                        small = TRUE )
                             Mb = mb$b.r
                             t2b = mb$mod_info$tau.sq  # this works because no userweights
                           }
                           
                           c(t2b, Mb)
                         }
        )  # end call to boot()
        
        # TEST ONLY
        mean( boot.res$t[,1] )
        head(boot.res$t)
        
        # mean estimates (using IPCW)
        t2hat.bt = mean( boot.res$t[,1] )
        muhat.bt = mean( boot.res$t[,2] ) # no longer necessary because we're IPCW-weighting the naive models
        t2.se.bt = sd(boot.res$t[,1])
        
        # in case there is some weird problem with 
        #  computing the boot CI
        tryCatch( {
          percCIs = get_boot_CIs(boot.res, "perc", n.ests = 2)
        }, error = function(err) {
          percCIs <<- list( c(NA, NA), c(NA, NA) )
        } )
        
      }  # end wtd vanilla bootstrap
      
    }  # end loop for boot.reps > 0
    

    if ( boot.reps == 0 ) {
      percCIs = list( c(NA, NA), c(NA, NA), c(NA, NA) )
      muhat.bt = NA
      t2hat.bt = NA
      Phat.bt = data.frame( Est = NA, 
                            lo = NA, 
                            hi = NA )
    }
    
    
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
    
    
    if ( p$orig.meta.model == "rma.uni" ) {
      # of course, with eta != 1, this won't work
      # and the SEs/coverage will be anticonservative with strong clustering
      meta.naive = rma.uni( yi = d$yi,
                     vi = d$vi,
                     knha = TRUE,
                     weights = d$weight * 1 / (d$vi + t2hat.bt),
                     method = "REML", 
                     control=list(maxiter=500, stepadj=0.5) )
      
      muhat.naive = meta.naive$b
      mu.se.naive = meta.naive$se  # used for Phat later
      t2.naive = meta.naive$tau2
      mu.lo.naive = meta.naive$ci.lb
      mu.hi.naive = meta.naive$ci.ub
      t2.lo.naive = tau_CI(meta.naive)[1]
      t2.hi.naive = tau_CI(meta.naive)[2]
      
      Phat.naive = list( lo = NA, hi = NA )
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
      # if we did bootstrapping to initialize t2hat, use that
      if ( boot.reps > 0 ) t2.guess = t2hat.bt
      
      # if we haven't bootstrapped, fit the lazy, biased rma.uni model
      #  to guess at t2
      if ( boot.reps == 0 ) {
        re = rma.uni( yi = d$yi,
                      vi = d$vi )
        t2.guess = re$tau2
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
    
    
    # sanity check on weighting
    # prop.table( table(d$weight) )
    
    
   
    ##### Compute Phat.bt #####
    # use point estimate and SE from naive model, 
    # but bootstrapped t2 and SE
    # but note that we aren't bootstrapping Phat itself
    if ( boot.reps > 0 ) {
      Phat.bt = suppressWarnings( prop_stronger2( q = p$q,
                                                  M = muhat.naive,
                                                  t2 = t2hat.bt,
                                                  se.M = mu.se.naive,
                                                  se.t2 = t2.se.bt,
                                                  tail = "above",
                                                  boot = "never",
                                                  always.analytical = TRUE ) )
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
                           # BELOW IS FOR MULTIPLE METHODS
                           # method of calculating CI: exponentiate logit or not?
                           Method = c( "Naive",
                                       "PercBT" ),
                           
                           # stats for mean estimate
                           # note that both boot CI methods use the same point estimate
                           MuEst = c( muhat.naive,
                                      muhat.bt ),  
                           
                           MuCover = c( covers( p$mu, mu.lo.naive, mu.hi.naive ),
                                        covers( p$mu, percCIs[[2]][1], percCIs[[2]][2] ) ),

                           MuLo = c( mu.lo.naive,
                                     percCIs[[2]][1] ),
                           
                           MuHi = c( mu.hi.naive,
                                     percCIs[[2]][2] ),
                           
                           # stats for t2 estimate
                           T2Est = c( t2.naive,
                                      t2hat.bt ),
                           
                           T2Cover = c( covers( p$V, t2.lo.naive, t2.hi.naive ),
                                        covers( p$V, percCIs[[1]][1], percCIs[[1]][2] ) ),
                           
                           T2Lo = c( t2.lo.naive,
                                     percCIs[[1]][1] ),
                           
                           T2Hi = c( t2.hi.naive,
                                     percCIs[[1]][2] ),
                           
                           # stats for Phat estimate
                           # the Phat.naive will only be filled in for the Vevea model
                           PEst = c( Phat.naive$Est,  
                                      Phat.bt$Est ),
                           
                           PCover = c( covers( p$trueP, Phat.naive$lo, Phat.naive$hi ),
                                       covers( p$trueP, Phat.bt$lo, Phat.bt$hi ) ),
                           
                           PLo = c( Phat.naive$lo,
                                     Phat.bt$lo ),
                           
                           PHi = c( Phat.bt$hi,
                                    Phat.bt$hi ),
                           
                           # observed sample size (after selection)
                           k.obs = rep( nrow(d), 2 ),
                           
                           # number of nonsignificant studies before selection
                           n.nonsig = rep( sum(d$pval > 0.05), 2 ),
                           
                           # FE mean in NS studies
                           nonsig.fe.mean = rep( nonsig.fe.mean,
                                             2 ),
                           
                           # empirical correlation between SE and true effects
                           SE.corr.emp = rep( d$SE.corr.emp[1], 
                                              2),
                           
                           # empirical correlation between SE and true effects
                           SE.mean.emp = rep( d$SE.mean.emp[1], 
                                              2),
                           
                           # E[Fi]
                           P.select.SE.emp = rep( d$P.select.SE.emp[1], 
                                              2),
                           
                           # E[Di]
                           P.publish.emp = rep( d$P.publish.emp[1], 
                                              2),
                  
                           # sanity check for normality of true effects
                           true.effect.shapiro.pval = rep( shapiro.test(d$mui)$p.value, 2 )
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




