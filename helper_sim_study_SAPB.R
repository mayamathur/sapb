

################################ FNs: FORMATTING ################################ 

# round while keeping trailing zeroes
my_round = function(x, digits) {
  formatC( round( x, digits ), format='f', digits=digits )
}

format_CI = function( lo, hi, digits ) {
  paste( "[", my_round( lo, digits ), ", ", my_round( hi, digits ), "]", sep="" )
}

format_sval = function( sval, digits ) {
  if( sval < 1 ) return("None")
  else return( my_round(sval, digits) )
}

format_pval = function(p) {
  if (p >= 0.01) return( my_round( p, 2 ) )
  if (p < 0.01 & p > 10^-5 ) return( formatC( p, format = "e", digits = 0 ) )
  if ( p < 10^-5 ) return("< 1e-05")
}

# check CI coverage
covers = function( truth, lo, hi ) {
  if ( is.na(lo) | is.na(hi) ) return(NA)
  return( (lo <= truth) & (hi >= truth) )
}


# gives CI for tau from meta-analysis fit in metafor
tau_CI = function( meta, z.to.r = FALSE ) {
  t2.lb = meta$tau2 - qnorm(1 - 0.05/2) * meta$se.tau2
  t2.ub = meta$tau2 + qnorm(1 - 0.05/2) * meta$se.tau2
  
  if ( t2.lb > 0 ) tau.lb = sqrt(t2.lb) else tau.lb = 0
  tau.ub = sqrt(t2.ub)
  
  if( z.to.r == TRUE ) {
    tau.lb = z_to_r(tau.lb)
    tau.ub = z_to_r(tau.ub)
  }
  
  if ( tau.lb < 0 ) tau.lb = 0
  
  return( c(tau.lb, tau.ub))
}


############################# FNs FOR MY VERSION OF VEVEA-WOODS (WEIGHTED SCORE)  #############################


# returns score fns evaluated at x = c(muhat, tau^2)
# individ.terms = should it return the vectors of score contributions?
get_score = function(x,
                 yi,
                 vi,
                 weights,
                 individ.terms = FALSE ) {
  
  # #### Unweighted version -- directly from Brockwell Appendix A #####
  # # partial derivative of LL wrt mu
  # F1 = sum( ( yi - x[1] ) / ( vi + x[2] ) )
  # 
  # # partial derivative of LL wrt T2
  # num = ( yi - x[1] )^2 - ( vi + x[2] )
  # denom = 2 * ( vi + x[2] )^2
  # F2 = sum( num/denom )
  
  ##### Weighted #####
  # partial derivative of LL wrt mu
  F1.contribs = weights * ( yi - x[1] ) / ( vi + x[2] )
  F1 = sum( F1.contribs )
  
  # partial derivative of LL wrt T2
  num = ( yi - x[1] )^2 - ( vi + x[2] )
  denom = 2 * ( vi + x[2] )^2
  F2.contribs = weights * (num/denom)
  F2 = sum( F2.contribs )
  
  if ( individ.terms == FALSE ) return( c(F1 = F1, F2 = F2) )
  if (individ.terms == TRUE ) return( matrix( c( F1.contribs, F2.contribs ),
                                              ncol = 2 ) )
} 

##### Weighted Jacobian #####
Jac = function(x,
               yi,
               vi,
               weights ) {
  J <- matrix(0,nrow=2,ncol=2)
  
  # unweighted version -- MATCHES METAFOR! :)
  J[1,1] = -sum( weights / (vi + x[2]) )
  J[1,2] = -sum( weights * ( yi - x[1] ) / ( vi + x[2] )^2 )
  J[2,1] = J[1,2]
  J[2,2] = sum( weights * ( 0.5*vi + 0.5*x[2] - yi^2 + 2*x[1]*yi - x[1]^2 ) /
                  ( vi + x[2] )^3 )
  
  J
}


# start.ests: where to start the optimizer's search
# weights: a vector
correct_est_re_score = function( yi,
                           vi,
                           weights,
                           start.ests = c(0,0) ) {

  library(nleqslv)
  solns = nleqslv(x = start.ests,
                  fn = function(x) get_score( x = x, 
                                     yi = yi, 
                                     vi = vi, 
                                     weights = weights ),
                  method = "Newton",
                  jac = function(x) Jac( x,
                                         yi = yi,
                                         vi = vi, 
                                         weights = weights ),
                  jacobian = TRUE,
                  global = "cline" )
  ( my.ests = solns$x )
  J.emp = solns$jac
  
  # termination code of 1 indicates algo converged
  if ( as.character(solns$termcd) != "1" ) {
    return( data.frame( 
      eta = max(weights),  # only return the one weight that isn't 1
      param = c("b", "t2"),
      est = c(NA, NA),
      se = c(NA, NA),
      lo = c(NA, NA),
      hi = c(NA, NA)
    )
    )
  }
  
  
  # inference
  # Wooldridge pg 13
  # get estimated score contributions of each observation
  ki = get_score( x = my.ests, 
                  yi = yi, 
                  vi = vi, 
                  weights = weights,
                  individ.terms = TRUE )
  B0 = t(ki) %*% ki   
  var = solve(-J.emp) %*% B0 %*% solve(-J.emp)
  
  # this is a 2 x 2 matrix:
  # suppress warnings about off-diagonal elements being NA
  ( se = suppressWarnings( sqrt(var) ) )
  
  bhat.adj.lo = my.ests[1] - se[1,1] * qnorm(0.975)
  bhat.adj.hi = my.ests[1] + se[1,1] * qnorm(0.975)
  
  t2.adj.lo = max( 0, my.ests[2] - se[2,2] * qnorm(0.975) )
  t2.adj.hi = max( 0, my.ests[2] + se[2,2] * qnorm(0.975) )
  
  # fix negative tau^2 estimates
  my.ests[2] = max( 0, my.ests[2] )
  
  return( data.frame( 
    eta = max(weights),  # only return the one weight that isn't 1
    param = c("b", "t2"),
    est = my.ests,
    se = c(se[1,1], se[2,2]),
    lo = c(bhat.adj.lo, t2.adj.lo),
    hi = c(bhat.adj.hi, t2.adj.hi)
  )
  )
}


################################ FN: FIXED-EFFECTS CORRECTED ESTIMATE ################################ 

# correct point estimate and CI for FE specification
correct_est_fe = function( yi,
                           vi,
                           eta,
                           selection.tails = 1 ) {
  
  #browser()
  
  # 2-sided p-value even if 1-tailed selection
  pval = 2 * ( 1 - pnorm( abs(yi) / sqrt(vi) ) )
  
  # significance indicator
  if ( selection.tails == 1 ) S = (pval < 0.05) & (yi > 0)
  if ( selection.tails == 2 ) S = (pval < 0.05)
  
  # get FE means by significance status
  library(dplyr)
  
  dat = data.frame( yi, vi, S )
  
  # FE mean and sum of weights stratified by S
  strat = dat %>% group_by(S) %>%
    summarise( nu = sum( 1 / vi ),
               ybar = sum( yi / vi ) )
  
  # components by significance status
  ybarN = strat$ybar[ strat$S == 0 ]
  ybarS = strat$ybar[ strat$S == 1 ]
  nuN = strat$nu[ strat$S == 0 ]
  nuS = strat$nu[ strat$S == 1 ]
  
  # corrected pooled point estimate
  est = ( eta * ybarN + ybarS ) / ( eta * nuN + nuS )
  
  # inference
  var = ( eta^2 * nuN + nuS ) / ( eta * nuN + nuS )^2
  lo = est - qnorm(.975) * sqrt(var)
  hi = est + qnorm(.975) * sqrt(var)
  
  return( list(est.adj = est,
               lo.adj = lo,
               hi.adj = hi ) )
}




############################# FN: GET BOOT CIs FOR A VECTOR OF ESTIMATES #############################


# list with first entry for b and second entry for t2
# n.ests: how many parameters were estimated?
get_boot_CIs = function(boot.res, type, n.ests) {
  bootCIs = lapply( 1:n.ests, function(x) boot.ci(boot.res, type = type, index = x) )
  
  # list with first entry for b and second entry for t2
  # the middle index "4" on the bootCIs accesses the stats vector
  # the final index chooses the CI lower (4) or upper (5) bound
  bootCIs = lapply( 1:n.ests, function(x) c( bootCIs[[x]][[4]][4],
                                             bootCIs[[x]][[4]][5] ) )
}


############################# FN: GENERATE CLUSTERED META-ANALYSIS DATA #############################


# p: row of parameters dataframe
# unlike first version of sim_data, this allows non-normal true effects and 
#  SEs correlated with true effects
# so takes two new parameters: 
# p$true.dist: "norm" or "expo"
# p$SE.corr: TRUE or FALSE
# keep.all.studies: should we return all studies, including the unpublished ones?

# unclustered exponential
sim_data_DEBUG = function(p) {
  N = p$k * p$per.cluster
  true.effect.var = p$V
  # set var using properties of exponential
  gam2i = rexp( n = p$k, rate = true.effect.var^(-1/2) )
  # shift to have mean of 0
  # use fact that var = mean^2 in exponential
  gam2i = gam2i - true.effect.var^(1/2)
  
  # individual study SEs
  sei = runif( n = N, min = p$sei.min, max = p$sei.max )
  vi = sei^2
  
  mui = p$mu + gam2i
  
  # individual study point estimates
  yi = rnorm( n = N, mean = mui, sd = sei )
  
  return(data.frame(yi, sei, vi, mui))
}


# THIS IS FROM THE PREVIOUS SIMULATIONS!!! WORKS
sim_data_old = function(p) {
  
  # # TEST ONLY
  # p = data.frame( k = 500,
  #                 per.cluster = 1,
  #                 mu = .5,
  #                 V = 1,
  #                 V.gam = 0,
  #                 sei.min = 1,
  #                 sei.max = 1.5,
  #                 eta = 1,
  #                 true.dist = "exp",
  #                 SE.corr = TRUE )
  
  N = p$k * p$per.cluster
  
  # generate cluster random intercepts
  # called "zeta" in paper
  # these are normal even when true effect dist is exponential
  gam1 = rnorm( n = p$k, mean = 0, sd = sqrt( p$V.gam ) )
  gam1i = rep( gam1, each = p$per.cluster )
  
  # generate individual-study random intrcepts
  # these are called "gamma" in the paper
  # these are either normal or exponential
  if ( p$true.dist == "norm" ) gam2i = rnorm( n = N, mean = 0, sd = sqrt( p$V - p$V.gam ) )
  
  if ( p$true.dist == "exp" ) {
    true.effect.var = p$V - p$V.gam
    # set var using properties of exponential
    gam2i = rexp( n = p$k, rate = true.effect.var^(-1/2) )
    # shift to have mean of 0
    # use fact that var = mean^2 in exponential
    gam2i = gam2i - true.effect.var^(1/2)
  }
  
  
  # individual study SEs
  sei = runif( n = N, min = p$sei.min, max = p$sei.max )
  
  # individual study means
  if ( p$SE.corr == TRUE ) {
    beta = 6
    mui = p$mu + beta * sei + gam1i + gam2i
    
    # recenter them to have desired mean of p$mu
    # because E[sei}]
    mui = mui - beta * mean(sei)
    
    cor(sei, mui)
  } else {
    mui = p$mu + gam1i + gam2i
  }
  
  # individual study point estimates
  yi = rnorm( n = N, mean = mui, sd = sei )
  
  browser()
  
  
  d = data.frame( cluster = rep(1:p$k, each = p$per.cluster),
                  Study.name = 1:N,
                  yi = yi,
                  sei = sei,
                  vi = sei^2,
                  pval = 2 * ( 1 - pnorm( abs(yi) / sei ) ),
                  SE.corr.emp = cor(sei, mui) )  # empirical correlation
  
  # 1-tailed publication bias
  signif = d$pval < 0.05 & d$yi > 0
  publish = rep( 1, nrow(d) )
  publish[ signif == FALSE ] = rbinom( n = sum(signif == FALSE), size = 1, prob = 1/p$eta )
  
  d$weight = 1
  d$weight[ signif == 0 ] = p$eta
  d = d[ publish == 1, ]
  
  return(d)
}


# NEW VERSION
sim_data2 = function(p,
                     keep.all.studies = FALSE) {

  # # TEST ONLY
  # p = data.frame( k = 500,
  #                 per.cluster = 1,
  #                 mu = .5,
  #                 V = 1,
  #                 V.gam = 0,
  #                 sei.min = 1,
  #                 sei.max = 1.5,
  #                 eta = 1,
  #                 true.dist = "exp",
  #                 SE.corr = TRUE )

  N = p$k * p$per.cluster

  # generate cluster random intercepts
  # called "zeta" in paper
  # these are normal even when true effect dist is exponential
  gam1 = rnorm( n = p$k, mean = 0, sd = sqrt( p$V.gam ) )
  gam1i = rep( gam1, each = p$per.cluster )

  # generate individual-study random intrcepts
  # these are called "gamma" in the paper
  # these are either normal or exponential
  if ( p$true.dist == "norm" ) gam2i = rnorm( n = N, mean = 0, sd = sqrt( p$V - p$V.gam ) )

    if ( p$true.dist == "exp" ) {
      true.effect.var = p$V - p$V.gam
      # set var using properties of exponential
      gam2i = rexp( n = p$k, rate = true.effect.var^(-1/2) )
      # shift to have mean of 0
      # use fact that var = mean^2 in exponential
      gam2i = gam2i - true.effect.var^(1/2)
    }


  # individual study SEs
  sei = runif( n = N, min = p$sei.min, max = p$sei.max )

  # individual study means
  if ( p$SE.corr == TRUE ) {
    beta = 6
    mui = p$mu + beta * sei + gam1i + gam2i

    # recenter them to have desired mean of p$mu
    # because E[sei}]
    mui = mui - beta * mean(sei)

    cor(sei, mui)
  } else {
    mui = p$mu + gam1i + gam2i
  }

  # individual study point estimates
  yi = rnorm( n = N, mean = mui, sd = sei )

  # NEW: make cluster variable
  # this way, if we generate (say) 5 studies per "cluster", but V.gam = 0 in order
  #  to effectively generate without clustering, the cluster variable we ultimately
  #  pass to robumeta correctly reflects the lack of clustering
  if ( p$V.gam == 0 ) cluster = 1:N else cluster = rep(1:p$k, each = p$per.cluster)
  
  d = data.frame( cluster = cluster,
                  Study.name = 1:N,
                  mui,
                  yi = yi,
                  sei = sei,
                  vi = sei^2,
                  pval = 2 * ( 1 - pnorm( abs(yi) / sei ) ),
                  SE.corr.emp = cor(sei, mui) )  # empirical correlation

  # 1-tailed publication bias
  signif = d$pval < 0.05 & d$yi > 0
  publish = rep( 1, nrow(d) )
  publish[ signif == FALSE ] = rbinom( n = sum(signif == FALSE), size = 1, prob = 1/p$eta )

  d$weight = 1
  d$weight[ signif == 0 ] = p$eta

  # optionally, also select for small SE
  # Fi as defined in Supplement
  # tertiles of empirical SEs prior to selection
  Fi = rep( 1, nrow(d) )  # so that it has no effect if we're not selecting based on SE
  if ( p$select.SE == TRUE ) {
    terts = p$sei.min + (p$sei.max - p$sei.min) * c(1/3, 2/3)
    # prob of inclusion based on SE
    P.Fi = rep(NA, nrow(d))
    P.Fi[ d$sei < terts[1] ] = 1
    #Fi[ d$sei >= terts[1] ] = 0
    P.Fi[ d$sei >= terts[1] & d$sei < terts[2] ] = .75
    P.Fi[ d$sei >= terts[2] ] = .5
    Fi = rbinom( n = nrow(d),
                 size = 1,
                 prob = P.Fi )
  }

  if ( keep.all.studies == FALSE ) d = d[ publish == 1 & Fi == 1, ]
  if ( keep.all.studies == TRUE ){
    # just add the indicators
    d$publish = publish
    d$Fi = Fi
  }

  # record the mean empirical SE (after publication bias) as sanity check on SE selection
  d$SE.mean.emp = mean(d$sei)
  d$P.select.SE.emp = mean(Fi)
  d$P.publish.emp = mean(publish)

  return(d)
}



# # p: row of parameters dataframe
sim_data = function(p) {

  N = p$k * p$per.cluster

  # generate cluster random intercepts
  gam1 = rnorm( n = p$k, mean = 0, sd = sqrt( p$V.gam ) )
  gam1i = rep( gam1, each = p$per.cluster )

  # generate individual-study random intercepts
  gam2i = rnorm( n = N, mean = 0, sd = sqrt( p$V - p$V.gam ) )

  # individual study means
  mui = p$mu + gam1i + gam2i
  sei = runif( n = N, min = p$sei.min, max = p$sei.max )
  yi = rnorm( n = N, mean = mui, sd = sei )

  d = data.frame( cluster = rep(1:p$k, each = p$per.cluster),
                  Study.name = 1:N,
                  yi = yi,
                  sei = sei,
                  vi = sei^2,
                  pval = 2 * ( 1 - pnorm( abs(yi) / sei ) ) )

  # 1-tailed publication bias
  signif = d$pval < 0.05 & d$yi > 0
  publish = rep( 1, nrow(d) )
  publish[ signif == FALSE ] = rbinom( n = sum(signif == FALSE), size = 1, prob = 1/p$eta )

  d$weight = 1
  d$weight[ signif == 0 ] = p$eta
  d = d[ publish == 1, ]

  return(d)
}

# ##### Sanity Check #####
# d = sim_data( data.frame( k = 5,
#                           per.cluster = 20,
#                           mu = .5,
#                           V = 1,
#                           V.gam = .5,
#                           sei.min = 0.01,
#                           sei.max = 0.01,
#                           eta = 1 ) )
#
# # check out the clusters
# # randomly choose 5 of them for visibility
# library(ggplot2)
# ggplot( data = d, aes(x = yi, fill = as.factor(cluster) ) ) +
#   theme_classic() +
#   geom_histogram()
# # looks good


########################### FN: STITCH RESULTS FILES ###########################

# given a folder path for results and a common beginning of file name of results files
#   written by separate workers in parallel, stitch results files together into a
#   single csv.

stitch_files = function(.results.singles.path, .results.stitched.write.path=.results.singles.path,
                        .name.prefix, .stitch.file.name="stitched_model_fit_results.csv") {
  
  # get list of all files in folder
  all.files = list.files(.results.singles.path, full.names=TRUE)
  
  # we only want the ones whose name includes .name.prefix
  keepers = all.files[ grep( .name.prefix, all.files ) ]
  
  # grab variable names from first file
  names = names( read.csv(keepers[1] )[-1] )
  
  # # DEBUGGING
  tables <- lapply(keepers, function(x) read.csv(x, header= TRUE) )
  #tables <- lapply(keepers, function(x) read.csv(x, quote = "", header= TRUE) )
  s <- do.call(rbind , tables)
  
  names(s) = names( read.csv(keepers[1], header= TRUE) )
  # 
  # bad = lapply( tables, function(df)
  #   any( grepl("[[:digit:]]", df$scen.name) ) )
  # which(bad == TRUE)
  # 
  # rows = lapply(tables, function(df) nrow(df))
  # table( unlist(rows) )
  # END DEBUGGING
  
  # slower(?) version
  # # initialize stitched dataframe
  # s = as.data.frame( matrix(nrow=1, ncol=length(names)) )
  # names(s) = names
  # 
  # # stitch the files
  # for ( i in 1:length(keepers) ) {
  #   new.chunk = read.csv(keepers[i])[,-1]
  #   if ( nchar( as.character(new.chunk$scen.name) ) > 5 ) browser()
  #   #print(nrow(new.chunk))
  #   s = rbind(s, new.chunk)
  # }
  # 
  if( is.na(s[1,1]) ) s = s[-1,]  # delete annoying NA row
  write.csv(s, paste(.results.stitched.write.path, .stitch.file.name, sep="/") )
  return(s)
}


################################ FN: COMPUTE PROPORTION OF EFFECTS STRONGER THAN THRESHOLD ################################ 

# q: threshold of scientific importance
# M: pooled point estimate of meta-analysis
# t2: estimated heterogeneity, tau^2
# se.M: estimated standard error of M
# se.t2: estimated standard error of t2
# CI.level: confidence interval level
# tail: "above" to compute proportion of effects above q; 
#  "below" for proportion below q

# VERSION THAT DOES BOOTSTRAPPING
# dat: dataframe of yi and vi (with those names) for bootstrapping
# R: bootstrap replicates 
# boot: "never" or "ifneeded"
# always.analytical: if FALSE, returns NAs for SE, lo, and hi if Phat is extreme
#  but boot is never. If TRUE, returns possibly bad analytical inference.
prop_stronger2 = function( dat = NA,
                           q,
                           M,
                           t2,
                           se.M = NA,
                           se.t2 = NA,
                           CI.level = 0.95,
                           tail = NA,
                           R = 2000,
                           boot = "ifneeded",
                           always.analytical = FALSE ) {
  
  
  ##### Check for Bad Input #####
  if ( t2 < 0 ) stop("Heterogeneity cannot be negative")
  
  # the second condition is needed for Shiny app:
  #  if user deletes the input in box, then it's NA instead of NULL
  if ( ! is.na(se.M) ) {
    if (se.M < 0) stop("se.M cannot be negative")
  }
  
  if ( ! is.na(se.t2) ) {
    if (se.t2 < 0) stop("se.t2 cannot be negative")
  }
  
  ##### Messages When Not All Output Can Be Computed #####
  if ( is.na(se.M) | is.na(se.t2) ) message("Cannot compute inference without se.M and \nse.t2.\n Returning only point estimates.")
  
  
  ##### Point Estimates #####
  # same regardless of tail
  Z = (q - M) / sqrt(t2)
  
  if ( tail == "above" ) phat = 1 - pnorm(Z)
  else if ( tail == "below" ) phat = pnorm(Z)
  
  extreme = (phat < 0.15 | phat > 0.85) 
  
  if ( extreme & always.analytical == FALSE ) {
    if ( boot == "ifneeded" ) {
      warning("The estimated proportion is close to 0 or 1,\n so the theoretical CI may perform poorly. Using \nBCa bootstrapping instead.")
      
      boot.res = suppressWarnings( my_boot( data = dat, 
                                            parallel = "multicore",
                                            R = R, 
                                            statistic = get_stat,
                                            # last two arguments are being passed to get_stat
                                            q = q,
                                            tail = tail ) )
      
      bootCIs = boot.ci(boot.res, type="bca")
      lo = round( bootCIs$bca[4], 2 )
      hi = round( bootCIs$bca[5], 2 )
      SE = NA  
    } 
    
    if ( boot == "never" ) {
      warning("The estimated proportion is close to 0 or 1,\n so the theoretical CI may perform poorly. Should use \nBCa bootstrapping instead.")
      SE = lo = hi = NA
    }
  }
  
  
  if ( !extreme | always.analytical == TRUE ) {
    # do inference only if given needed SEs
    if ( !is.na(se.M) & !is.na(se.t2) ){
      
      ##### Delta Method Inference on Original Scale #####
      term1.1 = se.M^2 / t2
      term1.2 = ( se.t2^2 * ( q - M )^2 ) / ( 4 * t2^3 )
      term1 = sqrt( term1.1 + term1.2 )
      
      SE = term1 * dnorm(Z)
      
      # confidence interval
      tail.prob = ( 1 - CI.level ) / 2
      lo = max( 0, phat + qnorm( tail.prob )*SE )
      hi = min( 1, phat - qnorm( tail.prob )*SE )
    } else {
      SE = lo = hi = NA
    }
  }
  
  # return results
  res = data.frame( Est = phat, 
                    SE = SE,
                    lo = lo, 
                    hi = hi ) 
  rownames(res) = NULL
  res
}



################################### GET P-HAT STAT FOR USE WITH BOOT PACKAGE ################################### 

# defaults are for applied example #1
get_stat = function( original,
                     indices,
                     q, 
                     tail,
                     measure = "SMD", 
                     method = "REML",
                     yi.name = "yi", 
                     vi.name = "vi" ) {
  
  b = original[indices,]
  
  library(metafor)
  
  # keep track of whether there was an error fitting the model to these data
  # (Fisher convergence problems)
  got.error <<- FALSE
  
  tryCatch( {
    # meta-analyze the bootstrapped data
    mb = rma.uni( yi = b[[yi.name]],
                  vi = b[[vi.name]],
                  measure=measure,
                  knha = TRUE,
                  method = method,
                  # last argument helps prevent convergence problems
                  #control=list(stepadj=0.5, maxiter=1000) )
                  control=list(stepadj=0.5) )
    
    Mb = mb$b
    t2b = mb$tau2
    
    phat = suppressMessages( prop_stronger2( q = q,
                                             M = Mb,
                                             t2 = t2b,
                                             CI.level = CI.level,
                                             tail = tail,
                                             boot = "never" )$Est )  # this time don't bootstrap b/c we only need point est
    
  }, error = function(err) {
    message(err)
    # needs to be superassignment because inside fn
    errors <<- c( errors, err$message )
    got.error <<- TRUE
  } )
  
  if(got.error == FALSE) return(phat)
  
}



################################### MODIFIED FROM BOOT PACKAGE ################################### 

# see section called "MM additions"
# I minimally modified the function so that it can proceed even if some of the bootstrap iterates run into errors
# (in this case, Fisher convergence issues) because the boot package version gets confused about dimension mismatches

# source internal boot package functions
# source("bootfuns.R")

my_boot = function (data, statistic, R, sim = "ordinary", stype = c("i", 
                                                                    "f", "w"), strata = rep(1, n), L = NULL, m = 0, weights = NULL, 
                    ran.gen = function(d, p) d, mle = NULL, simple = FALSE, ..., 
                    parallel = c("no", "multicore", "snow"), ncpus = getOption("boot.ncpus", 
                                                                               1L), cl = NULL) 
{
  
  
  call <- match.call()
  stype <- match.arg(stype)
  if (missing(parallel)) 
    parallel <- getOption("boot.parallel", "no")
  parallel <- match.arg(parallel)
  have_mc <- have_snow <- FALSE
  if (parallel != "no" && ncpus > 1L) {
    if (parallel == "multicore") 
      have_mc <- .Platform$OS.type != "windows"
    else if (parallel == "snow") 
      have_snow <- TRUE
    if (!have_mc && !have_snow) 
      ncpus <- 1L
    loadNamespace("parallel")
  }
  if (simple && (sim != "ordinary" || stype != "i" || sum(m))) {
    warning("'simple=TRUE' is only valid for 'sim=\"ordinary\", stype=\"i\", n=0', so ignored")
    simple <- FALSE
  }
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
    runif(1)
  seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  n <- NROW(data)
  if ((n == 0) || is.null(n)) 
    stop("no data in call to 'boot'")
  temp.str <- strata
  strata <- tapply(seq_len(n), as.numeric(strata))
  t0 <- if (sim != "parametric") {
    if ((sim == "antithetic") && is.null(L)) 
      L <- empinf(data = data, statistic = statistic, stype = stype, 
                  strata = strata, ...)
    if (sim != "ordinary") 
      m <- 0
    else if (any(m < 0)) 
      stop("negative value of 'm' supplied")
    if ((length(m) != 1L) && (length(m) != length(table(strata)))) 
      stop("length of 'm' incompatible with 'strata'")
    if ((sim == "ordinary") || (sim == "balanced")) {
      if (isMatrix(weights) && (nrow(weights) != length(R))) 
        stop("dimensions of 'R' and 'weights' do not match")
    }
    else weights <- NULL
    if (!is.null(weights)) 
      weights <- t(apply(matrix(weights, n, length(R), 
                                byrow = TRUE), 2L, normalize, strata))
    if (!simple) 
      i <- index.array(n, R, sim, strata, m, L, weights)
    original <- if (stype == "f") 
      rep(1, n)
    else if (stype == "w") {
      ns <- tabulate(strata)[strata]
      1/ns
    }
    else seq_len(n)
    t0 <- if (sum(m) > 0L) 
      statistic(data, original, rep(1, sum(m)), ...)
    else statistic(data, original, ...)
    rm(original)
    t0
  }
  else statistic(data, ...)
  pred.i <- NULL
  fn <- if (sim == "parametric") {
    ran.gen
    data
    mle
    function(r) {
      dd <- ran.gen(data, mle)
      statistic(dd, ...)
    }
  }
  else {
    if (!simple && ncol(i) > n) {
      pred.i <- as.matrix(i[, (n + 1L):ncol(i)])
      i <- i[, seq_len(n)]
    }
    if (stype %in% c("f", "w")) {
      f <- freq.array(i)
      rm(i)
      if (stype == "w") 
        f <- f/ns
      if (sum(m) == 0L) 
        function(r) statistic(data, f[r, ], ...)
      else function(r) statistic(data, f[r, ], pred.i[r, 
                                                      ], ...)
    }
    else if (sum(m) > 0L) 
      function(r) statistic(data, i[r, ], pred.i[r, ], 
                            ...)
    else if (simple) 
      function(r) statistic(data, index.array(n, 1, sim, 
                                              strata, m, L, weights), ...)
    else function(r) statistic(data, i[r, ], ...)
  }
  RR <- sum(R)
  res <- if (ncpus > 1L && (have_mc || have_snow)) {
    if (have_mc) {
      parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus)
    }
    else if (have_snow) {
      list(...)
      if (is.null(cl)) {
        cl <- parallel::makePSOCKcluster(rep("localhost", 
                                             ncpus))
        if (RNGkind()[1L] == "L'Ecuyer-CMRG") 
          parallel::clusterSetRNGStream(cl)
        res <- parallel::parLapply(cl, seq_len(RR), fn)
        parallel::stopCluster(cl)
        res
      }
      else parallel::parLapply(cl, seq_len(RR), fn)
    }
  }
  else lapply(seq_len(RR), fn)
  #t.star <- matrix(, RR, length(t0))  # ~~~ MM commented out
  
  
  # ~~~~~ MM added
  # number of non-NULL elements of the results vector
  #browser()
  RR = length(unlist(res))
  nulls = sapply( res, is.null)
  res = res[ !nulls ]
  t.star <- matrix(, RR, length(t0))
  
  # without this, boot.CI gets confused about number of replicates
  R = RR
  browser()
  # ~~~~~ end of MM additions
  
  
  for (r in seq_len(RR)) t.star[r, ] <- res[[r]]
  if (is.null(weights)) 
    weights <- 1/tabulate(strata)[strata]
  boot.return(sim, t0, t.star, temp.str, R, data, statistic, 
              stype, call, seed, L, m, pred.i, weights, ran.gen, mle)
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                      FNs FOR MODEL 2 (REGULAR RANDOM EFFECTS)                       #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


################################ FN: RANDOM-EFFECTS CORRECTED ESTIMATE ################################ 

# minially modified from Hedges' package to also provide SEs with known selection probability
correct_est_re = function( effect,
                          v,
                          weights ) {
  
  res = weightfunct_mm( effect = effect, 
                        v = v,
                        weights = weights )
  
  
  H = res[[2]]$hessian
  
  # remove the parts corresponding to the fixed p parameter since its entries are 0
  H = H[1:2,1:2]
  
  SEs = sqrt(diag(solve(H)))
  
  bhat.adj = res[[2]]$par[2]
  t2.adj = res[[2]]$par[1]
  
  bhat.adj.lo = bhat.adj - qnorm(.975) * SEs[2]
  bhat.adj.hi = bhat.adj + qnorm(.975) * SEs[2]
  
  t2.adj.lo = t2.adj - qnorm(.975) * SEs[1]
  t2.adj.hi = t2.adj + qnorm(.975) * SEs[1]
  
  return( data.frame( 
    p.select = weights[2],  # only return the one weight that isn't 1
    param = c("b", "t2"),
    est = c(bhat.adj, t2.adj),
    se = c(SEs[2], SEs[1]),
    lo = c(bhat.adj.lo, t2.adj.lo),
    hi = c(bhat.adj.hi, t2.adj.hi)
  )
  )
}





########################### SLURM FUNCTIONS ###########################

# These just generate the sbatch files

# DO NOT CHANGE THE INDENTATION IN THE BELOW OR ELSE SLURM 
#  WILL SILENTLY IGNORE THE BATCH COMMANDS DUE TO EXTRA WHITESPACE!!
sbatch_skeleton <- function() {
  return(
"#!/bin/bash
#################
#set a job name  
#SBATCH --job-name=JOBNAME
#################  
#a file for job output, you can check job progress
#SBATCH --output=OUTFILE
#################
# a file for errors from the job
#SBATCH --error=ERRORFILE
#################
#time you think you need; default is one hour
#SBATCH --time=JOBTIME
#################
#quality of service; think of it as job priority
#SBATCH --qos=QUALITY
#################
#submit to both owners and normal partition
#SBATCH -p normal,owners
#################
#number of nodes you are requesting
#SBATCH --nodes=NODENUMBER
#################
#memory per node; default is 4000 MB
#SBATCH --mem=MEMPERNODE
#you could use --mem-per-cpu; they mean what we are calling cores
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=MAILTYPE
#################
#who to send email to; please change to your email
#SBATCH  --mail-user=USER_EMAIL
#################
#task to run per node; each node has 16 cores
#SBATCH --ntasks=TASKS_PER_NODE
#################
#SBATCH --cpus-per-task=CPUS_PER_TASK
#now run normal batch commands

ml load R
R -f PATH_TO_R_SCRIPT ARGS_TO_R_SCRIPT")
}



generateSbatch <- function(sbatch_params, runfile_path = NA, run_now = F) {
  
  #sbatch_params is a data frame with the following columns
  #jobname: string, specifies name associated with job in SLURM queue
  #outfile: string, specifies the name of the output file generated by job
  #errorfile: string, specifies the name of the error file generated by job
  #jobtime: string in hh:mm:ss format, max (maybe soft) is 48:00:00 
  #specifies the amoung of time job resources should be allocated
  #jobs still running after this amount of time will be aborted
  #quality: kind of like priority, normal works
  #node_number, integer: the number of nodes (computers w/16 cpus each) to allocate 
  #mem_per_node, integer: RAM, in MB, to allocate to each node
  #mailtype, string: ALL, BEGIN, END, FAIL: what types of events should you be notified about via email
  #user_email string: email address: email address to send notifications
  #tasks_per_node: integer, number of tasks, you should probably use 1
  #cpus_per_task: integer, 1-16, number of cpus to use, corresponds to number of available cores per task
  #path_to_r_script: path to r script on sherlock
  #args_to_r_script: arguments to pass to r script on command line
  #write_path: where to write the sbatch file
  #server_sbatch_path: where sbatch files will be stored on sherlock
  #runfile_path is a string containing a path at which to write an R script that can be used to run
  #the batch files generated by this function. 
  #if NA, no runfile will be written
  #run_now is a boolean specifying whether batch files should be run as they are generated
  
  sbatches <- list()
  if (!is.na(runfile_path)) {
    outfile_lines <- c(paste0("# Generated on ",  Sys.time()))
  }
  for (sbatch in 1:nrow(sbatch_params) ) {
    gen_batch <- sbatch_skeleton()
    #set job name
    if (is.null(sbatch_params$jobname[sbatch])) { 
      gen_batch <- gsub("JOBNAME", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("JOBNAME", sbatch_params$jobname[sbatch], gen_batch) 
    }
    #set outfile name
    if (is.null(sbatch_params$outfile[sbatch])) { 
      gen_batch <- gsub("OUTFILE", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("OUTFILE", sbatch_params$outfile[sbatch], gen_batch) 
    }
    #set errorfile name
    if (is.null(sbatch_params$errorfile[sbatch])) { 
      gen_batch <- gsub("ERRORFILE", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("ERRORFILE", sbatch_params$errorfile[sbatch], gen_batch) 
    }
    #set jobtime
    if (is.null(sbatch_params$jobtime[sbatch])) { 
      gen_batch <- gsub("JOBTIME", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("JOBTIME", sbatch_params$jobtime[sbatch], gen_batch) 
    }
    #set quality
    if (is.null(sbatch_params$quality[sbatch])) { 
      gen_batch <- gsub("QUALITY", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("QUALITY", sbatch_params$quality[sbatch], gen_batch) 
    }
    #set number of nodes
    if (is.null(sbatch_params$node_number[sbatch])) { 
      gen_batch <- gsub("NODENUMBER", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("NODENUMBER", sbatch_params$node_number[sbatch], gen_batch) 
    }
    #set memory per node
    if (is.null(sbatch_params$mem_per_node[sbatch])) { 
      gen_batch <- gsub("MEMPERNODE", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("MEMPERNODE", sbatch_params$mem_per_node[sbatch], gen_batch) 
    }
    #set requested mail message types
    if (is.null(sbatch_params$mailtype[sbatch])) { 
      gen_batch <- gsub("MAILTYPE", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("MAILTYPE", sbatch_params$mailtype[sbatch], gen_batch) 
    }
    #set email at which to receive messages
    if (is.null(sbatch_params$user_email[sbatch])) { 
      gen_batch <- gsub("USER_EMAIL", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("USER_EMAIL", sbatch_params$user_email[sbatch], gen_batch) 
    }
    #set tasks per node
    if (is.null(sbatch_params$tasks_per_node[sbatch])) { 
      gen_batch <- gsub("TASKS_PER_NODE", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("TASKS_PER_NODE", sbatch_params$tasks_per_node[sbatch], gen_batch) 
    }
    #set cpus per task
    if (is.null(sbatch_params$cpus_per_task[sbatch])) { 
      gen_batch <- gsub("CPUS_PER_TASK", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("CPUS_PER_TASK", sbatch_params$cpus_per_task[sbatch], gen_batch) 
    }
    #set path to r script
    if (is.null(sbatch_params$path_to_r_script[sbatch])) { 
      gen_batch <- gsub("PATH_TO_R_SCRIPT", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("PATH_TO_R_SCRIPT", sbatch_params$path_to_r_script[sbatch], gen_batch) 
    }
    #set args to r script
    if (is.null(sbatch_params$args_to_r_script[sbatch])) { 
      gen_batch <- gsub("ARGS_TO_R_SCRIPT", "unnamed", gen_batch) 
    } else { 
      gen_batch <- gsub("ARGS_TO_R_SCRIPT", sbatch_params$args_to_r_script[sbatch], gen_batch) 
    }
    
    #write batch file
    if (is.null(sbatch_params$write_path[sbatch])) { 
      cat(gen_batch, file = paste0("~/sbatch_generated_at_", gsub(" |:|-", "_", Sys.time()) ), append = F)
    } else { 
      cat(gen_batch, file = sbatch_params$write_path[sbatch], append = F)
    }
    
    if (!is.na(sbatch_params$server_sbatch_path[sbatch])) {
      outfile_lines <- c(outfile_lines, paste0("system(\"sbatch ", sbatch_params$server_sbatch_path[sbatch], "\")"))
    } 
    sbatches[[sbatch]] <- gen_batch
  }
  if (!is.na(runfile_path)) {
    cat(paste0(outfile_lines, collapse = "\n"), file = runfile_path)
  }
  if(run_now) { system(paste0("R -f ", runfile_path)) } 
  
  return(sbatches)
}




########################### FN: RETURN FILES THAT AREN'T COMPLETED ###########################

# Given a folder path for results and a common beginning of file name of results files
#   written by separate workers in parallel, stitch results files together into a
#   single csv.

# .max.sbatch.num: If not passed, defaults to largest number in actually run jobs.

sbatch_not_run = function(.results.singles.path,
                          .results.write.path,
                          .name.prefix,
                          .max.sbatch.num = NA ) {
  
  # get list of all files in folder
  all.files = list.files(.results.singles.path, full.names=TRUE)
  
  # we only want the ones whose name includes .name.prefix
  keepers = all.files[ grep( .name.prefix, all.files ) ]
  
  # extract job numbers
  sbatch.nums = as.numeric( unlist( lapply( strsplit( keepers, split = "_"), FUN = function(x) x[5] ) ) )
  
  # check for missed jobs before the max one
  if ( is.na(.max.sbatch.num) ) .max.sbatch.num = max(sbatch.nums)
  all.nums = 1 : .max.sbatch.num
  missed.nums = all.nums[ !all.nums %in% sbatch.nums ]
  
  # give info
  print( paste("The max job number is: ", max(sbatch.nums) ) )
  print( paste( "Number of jobs that weren't run: ",
                ifelse( length(missed.nums) > 0, length(missed.nums), "none" ) ) )
  
  if( length(missed.nums) > 0 ) {
    setwd(.results.write.path)
    write.csv(missed.nums, "missed_job_nums.csv")
  }
  
  return(missed.nums)
  
}
