
# overleaf.dir = "~/Dropbox/Apps/Overleaf/SAPB manu and appendix/R_objects"
# results.dir = "~/Dropbox/Personal computer/Independent studies/Sensitivity analysis for publication bias (SAPB)/Linked to OSF (SAPB)/Simulation study/Results"
# data.dir = "~/Dropbox/Personal computer/Independent studies/Sensitivity analysis for publication bias (SAPB)/Linked to OSF (SAPB)/Simulation study/Results/*2019-1-6 in manuscript"




overleaf.dir = "~/Dropbox/Apps/Overleaf/SAPB manu and appendix/R_objects"
data.dir = "~/Dropbox/Personal computer/Independent studies/Sensitivity analysis for publication bias (SAPB)/Linked to OSF (SAPB)/Simulation study/Results/*2020-3-21 Expo with selection on SE"
results.dir = data.dir

################################### HELPER FN ###################################

plot_group = function( 
  .title = NA,
  .ylab = NA,
  .legend = TRUE,
  .y.name,
  .refline.y = NA,  # fixed at a specific value
  .limits = NA, 
  .breaks = NA,
  .V = NA,
  .V.gam = NA
  ) {  # value of group var
  
  if ( !is.na(.V) & !is.na(.V.gam) ) {
    temp = agg[ agg$V == .V & agg$V.gam == .V.gam, ] 
  }

  # with different reference lines in each facet
  if ( .y.name == "MuEst" ) {
    p =   ggplot( temp, aes_string( x="eta", y=.y.name, color="orig.pretty" ) ) +
      geom_hline(aes(yintercept=mu), lty = 2, color="red") +
      geom_line(lwd=1) +
      geom_point(size=2) +
      theme_bw() +
      scale_color_manual(values=colors) +
      
      
      facet_grid( k.pretty ~ mu,
                  labeller = label_bquote( cols = mu ~ "=" ~ .(mu) ) ) +
      
      guides(color=guide_legend(title="Model"))
  }

  if ( .y.name != "MuEst" ) {
    p = 
      ggplot( temp, aes_string( x="eta", y=.y.name, color="orig.pretty" ) ) +
      geom_line(lwd=1) +
      geom_point(size=2) +
      theme_bw() +
      scale_color_manual(values=colors) +
      
      facet_grid( k.pretty ~ mu,
                  labeller = label_bquote( cols = mu ~ "=" ~ .(mu) ) ) +
      
      guides(color=guide_legend(title="Model")) +
      
      scale_x_continuous( limits = c(min(temp$eta),
                                     max(temp$eta)), 
                          breaks = unique(temp$eta) ) +
      xlab( bquote(eta) ) +
      ggtitle( bquote( tau^2 ~ "=" ~ .(.V) ~ ", Var" ~ zeta ~ "=" ~ .(.V.gam) ) )
  }
  
  # log scale 
  if ( .y.name == "MuCIWidth") p = p + scale_y_continuous(trans='log10')

  # y-label  
  if ( .y.name == "MuCover" ) p = p + ylab( bquote(hat(mu) ~ " coverage") )
  if ( .y.name == "MuCIWidth" ) p = p + ylab( "Median CI width" ) 
  if ( .y.name == "MuEst" ) p = p + ylab( bquote( "Median " ~ hat(mu) ) ) 
  
    # if we have a fixed reference line
    if ( !is.na(.refline.y) ) {
      p = p + geom_hline(aes(yintercept=.refline.y), lty = 2, color="red")
    }
  
  if ( !is.na(.limits[1]) & !is.na(.breaks[1]) ) {
    p = p + scale_y_continuous( limits=.limits, breaks=.breaks )
  }
  
  if ( !is.na(.ylab) ) {
    p = p + ylab(.ylab)
  }
  
  if ( !is.na(.title) ) {
    p = p + ggtitle(.title)
  }
  
  if ( .legend == TRUE ) {
    return(p)
  } else {
    return(p + theme(legend.position="none"))
  }
  
}  


################################### READ IN DATA ################################### 


# set up a dataframe for use with Overleaf
sim.res = list()

# read in stitched results (1 row per simulation iterate)
library(tidyverse)
setwd(data.dir)
d = read_csv( "stitched.csv" )

# remove idiosyncratic failed rows
good.rows = grep("scen", d$scen.name)
all.rows = 1:nrow(d)
bad.rows = all.rows[ !all.rows %in% good.rows ]
length(bad.rows)
#View( d[bad.rows,] )
if( length(bad.rows) > 0 ) d = d[ -bad.rows, ]

# also remove bootstrap rows since not in use
d = d[ d$Method == "Naive", ]

# should always be "exp"
table(d$true.dist)

# should be half TRUE and half FALSE
table(d$select.SE)
d %>% group_by(select.SE) %>%
  summarise( mean(SE.mean.emp, na.rm = TRUE) )


# other NA rows
( cant.be.na = names(d)[3:16] )
options(scipen=999)
apply( d[,cant.be.na], 2,
       function(x)prop.table(table(is.na(x))))

# remove very small proportion of missing data
d = d[ complete.cases( d[ , cant.be.na ] ), ]

d = d[ !is.na(d$scen.name), ]

d = droplevels(d)

# how many sim reps per scenario?
table(d$scen.name)

# min successful reps per scenario
summary( as.numeric(table(d$scen.name) ) )  # 864
sim.res["min.sim.reps"] = min(as.numeric(table(d$scen.name) ))



################################### INFO ABOUT N STUDIES AND RUNTIMES ################################### 

##### Number of Published Studies #####
# number of observed studies should decline with larger eta
( n.studies = d %>% group_by(k, V, V.gam, eta, mu, select.SE) %>% summarise( k.obs = mean(k.obs),
                                         n.nonsig = median(n.nonsig),
                                         SE.corr = median(SE.corr.emp) ) )
setwd(results.dir)
write.csv(n.studies, "n_studies.csv", row.names = FALSE)

##### Max Runtime by Scenario #####
( max.runtimes = d %>% group_by(scen.name) %>%
    summarise( max.doParallel.min = max(doParallel.min) ) )
write.csv(max.runtimes, "max_runtimes.csv", row.names = FALSE)

##### Missing Data by Scenario #####
# I think this should only happen for wtd.score when it fails to converge
# since we already removed idiosyncratic missing data
( prop.missing = d %>% group_by(scen.name) %>%
    summarise( prop.missing = mean( is.na(MuCover) ) ) )
write.csv(prop.missing, "prop_missing.csv", row.names = FALSE)


################################### AGGREGATE AND RECODE CERTAIN BRATTY VARIABLES ################################### 


# quick check
View( d %>% group_by(select.SE) %>%
        summarise( mean(MuEst, na.rm = TRUE),
                   median(MuEst, na.rm = TRUE),
                   mean(MuCover, na.rm = TRUE) )
)

# bm
( analysis.vars = c(names(d)[ 19:38 ], "MuCIWidth") )

# these are the scenario variables; taken from genSbatch
grouping.vars = c("k",
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
                  "select.SE" )


# cast analysis variables as numeric
want.numeric = which( !names(d) %in% c("scen.name",
                                        "orig.meta.model",
                                         "bt.meta.model",
                                         "Method",
                                          "MuCover",
                                       "T2Cover",
                                       "PCover") )
d[,want.numeric] = sapply( d[,want.numeric], as.numeric )

# can't include these in the is.numeric above since they become 1 or 2
d$MuCover = as.logical(d$MuCover)
d$T2Cover = as.logical(d$T2Cover)
d$PCover = as.logical(d$PCover)

# aggregate statistical metrics of interest by scenario
# use median for some variables and mean for others
#take.median = c("MuCIWidth", "MuEst", "n.nonsig", "k.obs")
take.median = c("MuCIWidth", "n.nonsig", "k.obs")
( agg = d %>% group_by_at( vars( one_of( grouping.vars ) ) ) %>%
    summarise_at( vars( one_of(analysis.vars[ ! analysis.vars %in% take.median ]) ),
                  mean, na.rm = TRUE ) )
dim(agg) # should equal the number of unique scenarios (unless sims aren't done running)

# add variables that need median
agg2 = d %>% group_by_at( vars( one_of( grouping.vars ) ) ) %>%
  summarise_at( vars( take.median ),
                median, na.rm = TRUE )
agg$MuCIWidth = agg2$MuCIWidth
#agg$MuEst = agg2$MuEst
agg$n.nonsig = agg2$n.nonsig


# for plotting joy
agg$k.pretty = as.factor( paste( agg$per.cluster * agg$k,
                                 " latent studies", sep="") )
levels(agg$k.pretty)
# for ggplot ordering
# ~~~ IF YOU'VE CHANGED THE SCENARIOS, MAKE SURE NUMBERS HERE ARE STILL RIGHT
agg$k.pretty = factor( agg$k.pretty, levels = c("100 latent studies",
                                                "200 latent studies",
                                                "400 latent studies",
                                                "1000 latent studies") )
levels(agg$k.pretty)

agg$mu.pretty = paste( "mu = ", agg$mu, sep="")


agg$orig.pretty = rep(NA, nrow(agg))
agg$orig.pretty[ agg$orig.meta.model == "fixed" ] = "Common-effects"
agg$orig.pretty[ agg$orig.meta.model == "wtd.score" ] = "Wtd. score"
agg$orig.pretty[ agg$orig.meta.model == "robumeta" ] = "Robust (score)"
agg$orig.pretty[ agg$orig.meta.model == "robumeta.lazy" &
                   agg$V.gam == 0 ] = "Robust independent"
agg$orig.pretty[ agg$orig.meta.model == "robumeta.lazy" &
                   agg$V.gam > 0 ] = "Robust clustered"


agg$select.SE.pretty = rep(NA, nrow(agg))
agg$select.SE.pretty[ agg$select.SE == TRUE ] = "Selected for small SE"
agg$select.SE.pretty[ agg$select.SE == FALSE ] = "Not selected on SE"

agg = droplevels(agg)


# #### ~~~ SANITY CHECK
# 
# fake = d[ d$k == 20 & d$per.cluster == 5 & d$mu == 0.8 & d$eta == 5 & d$V == 1 & d$V.gam == 0,]

#View( agg[ agg$k == 20 & agg$per.cluster == 5 & agg$mu == 0.8 & agg$eta == 10 & agg$V == 1 & agg$V.gam == 0,] )

################################### PLOTS FOR MAIN TEXT ################################### 

# keep only recommended combinations
agg2 = agg[ agg$orig.meta.model %in% c("fixed", "robumeta.lazy"), ]
table(agg2$orig.meta.model, agg2$V.gam)
agg2 = droplevels(agg2)

colors = c("orange", "red", "black", "blue")
shapes = c(16,2)

##### Coverage #####
( coverage.table = agg2 %>% group_by(orig.meta.model, V.gam, select.SE) %>%
  summarise( MeanCover = mean(MuCover),
             MinCover = min(MuCover) ) )

( coverage.table = agg2 %>% group_by(orig.meta.model, V.gam) %>%
    summarise( MeanCover = mean(MuCover),
               MinCover = min(MuCover) ) )




##### Point Estimate #####
ggplot( agg2, aes_string( x="eta",
                          y="MuEst",
                          color="orig.pretty",
                          shape = "select.SE.pretty") ) +
  
  geom_hline(aes(yintercept=mu), lty = 2, color="red") +
  
  geom_line(lwd=1) +
  geom_point(size=4) +
  theme_bw() +
  scale_color_manual(values=colors) +
  scale_shape_manual( values = shapes)+
  
  #facet_grid( k.pretty ~ mu.pretty ) +
  
  facet_grid( k.pretty ~ mu,
              labeller = label_bquote( cols = mu ~ "=" ~ .(mu) ) ) +
  
  guides( color=guide_legend(title="Model"),
          shape = guide_legend(title="Additional selection on SE")) +
  
  scale_y_continuous( limits = c(0, 1),
                      breaks = seq(0, 1, .2)) +
  
  scale_x_continuous( limits = c(min(agg2$eta),
                                 max(agg2$eta)), 
                      breaks = unique(agg2$eta) ) +
  xlab( bquote(eta) ) +
  ylab( bquote( "Mean " ~ hat(mu) ) ) 

# bm
  
setwd(overleaf.dir)
ggsave( "muhat_exponential.pdf",
        height = 10,
        width = 8,
        units = "in")
setwd(results.dir)
ggsave( "muhat_exponential.pdf",
        height = 10,
        width = 8,
        units = "in")


##### CI Width #####
ggplot( agg2, aes_string( x="eta",
                          y="MuCIWidth",
                          color="orig.pretty",
                          shape = "select.SE.pretty") ) +
  
  geom_line(lwd=1) +
  geom_point(size=3) +
  theme_bw() +
  scale_color_manual(values=colors) +
  scale_shape_manual( values = shapes)+
  
  #facet_grid( k.pretty ~ mu.pretty ) +
  
  facet_grid( k.pretty ~ mu,
              labeller = label_bquote( cols = mu ~ "=" ~ .(mu) ) ) +
  
  guides(color=guide_legend(title="Model")) +
  
  
  scale_x_continuous( limits = c(min(agg2$eta),
                                 max(agg2$eta)), 
                      breaks = unique(agg2$eta) ) +
  scale_y_continuous( trans = "log10" ) + 
  xlab( bquote(eta) ) +
  ylab( "Median CI width" ) 

setwd(overleaf.dir)
ggsave( "width_exponential.pdf",
        height = 10,
        width = 8,
        units = "in")
setwd(results.dir)
ggsave( "width_exponential.pdf",
        height = 10,
        width = 8,
        units = "in")


# among scenarios with at least 10 nonsigs
summary( agg$MuCIWidth[ agg$n.nonsig >= 10 ] )
summary( agg$MuCIWidth[ agg$n.nonsig < 10 ] )


sim.res["median.width.big"] = summary( agg$MuCIWidth[ agg$n.nonsig >= 10 ] )["Median"]
sim.res["lo.width.big"] = summary( agg$MuCIWidth[ agg$n.nonsig >= 10 ] )["1st Qu."]
sim.res["hi.width.big"] = summary( agg$MuCIWidth[ agg$n.nonsig >= 10 ] )["3rd Qu."]

sim.res["median.width.small"] = summary( agg$MuCIWidth[ agg$n.nonsig < 10 ] )["Median"]
sim.res["lo.width.small"] = summary( agg$MuCIWidth[ agg$n.nonsig < 10 ] )["1st Qu."]
sim.res["hi.width.small"] = summary( agg$MuCIWidth[ agg$n.nonsig < 10 ] )["3rd Qu."]


##### Number of Nonsig Studies #####

ggplot( agg2, aes_string( x="eta", y="n.nonsig", color="orig.pretty",
                          shape = "select.SE.pretty") ) +
  
  geom_line(lwd=1) +
  geom_point(size=4) +
  theme_bw() +
  scale_color_manual(values=colors) +
  scale_shape_manual( values = shapes)+
  
  facet_grid( k.pretty ~ mu,
              labeller = label_bquote( cols = mu ~ "=" ~ .(mu) ) ) +
  
  guides(color=guide_legend(title="Model")) +
  
  scale_y_continuous( trans="log10") +
  
  scale_x_continuous( limits = c(min(agg2$eta),
                                 max(agg2$eta)), 
                      breaks = unique(agg2$eta) ) +
  xlab( bquote(eta) ) +
  ylab( "Median no. published nonaffirmative studies" ) 

setwd(overleaf.dir)
ggsave( "nnonsig_exponential.pdf",
        height = 10,
        width = 8,
        units = "in")
setwd(results.dir)
ggsave( "nnonsig_exponential.pdf",
        height = 10,
        width = 8,
        units = "in")


################################### PLOTS FOR APPENDIX ################################### 

library(ggplot2)

agg = agg[ !is.na(agg$mu), ]

##### Coverage #####
cover.limits = c(0.5, 1)
cover.breaks = seq(0.5, 1, .1)
# plot_group(.y.name = "MuCover",
#            .V = 0,
#            .V.gam = 0,
#            .refline.y = 0.95,
#            .limits = cover.limits, 
#            .breaks = cover.breaks )

plot_group(.y.name = "MuCover",
           .V = 1,
           .V.gam = 0.5,
           .refline.y = 0.95 )

setwd(overleaf.dir)
ggsave( "cover_cluster_appendix.pdf",
        height = 10,
        width = 8,
        units = "in")
setwd(results.dir)
ggsave( "cover_cluster_appendix.pdf",
        height = 10,
        width = 8,
        units = "in")



plot_group(.y.name = "MuCover",
           .V = 1,
           .V.gam = 0,
           .refline.y = 0.95 )

setwd(overleaf.dir)
ggsave( "cover_nocluster_appendix.pdf",
        height = 10,
        width = 8,
        units = "in")
setwd(results.dir)
ggsave( "cover_nocluster_appendix.pdf",
        height = 10,
        width = 8,
        units = "in")


##### Width #####
# note that the y-axis here is on log-10 scale, though the numbers are 
#  on regular scale
# plot_group(.y.name = "MuCIWidth",
#            .V = 0,
#            .V.gam = 0 )

plot_group(.y.name = "MuCIWidth",
           .V = 1,
           .V.gam = 0.5 )

setwd(overleaf.dir)
ggsave( "width_cluster_appendix.pdf",
        height = 10,
        width = 8,
        units = "in")
setwd(results.dir)
ggsave( "width_cluster_appendix.pdf",
        height = 10,
        width = 8,
        units = "in")
     
plot_group(.y.name = "MuCIWidth",
           .V = 1,
           .V.gam = 0 )

setwd(overleaf.dir)
ggsave( "width_nocluster_appendix.pdf",
        height = 10,
        width = 8,
        units = "in")
setwd(results.dir)
ggsave( "width_nocluster_appendix.pdf",
        height = 10,
        width = 8,
        units = "in")


##### Point Estimate #####
# plot_group(.y.name = "MuEst",
#            .V = 0,
#            .V.gam = 0)
# 
# plot_group(.y.name = "MuEst",
#            .V = 1,
#            .V.gam = 0)
# 
# plot_group(.y.name = "MuEst",
#            .V = 1,
#            .V.gam = 0.5)




#################### SAVE RESULTS ####################

here("Results")
write.csv( as.data.frame(sim.res), "simres.csv" )

write.csv( coverage.table, "coveragetable.csv" )

write.csv( as.data.frame(sim.res), "agg.csv" )


