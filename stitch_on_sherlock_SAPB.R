

# run this interactively in ml load R

path = "/home/groups/manishad/SAPB"
setwd(path)
source("helper_sim_study_SAPB.R")

######## FOR STITCHING LONG FILES ########

# find out how many files have run
setwd("/home/groups/manishad/SAPB/results/long")
length(list.files())

s = stitch_files(.results.singles.path = "/home/groups/manishad/SAPB/results/long",
                 .results.stitched.write.path = "/home/groups/manishad/SAPB/results/overall_stitched",
                 .name.prefix = "long",
                 .stitch.file.name="stitched.csv")

# # for stitching long files on my Desktop
# s = stitch_files(.results.singles.path = "~/Desktop",
#                  .results.stitched.write.path = "~/Desktop",
#                  .name.prefix = "long",
#                  .stitch.file.name="stitched.csv")

table(s$scen)

table(s$orig.meta.model)

######## SNEAK PEEK AT RESULTS ########

library(crayon, lib.loc = "/home/groups/manishad/Rpackages/") # for dplyr
library(dplyr, lib.loc = "/home/groups/manishad/Rpackages/")
library(foreach, lib.loc = "/home/groups/manishad/Rpackages/")
library(doParallel, lib.loc = "/home/groups/manishad/Rpackages/")
library(boot, lib.loc = "/home/groups/manishad/Rpackages/")
library(metafor, lib.loc = "/home/groups/manishad/Rpackages/")
library(data.table, lib.loc = "/home/groups/manishad/Rpackages/")
library(robumeta, lib.loc = "/home/groups/manishad/Rpackages/")
library(fansi, lib.loc = "/home/groups/manishad/Rpackages/")  # for dplyr
library(utf8, lib.loc = "/home/groups/manishad/Rpackages/")  # for dplyr
library(cli, lib.loc = "/home/groups/manishad/Rpackages/")  # for dplyr
library(dplyr)

# my.vars = c("MuEst", "MuCover", "T2Est", "T2Cover", "PEst", "PCover")
# s %>% group_by(sei.max, q, k, V.gam) %>%
#   filter( Method == "Naive" ) %>%
#   summarise_at( vars(my.vars), mean )

# WHY ARE THERE WEIRD ROWS??
valid.names = paste("scen", 1:9999, sep=".")
bad = which( !s$scen.name %in% valid.names )
length(bad)
#View( d[bad,] )
if( length(bad) > 0 ) s = s[ -bad, ]
s = droplevels(s)

table(s$scen.name)

# sanity check -- number of observed studies should decline with larger eta
s %>% group_by(k, eta) %>% summarise( k.obs = mean(k.obs),
                                         n.nonsig = mean(n.nonsig),
                                      SE.corr = mean(SE.corr.emp) )



my.vars = c("MuEst", "MuCover", "MuLo", "MuHi",
            "T2Est", "T2Cover", "T2Lo", "T2Hi",
            "PEst", "PCover", "PLo", "PHi",
            "k.obs", "n.nonsig", "SE.corr.emp")

# mean performance among non-NA reps

s$MuCover = as.numeric(s$MuCover)

( agg = s %>% group_by(scen.name, Method, V, k ) %>%
    summarise_at( vars(my.vars), mean, na.rm = TRUE ) )

round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  
  df[,nums] <- round(df[,nums], digits = digits)
  
  (df)
}

# round for easy viewing
agg2 = round_df(agg, digits = 2)
as.data.frame(agg2)


# move it to Desktop
scp mmathur@login.sherlock.stanford.edu:/home/groups/manishad/SAPB/results/overall_stitched/stitched.csv ~/Desktop
Vegemite2017
scp mmathur@login.sherlock.stanford.edu:/home/groups/manishad/SAPB/scen_params.csv ~/Desktop
Vegemite2017

