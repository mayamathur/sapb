

####################### CLEAN UP ####################### 

# clean up the directory
rm /home/groups/manishad/SAPB/results/*
  rm /home/groups/manishad/SAPB/results/short/*
  rm /home/groups/manishad/SAPB/results/long/*
  
  rm /home/groups/manishad/SAPB/results/overall_stitched/*
  rm /home/groups/manishad/SAPB/sbatch_files/rm*
  rm /home/groups/manishad/SAPB/sbatch_files/slurm*
  
  rm -r /home/groups/manishad/SAPB/sbatch_files/*

  
########################### SET SIMULATION PARAMETERS MATRIX ###########################

# FOR CLUSTER USE
path = "/home/groups/manishad/SAPB"
setwd(path)

k = c(20, 40, 80, 300)  # number of clusters prior to selection
per.cluster = c(5)  # studies per cluster
mu = c(-99)  # RE distribution mean (will be ignored and overwritten since using exponential)
V = c(2, 1, 0)  # RE heterogeneity
V.gam = c(0, 0.5)  # variance of random intercepts (can't be > V because that's total heterogeneity!)
sei.min = 1  # runif lower bound for study SEs
sei.max = c(1.5)  # runif upper bound for study SEs
eta = c(100, 50, 20, 10, 1)  # selection prob
q = c(1) # no longer relevant because not using Phat
#q = c(2.2, 1.2)
boot.reps = c(0)
bt.meta.model = c("rma.uni")
bt.type = c("wtd.vanilla")
orig.meta.model = rev( c("robumeta.lazy") )
true.dist = "exp"
SE.corr = FALSE


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
                          SE.corr )

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
                       "SE.corr" )


# remove scenarios with V = 0 but V.gam > 0
scen.params = scen.params[ !(scen.params$V == 0 & scen.params$V.gam > 0), ]

# remove scenarios with nonsense combinations of methods and parameters
stupid = rep( FALSE, nrow(scen.params) )
stupid[ scen.params$V > 0 & scen.params$orig.meta.model == "fixed" ] = TRUE
stupid[ scen.params$V == 0 & scen.params$orig.meta.model != "fixed" ] = TRUE
scen.params = scen.params[ stupid == FALSE, ]

# check it
table(scen.params$V, scen.params$V.gam, scen.params$orig.meta.model)

# if true.dist is exponential, then mu is predetermined by the variance
# this will not be used but is meant to record what the truth is
scen.params$mu[ scen.params$true.dist == "exp" ] = scen.params$V[ scen.params$true.dist == "exp" ] - scen.params$V.gam[ scen.params$true.dist == "exp" ]


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


# write the csv file of params (to Sherlock)
write.csv( scen.params, "scen_params.csv", row.names = FALSE )


########################### GENERATE SBATCHES ###########################

# load functions for generating sbatch files
path = "/home/groups/manishad/SAPB"
setwd(path)
source("helper_sim_study_SAPB.R")

# number of sbatches to generate (i.e., iterations within each scenario)
n.reps.per.scen = 1000
n.reps.in.doParallel = 200
n.files.per.scen = n.reps.per.scen / n.reps.in.doParallel # this is also the expected number of rows per long results file
( n.files = n.files.per.scen * n.scen )

scen.name = rep( scen.params$scen.name, each = ( n.files / n.scen ) )
jobname = paste("job", 1:n.files, sep="_")
outfile = paste("rm_", 1:n.files, ".out", sep="")
errorfile = paste("rm_", 1:n.files, ".err", sep="")
write_path = paste(path, "/sbatch_files/", 1:n.files, ".sbatch", sep="")
runfile_path = paste(path, "/testRunFile.R", sep="")

sbatch_params <- data.frame(jobname,
                            outfile,
                            errorfile,
                            jobtime = "2:00:00",  # had used 1:00:00 for main results
                            quality = "normal",
                            node_number = 1,
                            mem_per_node = 64000,
                            mailtype =  "NONE",
                            user_email = "mmathur@stanford.edu",
                            tasks_per_node = 16,
                            cpus_per_task = 1,
                            path_to_r_script = paste(path, "/doParallel_SAPB.R", sep=""),
                            args_to_r_script = paste("--args", jobname, scen.name, sep=" "),
                            write_path,
                            stringsAsFactors = F,
                            server_sbatch_path = NA)

generateSbatch(sbatch_params, runfile_path)

n.files

# 800 total
path = "/home/groups/manishad/SAPB"
setwd( paste(path, "/sbatch_files", sep="") )
for (i in 3:800) {
  system( paste("sbatch -p owners /home/groups/manishad/SAPB/sbatch_files/", i, ".sbatch", sep="") )
  
  # max hourly submissions seems to be 300, which is 12 seconds/job
  if ( i %% 300 == 0 ) Sys.sleep(60)  # delay in seconds
}

# run just one
# sbatch -p owners /home/groups/manishad/SAPB/sbatch_files/1.sbatch


# ######## If Running Only Some Jobs To Fill Gaps ########
# 
# run in Sherlock ml load R
path = "/home/groups/manishad/S
APB"
setwd(path)
source("helper_sim_study_SAPB.R")

missed.nums = sbatch_not_run( "/home/groups/manishad/SAPB/results/long",
                              "/home/groups/manishad/SAPB/results",
                              .name.prefix = "long_results",
                              .max.sbatch.num = 560 )



setwd( paste(path, "/sbatch_files", sep="") )
for (i in missed.nums) {
  system( paste("sbatch -p owners /home/groups/manishad/SAPB/sbatch_files/", i, ".sbatch", sep="") )
}