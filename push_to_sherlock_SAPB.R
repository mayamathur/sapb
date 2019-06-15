

####################### CHECK IN ####################### 
# see the sbatches
cd /home/groups/manishad/SAPB/sbatch_files

sbatch -p manishad /home/groups/manishad/SAPB/sbatch_files/1.sbatch

# check on my running or pending jobs
squeue -u mmathur -t RUNNING
squeue -u mmathur -t PENDING


# see the datasets
cd /home/groups/manishad/SAPB/results/long
ls -l . | egrep -c '^-'

# see the errors
nano /home/groups/manishad/SAPB/sbatch_files/slurm*
  nano /home/groups/manishad/SAPB/sbatch_files/rm_1.err

# see the scen parameters
nano /home/groups/manishad/SAPB/scen_params.csv

# see the stitched results
nano /home/groups/manishad/SAPB/results/overall_stitched/sti*
  
  
  
####################### CODE -> SHERLOCK ####################### 

# push all the individual files
scp /Users/mmathur/Dropbox/Personal\ computer/Independent\ studies/Sensitivity\ analysis\ for\ publication\ bias\ \(SAPB\)/Simulation\ study/Code/* mmathur@login.sherlock.stanford.edu:/home/groups/manishad/SAPB

# Sherlock -> Desktop
scp mmathur@login.sherlock.stanford.edu:/home/groups/manishad/SAPB/results/overall_stitched/stitched.csv ~/Desktop


scp mmathur@login.sherlock.stanford.edu:/home/groups/manishad/SAPB/* ~/Desktop


####################### SHERLOCK -> DESKTOP (DEBUGGING) ####################### 

# move error file to Desktop
scp -r mmathur@sherlock:/home/groups/manishad/SAPB/sbatch_files/rm_19.err ~/Desktop
scp -r mmathur@sherlock:/home/groups/manishad/SAPB/sbatch_files/rm_19.out ~/Desktop

# move one sbatch file to Desktop
scp -r mmathur@login.sherlock.stanford.edu:/home/groups/manishad/SAPB/sbatch_files/2296.sbatch ~/Desktop


####################### RUN SBATCH ####################### 

# run one of them on Manisha's nodes
sbatch -p manishad /home/groups/manishad/SAPB/sbatch_files/1.sbatch
# not on Manisha's nodes
sbatch -p normal,owners /home/groups/manishad/SAPB/sbatch_files/1.sbatch




####################### RESULTS -> DESKTOP FOR ANALYSIS ####################### 

scp mmathur@login.sherlock.stanford.edu /home/groups/manishad/SAPB/results/overall_stitched/stitched.csv ~/Desktop

####################### CLEAN UP ####################### 

# clean up the directory
rm /home/groups/manishad/SAPB/results/*
  rm /home/groups/manishad/SAPB/results/short/*
  rm /home/groups/manishad/SAPB/results/long/*
  
  rm /home/groups/manishad/SAPB/results/overall_stitched/*
  rm /home/groups/manishad/SAPB/sbatch_files/rm*
  rm /home/groups/manishad/SAPB/sbatch_files/slurm*
  
  rm -r /home/groups/manishad/SAPB/sbatch_files/*
  