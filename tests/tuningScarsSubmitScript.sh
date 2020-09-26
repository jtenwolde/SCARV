#!/bin/bash
#!
#! Example SLURM job script for Darwin (Sandy Bridge, ConnectX3)
#! Last updated: Sat Apr 18 13:05:53 BST 2015
#!

#!#############################################################
#!#### Modify the options in this section as appropriate ######
#!#############################################################


#! sbatch directives begin here ###############################
#! Name of the job:
#SBATCH -J darwinjob
#! Which project should be charged:
#SBATCH -A MRC-BSU-SL2-CPU
#! How many whole nodes should be allocated?
#SBATCH --nodes=1
#! How many (MPI) tasks will there be in total? (<= nodes*16)
##SBATCH --ntasks=32
#! How much wallclock time will be required?
#SBATCH --time=06:00:00
#! How much memory in MB is required _per node_? Not setting this
#! will lead to a default of (1/16)*total memory per task.
#! Setting a larger amount per task increases the number of cores.
#SBATCH --mem=200000
#! What types of email messages do you wish to receive?
#SBATCH --mail-type=FAIL
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
##SBATCH --no-requeue

#! Do not change:
#SBATCH -p skylake-himem

#! sbatch directives end here (put any additional directives above this line)

#! Notes:
#! Charging is determined by core number*walltime.
#! The --ntasks value refers to the number of tasks to be launched by SLURM only. This
#! usually equates to the number of MPI tasks launched. Reduce this from nodes*16 if
#! demanded by memory requirements, or if OMP_NUM_THREADS>1.
#! Each task is allocated 1 core by default, and each core is allocated 3994MB. If this
#! is insufficient, also specify --cpus-per-task and/or --mem (the latter specifies
#! MB per node).

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')
#! ############################################################
#! Modify the settings below to specify the application's environment, location
#! and launch method:

#! Optionally modify the environment seen by the application
#! (note that SLURM reproduces the environment at submission irrespective of ~/.bashrc):
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
                   # REQUIRED - loads the basic environment

module load xz/5.2.2
module load pcre/8.38

source activate tf-gpu

export TMPDIR=/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/tmp
ipython /home/jwt44/git/constraint/packaging_scars/tests/tuningScars.py 125 very_common singleton /home/jwt44/pathoPercentiles_125_verycommon_rare.txt
