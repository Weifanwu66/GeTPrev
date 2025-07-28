#!/bin/bash
#$ -N name
#$ -q queue
#$ -l h_rt=runtime
#$ -l mem_free=RAM
#$ -pe smp hpctasks
#$ -cwd
#$ -P account
##$ -S /bin/bash #smrt.q only?
export PATH="$HOME/.conda/envs/getprev/bin:$PATH" #Pass conda environment to the scheduler
vars
command
