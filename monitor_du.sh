START_TIME=$(date +%s); echo "starting getprev disk usage monitoring"

if [[ $(squeue -u $(whoami) | grep slurm2 | wc -l) -eq 1 ]]; then echo \
"single getprev instance detected on a slurm system. logging disk usage"; else echo \
"monitor_du.sh will not run. make sure you have a single getprev instance running on a slurm system."
fi

while [[ $(squeue -u $(whoami) | grep slurm2 | wc -l) -eq 1 ]]
do du database | tail -n1 >> disklog.txt; sleep 30; done

if [[ -f "disklog.txt" ]]; then echo "a disklog.txt file was detected"; fi

END_TIME=$(date +%s); echo "elapsed time for monitoring disk usage of getprev $(($END_TIME - $START_TIME)) seconds."