#PBS -l walltime=5:00:00
#PBS -l select=1:ncpus=1:mem=5gb
#PBS -N dictionaries

# change directory to repo folder
cd $PBS_O_WORKDIR 
while [[ $PWD != '/' && ${PWD##*/} != 'translational-data-science-2024-2025-Group12' ]]
do 
    cd ..
done

# print job info 
echo "Running $PBS_JOBNAME from $PWD"
echo "Current branch: $(git branch --show-current)"
echo "Latest commit: $(git rev-parse --short HEAD)"
(git diff-index --quiet HEAD -- && echo "No untracked changes") || echo "UNTRACKED CHANGES"

# load path constants
source config/hpc_paths.sh 

# activate virtual environment
module load anaconda3/personal
source activate tdsgrp12

# run script
cd scripts
rscript_filepath="01-preprocess_proteins/01-dictionaries.R" # specify script to run

if [ ! -f $rscript_filepath ]; then # check if file exists
    echo "File not found in $PWD - $rscript_filepath"
else 
    # create output directory if it does not exist
    rscript_output_dir="$OUTPUT_DIR/$(dirname $rscript_filepath)"
    mkdir -p $rscript_output_dir

    # run script 
    Rscript $rscript_filepath $DATA_DIR $rscript_output_dir
fi

source utils/mod_folder_permissions.sh

# log completion time
echo -e "\n=========="
echo "Completed: $(date +%F_%T)"
secs=$SECONDS
printf "Time taken: %02dh:%02dm:%02ds\n" $((secs/3600)) $((secs%3600/60)) $((secs%60))
