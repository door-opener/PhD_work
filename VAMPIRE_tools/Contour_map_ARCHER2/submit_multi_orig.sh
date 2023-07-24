#!/bin/bash
# Slurm job options (job-name, compute nodes, job time)
#SBATCH --job-name=MultiSerialOnComputes
#SBATCH --time=48:00:00
#SBATCH --nodes=8
#SBATCH --tasks-per-node=128
#SBATCH --cpus-per-task=1

# Replace [budget code] below with your budget code (e.g. t01)
#SBATCH --account=e05-bulk-smw
#SBATCH --partition=standard
#SBATCH --qos=long

module load PrgEnv-cray/8.0.0
module load epcc-setup-env
module load load-epcc-module

module load cray-python
module load genmaskcpu

export OMP_NUM_THREADS=1

# Method to get list of directories into an array #
shopt -s nullglob
list=(NODE*)
shopt -u nullglob

# Getting list of node ids #
list2=( $(scontrol show hostnames $SLURM_JOB_NODELIST) )

for i in "${!list[@]}"; do

    cd ${list[$i]}

    # List of inner dirs #
    shopt -s nullglob
    inner_list=(STEP*)
    shopt -u nullglob

    for j in "${!inner_list[@]}"; do

        pass=$((j + 1))
        maskcpu=$( genmaskcpu 128 $pass 1 1 )

        cd ${inner_list[$j]}
        echo "genmask cpu 128 $pass 1 1"
        echo "nodeidx -> $i"
        echo "$pass   -> $maskcpu"
        echo "nodeid  -> ${list2[$i]}"
        #srun --cpu-bind=mask_cpu:${maskcpu} --nodelist=${list2[$i]} --nodes=1 --ntasks=1 --tasks-per-node=1 \
        #--oversubscribe --mem=1500M /work/e05/e05/isa/bin/vampire-serial > out_${i}_${j}.out &
        cd ../

    done
    cd ../
done

wait
