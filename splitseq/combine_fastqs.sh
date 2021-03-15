#!/bin/bash
#SBATCH --job-name=fastq
#SBATCH -p standard
#SBATCH -A vswarup_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --error=slurm-%J.err
#SBATCH --mem 8
#SBATCH --array=0-15

let index="$SLURM_ARRAY_TASK_ID"
let index="0"

fastqs="/dfs3b/swaruplab/smorabit/data/FIRE_mouse_2021/expdata/"
outdir="/dfs3b/swaruplab/smorabit/data/FIRE_mouse_2021/expdata_combined/"

sublibraries=($(ls $fastqs | cut -d '-' -f 5-6 | sort | uniq))
lib=${sublibraries[$index]}
echo $lib
zcat $fastqs*$lib* | gzip > $outdir$lib.fastq.gz
