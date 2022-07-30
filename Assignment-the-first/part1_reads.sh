#!/bin/bash
#SBATCH --partition=bgmp
#SBATCH --job-name=reads
#SBATCH --cpus-per-task=16
#SBATCH --account=bgmp 
#SBATCH --time=0-24:00:00

#Files
R1=/gpfs/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz
R2=/gpfs/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz
R3=/gpfs/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
R4=/gpfs/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz


conda activate bgmp_py310

/usr/bin/time \
  ./part1.py -f $R1 -o read1_quality.png -l 101 -t "Read 1"

  /usr/bin/time \
  ./part1.py -f $R4 -o read2_quality.png -l 101 -t "Read 2"

