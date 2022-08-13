#!/bin/bash

#SBATCH --partition=bgmp
#SBATCH --job-name=demux
#SBATCH --cpus-per-task=1
#SBATCH --account=bgmp

conda activate bgmp_py310

R1="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
I1="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
I2="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
R2="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
indexes="/projects/bgmp/shared/2017_sequencing/indexes.txt"
demultiplexed="/projects/bgmp/pgenge/bioinfo/Bi622/Demultiplex/Assignment-the-third/demultiplexed"

/usr/bin/time -v \
/projects/bgmp/pgenge/bioinfo/Bi622/Demultiplex/Assignment-the-third/demux.py \
$R1 $I1 $I2 $R2 $indexes -o $demultiplexed