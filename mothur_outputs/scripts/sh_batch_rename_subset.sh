#!/bin/bash
#SBATCH --job-name=mothur_miseq_subsample_and_prep
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=180g
#SBATCH --time=36:00:00
#SBATCH --account=vdenef0
#SBATCH --partition=standard
#SBATCH --mail-user=bakerdyl@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL

#  Show list of CPUs you ran on, if you're running under PBS
echo $SLURM_JOB_NODELIST

#  Change to the directory you submitted from
if [ -n "$SLURM_SUBMIT_DIR" ]; then cd $SLURM_SUBMIT_DIR; fi
pwd

source ~/.bashrc
module load Bioinformatics mothur python3.9-anaconda fasttree
conda init
#renames files to replace first underscore ("_") with the characters "sep"
#purpose: keeps d3/d31 classification in group name when fastq is read into mothur
for filename in *.fastq; do mv "./$filename" "./$(echo "$filename" | sed -e 's/_S.*_L001_/_/1')";  done
for filename in *.fastq; do mv "./$filename" "./$(echo "$filename" | sed -e 's/_001//1')";  done

#subsets data, keeping 1000 randomly selected reads
#purpose: speeds analysis throughout the rest of the pipeline, but does represent some data loss
conda activate /nfs/turbo/lsa-dudelabs/conda_envs/miniconda/envs/seqtk/
for fq in *R1.fastq; do seqtk sample -s100 $fq 10000 > $fq\.sub.fastq; done
for fq in *R2.fastq; do seqtk sample -s100 $fq 10000 > $fq\.sub.fastq; done

#Move subsetted sequences to another folder
find ./ -name  "*sub*" -exec mv {} ../miseq_subsample \;
cd ../miseq_subsample

#Remove extra underscores in sample names that interfere with "make.file" command in mothur
find . -name 'J*' -exec bash -c 'echo mv $0 ${0/_/}' {} \;
find . -name 'J*' -exec bash -c 'echo mv $0 ${0/_/}' {} \;
#update naming conventions for easier down-the-line analysis
find . -name '[[:digit:]]*_[[:digit:]]*.fastq' -exec bash -c 'mv $0 ${0/_/point}' {} \;
find . -name '[[:digit:]]*.fastq' -exec bash -c 'echo mv $0 ${0/DF/D31}' {} \;

#run analysis of subsampled sequences in mothur once everything you need (silva databases, sanger sequences, miseq sequences, etc.) are in the folder "miseq_subsample"
mothur mothur_miseq.batch.taxass

#Make tree of isolates, some renaming will need to be done between the last step of mothur pipeline and tree making, but I will leave that up to you.
fasttree -gtr -nt -spr 4 -mlacc 2 -slownni isolates_asv.fasta > isolates.asv.ruben.nwk