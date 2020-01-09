#!/bash


#SBATCH -p main         # Partition (job queue)
#SBATCH --job-name=Stau2      # Assign an 8-character name to your job, no spaces, no special characters
#SBATCH --nodes=1                # Number of compute nodes
#SBATCH --ntasks=1              # Number of tasks to run (often = cores) on each node
#SBATCH --mem=128000             # Total real memory required (MB) for each node
#SBATCH--cpus-per-task=24
#SBATCH --time=12:00:00
#SBATCH --output=slurm.%N.%j.out
#SBATCH --error=slurm.%N.%j.out
#SBATCH --export=ALL

##########for ENCODE eCLIP, skip step1 and step2, start with bam file 

scrdir=/scratch/ww346/analyze/Str_data/Stau2_eCLIP/script

### python step1download.py --rootdir /scratch/ww346/analyze/Str_data/  \

### python step2mapping.py --rootdir /scratch/ww346/analyze/Str_data/  \


########1. check strand and 2. get bed from bam files

# python step3strand.py --rootdir /scratch/ww346/analyze/Str_data/  \

				
########1. do first position alignment in bed file, 2. save as bedgraphand, and 3. sort bedgraph file
# python step4counts.py --rootdir /scratch/ww346/analyze/Str_data/  \
				# 	--project Stau2_eCLIP  \
				# --sradir /scratch/ww346/analyze/APAanalyzer/sra/sra/  \
				# --genodir /scratch/ww346/analyze/hg19_star/  \
				# --refdir /scratch/ww346/analyze/APAanalyzer/REF  \
				# --genome hg19  \
				# --threads 24	


#####make bigwig files from bedgraph files
# srun ./step5bw.sh

######normalization of eclip data by input data
# Rscript $scrdir/step6norm.R	


######sort bg and make bw for delta rpm samples
# srun ./step7bwdrpm.sh


######use deeptools to make matrix for metagene plots
# srun ./Metamatrix.sh

#--genodir /scratch/ww346/analyze/annotation/bowtie2/hg19_transcriptome_ww/  \     #use this path if map to transcriptome


# #######to count distribution of reads/or clusters in autr vs. cutr vs. cds
Rscript $scrdir/distribution.R 


######to clustering for iclip sam file (uniquely mapped sam from bam using ./bw.sh)
# conda activate py3
# python cluster.py
# conda deactivate


	



				