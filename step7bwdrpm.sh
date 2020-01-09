

#!/bin/bash

###1. sort bedgraph, 2.create bigwig files which can be used in visualization in UCSC genome browser: 


root_dir1=/scratch/ww346
root_dir=//scratch/ww346/analyze/Str_data
geno=hg19
study_name=Stau2_eCLIP

samples="Hep_Stau2_1_drpm Hep_Stau2_2_drpm"

PROG_DIR=$root_dir1/analyze/UCSC
cd $root_dir/$study_name/rawout/
#genofasta=/scratch/ww346/analyze/annotation/data/ucsc/genomes/$geno/$geno.fa
#fext='.clipped' # '/accepted_hits' for tophat ; '' #for novoalign; .sam.unique ; .BWASW 
#fext2=".pAreads" #'.pAreads' or ''


for sample in $samples; do
 #totalReadNum=`wc -l $sample$fext.sam$fext2 | sed s/[[:blank:]].*//`
 echo step7, ____ $sample, convert to bw ____
 
 bedtools sort -faidx $PROG_DIR/$geno.chrom.sizes -i $sample.p.bedGraph > $sample.sorted.p.bedGraph    
 bedtools sort -faidx $PROG_DIR/$geno.chrom.sizes -i $sample.m.bedGraph > $sample.sorted.m.bedGraph 
# rm $sample.p.bedGraph    ##this step finish in bash
# rm $sample.m.bedGraph    ##this step finish in bash
 $PROG_DIR/bedGraphToBigWig $sample.sorted.p.bedGraph $PROG_DIR/$geno.chrom.sizes $sample.sorted.p.bw 
 $PROG_DIR/bedGraphToBigWig $sample.sorted.m.bedGraph $PROG_DIR/$geno.chrom.sizes $sample.sorted.n.bw 

rm $sample.p.bedGraph
rm $sample.m.bedGraph
done



