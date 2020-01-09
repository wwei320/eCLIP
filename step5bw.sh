######this is a .sh script to convert bedgraph files to bw files 


#!/bin/bash

###create bigwig files which can be used in visualization in UCSC genome browser: 
root_dir1=/scratch/ww346
root_dir=//scratch/ww346/analyze/Str_data
geno=hg19
study_name=Stau2_eCLIP
# samples="HEK293_rpm_sorted"
samples="Hep_Stau2_1_count Hep_Stau2_2_count Hep_INPUT_count Hep_Stau2_1_rpm Hep_Stau2_2_rpm Hep_INPUT_rpm"

PROG_DIR=$root_dir1/analyze/UCSC
cd $root_dir/$study_name/rawsam/


for sample in $samples; do
 #totalReadNum=`wc -l $sample$fext.sam$fext2 | sed s/[[:blank:]].*//`
 echo step8, ____ $sample, convert to bw ____
 $PROG_DIR/bedGraphToBigWig $sample.sorted.p.bedGraph $PROG_DIR/$geno.chrom.sizes $sample.sorted.p.bw 
 $PROG_DIR/bedGraphToBigWig $sample.sorted.m.bedGraph $PROG_DIR/$geno.chrom.sizes $sample.sorted.n.bw 

done

 mv *.bedGraph ../rawout
 mv *.bw ../rawout