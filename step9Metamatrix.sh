#!/bin/bash

########use bed file from Metaplot.R and bw files to make matrix for metagene plot

root_dir=/scratch/ww346/analyze/Str_data
geno=hg19
study_name=TIA1_eCLIP
rawout_dir=$root_dir/$study_name/rawout
refpoint=TES
out_dir=$root_dir/$study_name/report

cd $rawout_dir

###########matrix for plus and minus strand seperately
# samples="Hep_TIA1_1_drpm_sorted Hep_TIA1_2_rpm_sorted Hep_INPUT_rpm_sorted"
samples="Hep_TIA1_1_drpm.sorted"

for sample in $samples; do
	computeMatrix reference-point \
       --referencePoint $refpoint \
       -b 500 -a 500 \
       -R TSS_TES.p.bed \
       -S $sample.p.bw  \
       --skipZeros \
       -o $out_dir/$sample.TES.p.gz \
       --outFileNameMatrix $out_dir/$sample.TES.p.tab \
       --outFileSortedRegions $out_dir/$sample.TES.p.bed

	computeMatrix reference-point \
       --referencePoint $refpoint \
       -b 500 -a 500 \
       -R TSS_TES.n.bed \
       -S $sample.n.bw  \
       --skipZeros \
       -o $out_dir/$sample.TES.n.gz \
       --outFileNameMatrix $out_dir/$sample.TES.n.tab \
       --outFileSortedRegions $out_dir/$sample.TES.n.bed

done
# plotHeatmap \
#  -m mBCELL_TSS.p.gz\
#  -out mBCELL_TSS.p.png \
#  --heatmapHeight 10  \
#  --refPointLabel TES \
#  --regionsLabel RPM \
#  --plotTitle 'TIA1_iCLIP'
