#===============================================================================
#         FILE: frombam.eCLIP.pip.py  step4counts.py 
#
#        USAGE: python frombam.eCLIP.pip.py --project --<OPTIONS> 
#
#  DESCRIPTION: Fully automatic pipline processing eCLIP data starting from *bam files
# 				step4: count eCLIP #reads with first position alignment
#
#      OPTIONS:  --help;                                       show the help message and exit
#                --project(required) <project_name>; define the project name;
#                --rootdir(required) <path_to_your_root_dir>; define the rootdir;  
#                --sradir(required) <path_to_your_sra_dir>; define the sra dir;
#                --refdir(required) <path_to_your_ref_dir>; define the reference dir;#                
#                --threads(required) <# of threads>; define # of threads;
#                --genome(Optional) <mm9(default)/hg19/rn5>; define genome version;
#
# 
#       AUTHOR: Wei Wang  	wwei320@gmail.com
#      
#      VERSION: 2.0
#      CREATED: 2019-Dec-05
#===============================================================================
####import section
import matplotlib; matplotlib.use('agg')
import itertools
import numpy as np
import pandas as pd
from pylab import *
import os
import sys
import ntpath
import re
import glob
import scipy
import collections
import subprocess
import time
# import multiprocessing
import argparse
import unittest

def strandinfo_SINGLE(logfile,outfile):
	f=open( logfile ).readlines()
	F = 1
	R = 1
	for line in f:
		if "++,--" in line:
			inputline=line.split( ": " )[1]
			F=float(inputline.replace('\n',''))
		elif "+-,-+" in line:
			inputline=line.split( ": " )[1]
			R=float(inputline.replace('\n',''))
	with open(outfile, 'w') as f:	
		if F/R>=3:
			f.write('forward\n')
			print str(F/R) + ' strand is forward'
		elif R/F>=3:
			f.write('invert\n')
			print str(F/R) + ' strand is invert'
		else:
			f.write('NONE\n')
			print str(F/R) + ' strand is none'
		f.close()	
	os.system('rm '+logfile)

def strandinfo_PAIR(logfile):
	f=open( logfile ).readlines()
	STRINFO = ''
	F = 1
	R = 1
	for line in f:
		if "1++,1--,2+-,2-+" in line:
			inputline=line.split( ": " )[1]
			F=float(inputline.replace('\n',''))
		elif "1+-,1-+,2++,2--" in line:
			inputline=line.split( ": " )[1]
			R=float(inputline.replace('\n',''))	
		if F/R>=3:
			STRINFO='FR' 
		elif R/F>=3:
			STRINFO='RF'		
		else:
			STRINFO='P_NONE'
	return STRINFO
	
def handel_PAIR(strKEY,fastqfile1,fastqfile2,bamfile):	
	if strKEY=='RF':		
		fastqfile=fastqfile1		
	else:
		fastqfile=fastqfile2
	return fastqfile

def downloadcheck(fqmdir,sample,seqtype,tail1,tail2):
	if seqtype=='SINGLE':
		exists = os.path.isfile(fqmdir+sample+'.fastq')
		if exists:
			print sample+"Downloading finished"
			DNCHECK='YES'
		else:
			print sample+"Downloading error, try re-download once"
			DNCHECK='NO'

	elif seqtype=='PAIRED':
		exists1 = os.path.isfile(fqmdir+sample+tail1+'.fastq')
		exists2 = os.path.isfile(fqmdir+sample+tail2+'.fastq')
		if exists1 and exists2:
			print sample+"Downloading finished"
			DNCHECK='YES'
		else:
			print sample+"Downloading error, try to re-download once"
			DNCHECK='NO'
	return 	DNCHECK
		

##arg setting parts###
parser = argparse.ArgumentParser(description="Use this to run RNAseq-APA automatic pipeline.")
parser.add_argument("--project", action="store", dest='project',default='', metavar='<name_of_project>', help="define the project name")
parser.add_argument("--rootdir", action="store", dest='rootdir',default='', metavar='<rootdir>', help="define the rootdir")
parser.add_argument("--sradir", action="store", dest='sradir',default='', metavar='<sradir>', help="define the sradir")
parser.add_argument("--refdir", action="store", dest='refdir',default='', metavar='<refdir>', help="define the refdir")
parser.add_argument("--genodir", action="store", dest='genodir',default='', metavar='<genodir>', help="define the genome dir for mapping")
parser.add_argument("--genome", action="store", dest='genome',default='mm9', metavar='<mm9(default)/hg19/rn4/...>', help="define the genome version")
parser.add_argument("--threads", action="store", dest='threads',default='', metavar='<threads>', help="define the number of threads")

args=parser.parse_args()

#####python setting #####
rootdir=args.rootdir
project=args.project
sradir=args.sradir
refdir=args.refdir
geno=args.genome
genoDir=args.genodir
CPUS=args.threads

CPUS=str(CPUS)

####automatic setting section####
scrdir=os.path.dirname(os.path.abspath(__file__))
mapper='_star'
tail1='_1'
tail2='_2'
UTRdir=refdir
fqmdir=os.path.join(rootdir, project+'/rawfastq/')
# genoDir=os.path.join(rootdir, 'data/ucsc/genomes/'+geno+mapper)
##toggled by ww: samoutdir=os.path.join(rootdir, project+'/rawsam/')
samoutdir=os.path.join(rootdir, project+'/rawsam')
rawoutdir=os.path.join(rootdir, project+'/rawout/')
testdir=os.path.join(rootdir, project+'/test/')
plotdir=os.path.join(rootdir, project+'/plot/')
reportdir=os.path.join(rootdir, project+'/report/')
samplefile=rootdir+project+'/sample_list.txt'
PROG_DIR=os.path.dirname("/scratch/ww346/analyze/UCSC/")

if not os.path.exists(fqmdir):
    os.makedirs(fqmdir)	
if not os.path.exists(samoutdir):
    os.makedirs(samoutdir)	
if not os.path.exists(rawoutdir):
    os.makedirs(rawoutdir)	
if not os.path.exists(testdir):
    os.makedirs(testdir)	
if not os.path.exists(plotdir):
    os.makedirs(plotdir)	
if not os.path.exists(reportdir):
    os.makedirs(reportdir)	
	
dfsample = pd.read_table(samplefile)
samples=dfsample.Run.unique().tolist()


# os.chdir(dndir)

for sample in samples:

	
 	seqtype=dfsample[dfsample.Run==sample].LibraryLayout.values[0]
 	samplename=dfsample[dfsample.Run==sample].samplename.values[0]
	

### step1 Mapping ####	
	samfile=samoutdir+'/'+sample+'.sam'
	samlog=samoutdir+'/'+sample+'.sam.log'
	bamfile=samoutdir+'/'+sample+'.Aligned.sortedByCoord.out.bam'


#### step2 checking strand ####	
	###samtools view -b -h -q 255 SRR1057941.Aligned.sortedByCoord.out.bam > SRR1057941.Aligned.sortedByCoord.out.unique.bam ##this step finish in bash
	bamfile=samoutdir+'/'+sample+'.Aligned.sortedByCoord.out.unique.bam'

	print "\n--------------Step2: Determing strand info("+seqtype+"),Sample:"+sample+"----------------"	
	logfile=testdir+'check.txt'
	outfile=testdir+'strand.txt'

# ####step3 Count reads using R####
	newbedfile=samoutdir+'/'+samplename+'.bed'
	bedfile=samoutdir+'/'+sample+'.bed'
	os.system('mv '+bedfile +' '+ newbedfile)
	
	print "\n--------------Step3: Counting reads("+seqtype+"),Sample:"+sample+"----------------"		
	prjdir=os.path.join(rootdir, project+'/')
	arg1=prjdir
	arg2=geno
	arg3=samoutdir+'/'
	arg4=UTRdir
	arg5=CPUS
	arg6=samplename
	###Rscriptfile0=scrdir+'/RNAseq.APAprofile.R'
	Rscriptfile0=scrdir+'/step4counts.R'
	print "\n--------------step4 counts----------------"
	####  RUN Rscript ###
	subprocess.call(['Rscript', Rscriptfile0, arg1, arg2, arg3,arg4,arg5, arg6])	
# 	os.system('rm '+newbamfile)
# 	if seqtype=='SINGLE':
# 		os.system('rm '+fqmdir+sample+'.*')
# 	elif seqtype=='PAIRED':
# 		os.system('rm '+fqmdir+sample+'_*')
	



	#######make rpm-based bedgraph 
	cmd31='bedtools sort -faidx '+PROG_DIR+'/'+geno+'.chrom.sizes -i '+samoutdir+'/'+samplename+'_rpm.p.bedGraph > '+samoutdir+'/'+samplename+'_rpm.sorted.p.bedGraph'     ##this step finish in bash
	print cmd31
	os.system(cmd31)

	cmd32='bedtools sort -faidx '+PROG_DIR+'/'+geno+'.chrom.sizes -i '+samoutdir+'/'+samplename+'_rpm.m.bedGraph > '+samoutdir+'/'+samplename+'_rpm.sorted.m.bedGraph'     ##this step finish in bash
	print cmd32
	os.system(cmd32)

	cmd41='rm '+samoutdir+'/'+samplename+'_rpm.p.bedGraph'    ##this step finish in bash
	print cmd41
	os.system(cmd41)

	cmd42='rm '+samoutdir+'/'+samplename+'_rpm.m.bedGraph'    ##this step finish in bash
	print cmd42
	os.system(cmd42)
		





	#######make count-based bedgraph 
	cmd31='bedtools sort -faidx '+PROG_DIR+'/'+geno+'.chrom.sizes -i '+samoutdir+'/'+samplename+'_count.p.bedGraph > '+samoutdir+'/'+samplename+'_count.sorted.p.bedGraph'     ##this step finish in bash
	print cmd31
	os.system(cmd31)

	cmd32='bedtools sort -faidx '+PROG_DIR+'/'+geno+'.chrom.sizes -i '+samoutdir+'/'+samplename+'_count.m.bedGraph > '+samoutdir+'/'+samplename+'_count.sorted.m.bedGraph'     ##this step finish in bash
	print cmd32
	os.system(cmd32)

	cmd41='rm '+samoutdir+'/'+samplename+'_count.p.bedGraph'    ##this step finish in bash
	print cmd41
	os.system(cmd41)

	cmd42='rm '+samoutdir+'/'+samplename+'_count.m.bedGraph'    ##this step finish in bash
	print cmd42
	os.system(cmd42)





