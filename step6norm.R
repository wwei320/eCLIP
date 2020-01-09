

#####to normalize eclip rpm with input rpm by subtraction formula
require('GenomicRanges')
require('rtracklayer')
require('dplyr')
require('rlist')

setwd('/scratch/ww346/analyze/Str_data/Stau2_eCLIP/rawout')

for (sample in c('Stau2_1','Stau2_2')){


	gr1p <- import(paste0('Hep_',sample,'_rpm.sorted.p.bw'),format="BigWig")

	gr1n <- import(paste0('Hep_',sample,'_rpm.sorted.n.bw'),format="BigWig")

	gr2p <- import('Hep_INPUT_rpm.sorted.p.bw',format="BigWig")

	gr2n <- import('Hep_INPUT_rpm.sorted.n.bw',format="BigWig")

	strand(gr1p)='+'
	strand(gr1n)='-'
	strand(gr2p)='+'
	strand(gr2n)='-'

	gr1= c(gr1p,gr1n)

	gr2= c(gr2p,gr2n)

	# gr2$rpminp=gr2$score
	# gr2$score=NULL

	#######keep input sample reads if eclip sample #reads>0
	hits <- findOverlaps(gr2, gr1,type='equal')
	  idx <- unique(subjectHits(hits))
	  idxq<-unique(queryHits(hits))
	  #ranges$rpmvivo=grvivo$rpm[idx]
	 noidxq=list.remove(seq(1,length(gr2)),idxq) 

	gr1$score[idx] =gr1$score[idx]-gr2$score[idxq]
	# gr2$score[noidxq] = 0-gr2$score[noidxq]
	# gr3=c(gr1, gr2[noidxq]) 

	  grdf=as.data.frame(gr1)

	grdf1<-grdf%>%filter(strand=='+')%>%select(seqnames,start,end,score)%>%mutate(end=end+1)

	options(scipen=999)
	write.table(grdf1,paste0('Hep_',sample,'_drpm.p.bedGraph'),sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)

	grdf2<-grdf%>%filter(strand=='-')%>%select(seqnames,start,end,score)%>%mutate(start=start-1)
	write.table(grdf2,paste0('Hep_',sample,'_drpm.m.bedGraph'),sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
	 
}
