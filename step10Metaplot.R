

# ########use this bed and bw files to make matrix for metagene plot, see Metamatrix.sh



######make a bed file with TES in the middle +/-500 or +/-200
# require(dplyr)
# setwd('/scratch/ww346/analyze/APAanalyzer/REF')

# refseq <- read.table(paste0("hg19.RefSeq.bed"), header = FALSE, sep = "\t",
#               quote = "", stringsAsFactor = FALSE)

# refseq2 <- refseq[,c(1:4,6:8)]

# names(refseq2) <- c('chromosome','TSS','TES','initiator','strand','cds_start','cds_end')

# bed1 <- refseq2%>%  
#   filter(strand=='+')%>%
# select(chromosome,TSS, TES, strand)

# bed2 <- refseq2%>% 	
# 	filter(strand=='-')%>%
# select(chromosome,TSS, TES, strand)

# setwd('/scratch/ww346/analyze/Str_data/TIA1_eCLIP/rawout')
# write.table(bed1, 'TSS_TES.p.bed', quote=F, sep="\t", row.names=F, col.names=F)
# write.table(bed2, 'TSS_TES.n.bed', quote=F, sep="\t", row.names=F, col.names=F)

########use this bed and bw files to make matrix for metagene plot, see Metamatrix.sh




###########after making matrix, (did normalization before making matrix) and plot metaplot 

#####1. normalization is to subtract input rpm from eclip rpm: result is similar to iclip result

#####2. normalization is to calculate log2(eclip_rpm / input_rpm): result has large noise
require(dplyr)
require(ggplot2)
setwd('/scratch/ww346/analyze/Str_data/TIA1_eCLIP/report')
ma1=read.table('Hep_TIA1_1_drpm.sorted.TES.p.tab',sep='\t',header=F, quote="", stringsAsFactor=FALSE,skip=3, row.names=NULL)

ma2=read.table('Hep_TIA1_1_drpm.sorted.TES.n.tab',sep='\t',header=F, quote="", stringsAsFactor=FALSE,skip=3, row.names=NULL)

ma=rbind(ma1,ma2)

ma[is.na(ma)] <- 0

mam=ma%>%
	select(contains("V"))%>%
	mutate(rsum=rowSums(.))%>%
	mutate_at(vars(contains("V")),all_vars(./rsum))

total = as.data.frame(apply(mam, 2, sum))
names(total)='csum'
total$cols= seq(1,101)



p1 <- ggplot() + 
  geom_line(aes(y = csum, x = cols), color = 'black',
                           data = total[1:100,], stat="identity")+
  #scale_y_continuous(limits = c(0, 1))+
  ggtitle(paste0("Hep_TIA1_r1_eCLIP_drpm_TES_metagene_plot \n#Genes ",nrow(ma))) + labs(x="-500nt TES +500nt", y="Sum of gene_wise normalized rpm")+
  # geom_hline(yintercept=0, color='black')+
  # geom_vline(xintercept=0, color='black')+
  # theme(axis.title.x = element_blank(),
  #       panel.background = element_blank())
  theme_classic()

setwd('/scratch/ww346/analyze/Str_data/TIA1_eCLIP/plot')
jpeg(paste0('Hep_TIA1_r1_eCLIP_drpm_TES_metagene.tiff'), res=300, width=2000, height=2000, type="cairo")

  p1
dev.off()




######make matrix column-wise csum for input and eclip data seperately and plot log2 fc 
######this method does not give result similar to iclip data

# require(dplyr)
# require(ggplot2)
# setwd('/scratch/ww346/analyze/Str_data/TIA1_eCLIP/report')

# for (sample in c('TIA1_1','INPUT')){

# ma1=read.table(paste0('Hep_',sample,'_rpm_sorted.TES.p.tab'),sep='\t',header=F, quote="", stringsAsFactor=FALSE,skip=3, row.names=NULL)

# ma2=read.table(paste0('Hep_',sample,'_rpm_sorted.TES.n.tab'),sep='\t',header=F, quote="", stringsAsFactor=FALSE,skip=3, row.names=NULL)

# ma=rbind(ma1,ma2)

# ma[is.na(ma)] <- 0

# mam=ma%>%
#   select(contains("V"))%>%
#   mutate(rsum=rowSums(.))%>%
#   mutate_at(vars(contains("V")),all_vars(./rsum))

# total = as.data.frame(apply(mam, 2, sum))
# names(total)='csum'
# total$cols= seq(1,101)

# assign(sample,total)
# }

# mergedf=merge(TIA1_1,INPUT,by='cols')
# mergedf$log2rpm=log2(mergedf$csum.x/mergedf$csum.y)

# p1 <- ggplot() + 
#   geom_line(aes(y = log2rpm, x = cols), color = 'black',
#                            data = mergedf[1:100,], stat="identity")+
#   #scale_y_continuous(limits = c(0, 1))+
#   ggtitle(paste0("Hep_TIA1_r1_eCLIP_log2rpm_TES_metagene_plot \n#Genes ",nrow(ma))) + labs(x="-500nt TES +500nt", y="Sum of gene_wise normalized rpm")+
#   # geom_hline(yintercept=0, color='black')+
#   # geom_vline(xintercept=0, color='black')+
#   # theme(axis.title.x = element_blank(),
#   #       panel.background = element_blank())
#   theme_classic()

# setwd('/scratch/ww346/analyze/Str_data/TIA1_eCLIP/plot')
# jpeg(paste0('Hep_TIA1_r1_eCLIP_log2rpm_TES_metagene.tiff'), res=300, width=2000, height=2000, type="cairo")

#   p1
# dev.off()











