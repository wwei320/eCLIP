########to count distribution of clip reads mapping to diffrent genomic regions
library(tidyr)
library(ggplot2)
library(gridExtra)
require(VennDiagram)
library(dplyr)
require(GenomicRanges)
library(AnnotationDbi)


####################functions for non-coding rna (optional), 3utr, 5utr, intron, intergenic (others), cds, cutr, autr ##########################################
create_noncoding_from_pAs <- function(pA.df,geno){
  if(geno == "mm9"){
  ncTx = read.csv("/scratch/ww346/analyze/Structure/dinghai/scriptfromDH/annotation/mm9_ncTx.csv", as.is = T)
  }else if (geno == "hg19"){
  ncTx = read.csv("/scratch/ww346/analyze/Structure/dinghai/scriptfromDH/annotation/hg19_ncTx.csv", as.is = T)
  }
  
  ncTx = GRanges(seqnames = ncTx$chrom, strand = ncTx$strand,
               ranges = IRanges(start = ncTx$txStart, end = ncTx$txEnd),
               gene_symbol = ncTx$name) #ncTx$name?

# remove duplicates
  ncTx = unique(ncTx)
  names(ncTx)=ncTx$gene_symbol
  ncTx               
}


create_5utr_from_pAs <- function(pA.df,geno){
  txdb = loadDb(paste0("/scratch/ww346/analyze/annotation/ucsc.txdb/",geno,".refGene.txdb.sqlite")) 
  if (geno=='mm9'){
   require(org.Mm.eg.db)
  IDDB <- org.Mm.eg.db 
  } else if (geno=='hg19'){
    require(org.Hs.eg.db)
    IDDB <- org.Hs.eg.db 
  }
  
  fiveUTRs = fiveUTRsByTranscript(txdb, use.names=T) 
  x = unlist(fiveUTRs)
  names(x) = mapIds(IDDB, keys = names(x), keytype =  "ACCNUM", column =  "SYMBOL")

  fiveUTRs = split(x, names(x))
  fiveUTRs = reduce(fiveUTRs)


  gr=unlist(fiveUTRs)                    
  gr$gene_symbol = names(gr)
  gr
}

create_3utr_from_pAs <- function(pA.df,geno){
  txdb = loadDb(paste0("/scratch/ww346/analyze/annotation/ucsc.txdb/",geno,".refGene.txdb.sqlite")) 
  if (geno=='mm9'){
   require(org.Mm.eg.db)
  IDDB <- org.Mm.eg.db 
  } else if (geno=='hg19'){
    require(org.Hs.eg.db)
    IDDB <- org.Hs.eg.db 
  }
  
  threeUTRs = threeUTRsByTranscript(txdb, use.names=T) 
  x = unlist(threeUTRs)
  names(x) = mapIds(IDDB, keys = names(x), keytype =  "ACCNUM", column =  "SYMBOL")

  threeUTRs = split(x, names(x))
  threeUTRs = reduce(threeUTRs)


  gr=unlist(threeUTRs)                    
  gr$gene_symbol = names(gr)
  gr
}


create_intron_from_pAs <- function(pA.df,geno){
  txdb = loadDb(paste0("/scratch/ww346/analyze/annotation/ucsc.txdb/",geno,".refGene.txdb.sqlite")) 
  if (geno=='mm9'){
   require(org.Mm.eg.db)
  IDDB <- org.Mm.eg.db 
  } else if (geno=='hg19'){
    require(org.Hs.eg.db)
    IDDB <- org.Hs.eg.db 
  }
  
  introns = intronsByTranscript(txdb, use.names=T)
  x = unlist(introns)
  names(x) = mapIds(IDDB, keys = names(x), keytype =  "ACCNUM", column =  "SYMBOL")

  introns = split(x, names(x))
  introns = reduce(introns)


  gr=unlist(introns)                    
  gr$gene_symbol = names(gr)
  gr
}


create_cds_from_pAs <- function(pA.df,geno){
  txdb = loadDb(paste0("/scratch/ww346/analyze/annotation/ucsc.txdb/",geno,".refGene.txdb.sqlite")) 
  if (geno=='mm9'){
   require(org.Mm.eg.db)
  IDDB <- org.Mm.eg.db 
  } else if (geno=='hg19'){
    require(org.Hs.eg.db)
    IDDB <- org.Hs.eg.db 
  }
  
  #exongr = exons(txdb, columns=c('EXONID','GENEID'), filter=list(gene_id=mergedPASS$gene_id))
  #rexongr = reduce(split(exongr,elementMetadata(exongr)$GENEID)) 
  #exonsby = exonsBy(txdb, by=c("tx", "gene"), columns=listColumns(txdb,"exon"), use.names=FALSE)
  CDSbygene = cdsBy(txdb, by=("gene"))
  x = unlist(CDSbygene)
  acc2sym = AnnotationDbi::select(IDDB, keys = names(x), keytype =  "ENTREZID", columns =  "SYMBOL")
  acc2sym[is.na(acc2sym$SYMBOL),]$SYMBOL = acc2sym[is.na(acc2sym$SYMBOL), ]$ENTREZID
  names(x) = acc2sym$SYMBOL
  CDSbygene = split(x, names(x))
  CDSbygene = reduce(CDSbygene)


  gr=unlist(CDSbygene)
  gr$gene_symbol = names(gr)
  gr
}

create_cUTR_from_pAs = function(pA.df,geno){
  refdir="/scratch/ww346/analyze/APAanalyzer/REF/"
  df= read.table(paste0(refdir,"SRS.",geno,".conserved.FL.txt"),header=TRUE)
  df$strand = df$Strand
  df$chr = df$Chrom
  df$start = as.numeric(ifelse(df$strand == "+", df$cdsend + 1, df$First))
  df$end = as.numeric(ifelse(df$strand == "+", df$First, df$cdsend - 1))   # DINGHAI'S CODE need to be checked by dinghai
  #df$end = as.numeric(ifelse(df$strand == "+", df$pos, df$cds_end - 1))   ###added by ww
  df = df[df$end >= df$start, ]
  df $ length <- df$end - df$start
  df = subset(df,select=-c(Strand, Chrom, First, Last))
  #df[, c("pos", "cds_end")] = NULL     ###added by ww
  
  makeGRangesFromDataFrame(df, keep.extra.columns=T)
}

####################functions for aUTR##########################################

create_aUTR_from_pAs = function(pA.df,geno){
  refdir="/scratch/ww346/analyze/APAanalyzer/REF/"
  df= read.table(paste0(refdir,"SRS.",geno,".conserved.FL.txt"),header=TRUE)
  df$strand = df$Strand
  df$chr = df$Chrom
  df$start = as.numeric(ifelse(df$strand == "+", df$First + 1, df$Last))
  df$end = as.numeric(ifelse(df$strand == "+", df$Last, df$First - 1))   # DINGHAI'S CODE need to be checked by dinghai
  #df$end = as.numeric(ifelse(df$strand == "+", df$pos, df$cds_end - 1))   ###added by ww
  df = df[df$end >= df$start, ]
  df $ length <- df$end - df$start
  df = subset(df,select=-c(Strand, Chrom, First, Last))
  
  makeGRangesFromDataFrame(df, keep.extra.columns=T)
}


#######distribution of #reads with no clustering
add_hg19_eCLIP_Stau2 = function(pA.df, sample){
  library(rtracklayer)
  # list.files("~/projects/fud/RIPiT/")
  # df_p = import(paste0(pA.df,sample,"_count.sorted.p.bw"), format = "BigWig") # hg19
  df_p = import(paste0(pA.df,sample,"_drpm.sorted.p.bw"), format = "BigWig") # hg19
  # summary(stau1_p$score)
  # hist(log10(stau1_p$score), breaks=50)
  strand(df_p)='+'

  # df_n = import(paste0(pA.df,sample,"_count.sorted.n.bw"), format = "BigWig") # hg19
  df_n = import(paste0(pA.df,sample,"_drpm.sorted.n.bw"), format = "BigWig") # hg19
  strand(df_n)='-'

  dfgr = c(df_p, df_n)

  #####calculate total count for rpm calculation

  total=sum(dfgr$score)

  noncoding = create_noncoding_from_pAs(pA.df,geno='hg19')
  utr5 = create_5utr_from_pAs(pA.df,geno='hg19')
  utr3 = create_3utr_from_pAs(pA.df,geno='hg19')

  intron=create_intron_from_pAs(pA.df,geno='hg19')

  CDS =create_cds_from_pAs(pA.df,geno='hg19')

  cUTR = create_cUTR_from_pAs(pA.df,geno='hg19')
  aUTR = create_aUTR_from_pAs(pA.df,geno='hg19')


  #################################3utr region
  ncread= as.data.frame(mergeByOverlaps(noncoding, dfgr)[, c("gene_symbol","score")]) %>%
    summarise(Stau2_eCLIP_count_noncoding = sum(score))

  utr5read= as.data.frame(mergeByOverlaps(utr5, dfgr)[, c("gene_symbol","score")]) %>%
    summarise(Stau2_eCLIP_count_5utr = sum(score))

  utr3read= as.data.frame(mergeByOverlaps(utr3, dfgr)[, c("gene_symbol","score")]) %>%
    summarise(Stau2_eCLIP_count_3utr = sum(score))

  intronread= as.data.frame(mergeByOverlaps(intron, dfgr)[, c("gene_symbol","score")]) %>%
    summarise(Stau2_eCLIP_count_intron = sum(score))

  CDSread= as.data.frame(mergeByOverlaps(CDS, dfgr)[, c("gene_symbol","score")]) %>%
    summarise(Stau2_eCLIP_count_cds = sum(score))


  cUTRread= as.data.frame(mergeByOverlaps(cUTR, dfgr)[, c("gene_symbol","score")]) %>%
    summarise(Stau2_eCLIP_count_cutr = sum(score))
  
  aUTRread= as.data.frame(mergeByOverlaps(aUTR, dfgr)[, c("gene_symbol","score")]) %>%
    summarise(Stau2_eCLIP_count_autr = sum(score))

  interread=total-utr5read[1,]-utr3read[1,]-intronread[1,]-CDSread[1,]-ncread[1,]

disdf = data.frame(region=c('5utr','CDS','3utr','intron','nc','intergenic','cUTR','aUTR'), 
  order=c(1,2,3,4,5,6,7,8),
  read_dis=c(utr5read[1,],CDSread[1,],utr3read[1,],intronread[1,],ncread[1,],interread,cUTRread[1,], aUTRread[1,]))

disdf
}





#######input rpm normalization for rep1 and rep2

sample = 'Hep_Stau2_1'
  
  max_RPM_df1=add_hg19_eCLIP_Stau2(pA.df, sample)
  print (max_RPM_df1)

  setwd('/scratch/ww346/analyze/Str_data/Stau2_eCLIP/plot')
  p1 <- ggplot(data=max_RPM_df1, aes(x=reorder(region,order),y=read_dis)) + 
    geom_bar(stat="identity")+
  #scale_y_continuous(limits = c(0, 1))+
    ggtitle(paste0("Stau2_eCLIP_distribution")) + labs(x="eCLIP rpm distribution", y="Sum of eCLIP rpm")+
  # geom_hline(yintercept=0, color='black')+
  # geom_vline(xintercept=0, color='black')+
  # theme(axis.title.x = element_blank(),
  #       panel.background = element_blank())
    theme_classic()

  jpeg(paste0(sample,'_eCLIP_rpm_distribution.tiff'), res=300, width=2000, height=2000, type="cairo")

    p1

  dev.off()

sample = 'Hep_Stau2_2'
  
  max_RPM_df2=add_hg19_eCLIP_Stau2(pA.df, sample)
  print (max_RPM_df2)
  

  setwd('/scratch/ww346/analyze/Str_data/Stau2_eCLIP/plot')
  p1 <- ggplot(data=max_RPM_df2, aes(x=reorder(region,order),y=read_dis)) + 
    geom_bar(stat="identity")+
  #scale_y_continuous(limits = c(0, 1))+
    ggtitle(paste0("Stau2_eCLIP_distribution")) + labs(x="eCLIP rpm distribution", y="Sum of eCLIP rpm")+
  # geom_hline(yintercept=0, color='black')+
  # geom_vline(xintercept=0, color='black')+
  # theme(axis.title.x = element_blank(),
  #       panel.background = element_blank())
    theme_classic()

  jpeg(paste0(sample,'_eCLIP_rpm_distribution.tiff'), res=300, width=2000, height=2000, type="cairo")

    p1

  dev.off()

  



