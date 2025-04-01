# conda activate te-g4-preeclamsia
# cd /Users/rintukutum/Documents/office/ashoka/Proposals-2025-GrandChallenges/g4hunter
# R

rm(list=ls())
ga_pe_43h_genes_ga_cor <- read.csv('../ga-preeclamsia-cfrna-43H-2022/ga-ga-cfrna-43H.csv')
ga_pe_43h_genes_cor <- unique(ga_pe_43h_genes_ga_cor$gene_symbol)

library(Biostrings)
library(org.Hs.eg.db)
library (BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#---- convert gene symbols to entrez ID
cols <- c("SYMBOL", "GENENAME", "ALIAS","ENSEMBL", "ENTREZID")
ga_pe_43h_gene_symbols <- ga_pe_43h_genes_cor
gene_maps <- select(org.Hs.eg.db, keys=ga_pe_43h_gene_symbols, columns=cols, keytype="SYMBOL")
gene_maps_ga_pe_43h_SYM_ID <- gene_maps[!duplicated(gene_maps$SYMBOL),
                                        c("SYMBOL","ENSEMBL","ENTREZID")]
# 9582
idx_na_entrez <- is.na(gene_maps_ga_pe_43h_SYM_ID$ENTREZID)
entrez_43h_af <- gene_maps_ga_pe_43h_SYM_ID$ENTREZID[!idx_na_entrez]
# 9366
source('function_G4Hunter.r')
trans2g4 <- function(input){
  #input <- exon.seqs.by.transcript[1]
  transcript <- as.character(input)
  k <- 25
  masked=5
  # chr_G4hk <- G4runmean(chr_trans[[1]])
  exon_bind <- paste(as.character(transcript),collapse='')
  chr_G4hk <- G4runmean(exon_bind)
  tgenome <- G4translate(Rle(strsplit(as.character(exon_bind),NULL)[[1]]))
  
  j <- 1
  chrCh <- Views(chr_G4hk, chr_G4hk<=(-j))
  chrGh <- Views(chr_G4hk, chr_G4hk>=j)
  
  IRC <- IRanges(start=start(chrCh),end=(end(chrCh)+k-1))
  IRG <- IRanges(start=start(chrGh),end=(end(chrGh)+k-1))
  nIRC <- reduce(IRC)
  nIRG <- reduce(IRG)
  
  
  # overlapping results on the same strand are fused
  nchrCh <- Views(chr_G4hk,start=start(nIRC),end=(end(nIRC)-k+1))
  nchrGh <- Views(chr_G4hk,start=start(nIRG),end=(end(nIRG)-k+1))
  nscoreC <- signif(mean(Views(tgenome,nIRC)),3)
  nscoreG <- signif(mean(Views(tgenome,nIRG)),3)
  
  nstraC <- Rle(rep('-',length(nIRC)))
  nstraG <- Rle(rep('+',length(nIRG)))
  nhlC <- Rle(rep(j,length(nIRC)))
  nhlG <- Rle(rep(j,length(nIRG)))
  nkC <- Rle(rep(k,length(nIRC)))
  nkG <- Rle(rep(k,length(nIRG)))
  maskC <- Rle(rep(masked,length(nIRC)))
  maskG <- Rle(rep(masked,length(nIRG)))
  
  if (length(nIRC)==0)
  {        
    nxC <- GRanges()
  }else{
    nxC <- GRanges(seqnames=Rle(exon_bind),
                   ranges=nIRC,
                   strand=nstraC,
                   score=nscoreC,hl=nhlC,k=nkC,mask=maskC)
  }
  
  if (length(nIRG)==0)
  {        
    nxG <- GRanges()
  }else{          
    nxG <- GRanges(seqnames=Rle(exon_bind),
                   ranges=nIRG,
                   strand=nstraG,
                   score=nscoreG,hl=nhlG,k=nkG,mask=maskG)
  }
  nx <- c(nxC,nxG)
  return(nx)
}
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
transcripts.txdb <- transcriptsBy(txdb, by="gene")
run_g4gene <- function(entrezID){
  if(any(names(transcripts.txdb) == entrezID)){
    transcripts.grL <- transcripts.txdb[entrezID]
    transcript.names <- mcols(unlist(transcripts.grL))$tx_name  #
    message("No. of trsncripts: ", length(transcript.names))
    #--------- INTRON
    # intronsByTx.grL <- intronsByTranscript(
    #   txdb,use.names=TRUE)[transcript.names]
    # 
    # intron.seqs <- getSeq(Hsapiens, unlist(intronsByTx.grL))
    # intron.seqs.by.transcript <- relist(
    #   intron.seqs,
    #   skeleton=intronsByTx.grL)
    # print(elementLengths(intron.seqs.by.transcript))
    
    #--------- EXON
    exonsByTx.grL <- exonsBy(
      txdb, by="tx",
      use.names=TRUE)[transcript.names]
    
    # print(elementLengths(exonsByTx.grL))
    
    exons.seq <- getSeq(Hsapiens, unlist(exonsByTx.grL))
    exon.seqs.by.transcript <- relist(exons.seq, skeleton=exonsByTx.grL)
    # print(elementLengths(exon.seqs.by.transcript))
    #-------------- transcript to G4 Propensity
    g4_out <- lapply(exon.seqs.by.transcript,trans2g4)
  }else{
    g4_out <- NA
  }
  return(g4_out)
}

library(doMC)
registerDoMC(10)
entrez_43h_af_out <- foreach(i=1:length(entrez_43h_af))%dopar%{
  print(i)
  run_g4gene(entrez_43h_af[i])
}
registerDoMC(1)
save(entrez_43h_af_out, file='entrez_43h_af_out.RData')
#-------
rm(list=ls())
load('entrez_43h_af_out.RData')
get_trans_score_per_gene <- function(x){
  if(length(x) == 0){
    score_df<-NA
  }else{
    score <- lapply(x,function(x){
      df_ <- data.frame(score = score(x))
      df_$g4 <- ifelse(score(x)>=1,1,0)
      return(df_)
    })
    score_df <- plyr::ldply(score, .id = 'trans.name')
  }
  return(score_df)
}
idx_na <- sapply(entrez_43h_af_out,function(x)is.na(x)[1])
entrez_43h_af__ok <- entrez_43h_af_out[!idx_na]

g4_score_af <- lapply(entrez_43h_af__ok,get_trans_score_per_gene)
g4_score_af_df <- plyr::ldply(g4_score_af)

g4_af <- plyr::ldply(lapply(g4_score_af, function(x){data.frame(table(x$g4))}))

save(g4_score_af_df, file='g4_score_af_df.RData')

rm(list=ls())
load('entrez_43h_af_out.RData')
load('g4_score_af_df.RData')
library(ggplot2)
library(ggpubr)
dir.create('figures', showWarnings = FALSE)
png('fig/ga-24h-adult-fetal-violin-G4score-G4hunter-RM2022.png',
    width = 600, 
    height = 400,res = 150)
ggplot(g4_score_af_df, aes(x=as.factor(g4), y=score)) +
  geom_violin(aes(col=as.factor(g4), fill=as.factor(g4))) +
  coord_flip() +
  scale_color_manual(
    'G4',
    values = c('#27548A', '#DDA853'),
    labels = c('absent(N=316110)', 'present (N=60471; 16%)')) +
  scale_fill_manual(
    'G4',
    values = c('#27548A80', '#DDA85380'),
    labels = c('absent(N=316110)', 'present (N=60471; 16%)')) +
  ylab('Predicted G4 score with G4hunter') +
  ggtitle('0.37 million regions of 9582 gene signatures\n[Rasmussen M, et al 2022]') +
  scale_x_discrete(labels = NULL) +
  xlab('') +
  theme_bw() +
  theme(legend.position = c(0.3,0.65))
dev.off()
  


