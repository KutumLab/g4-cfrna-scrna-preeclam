# conda activate te-g4-preeclamsia
# cd /Users/rintukutum/Documents/office/ashoka/Proposals-2025-GrandChallenges/g4hunter
# R
ga_ <- read.csv('../preeclamsia-cfrna-199-2022/sm3-ga-199-pe-de-genes-ensg.csv')[,1]
library(Biostrings)
library(org.Hs.eg.db)
library (BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#---- convert gene symbols to entrez ID
cols <- c("SYMBOL", "GENENAME", "ALIAS","ENSEMBL", "ENTREZID")
gene_maps <- select(org.Hs.eg.db, keys=ga_, columns=cols, keytype="SYMBOL")
gene_maps_ga_pe <- gene_maps[!duplicated(gene_maps$SYMBOL),
                                        c("SYMBOL","ENSEMBL","ENTREZID")]
idx_na_entrez <- is.na(gene_maps_ga_pe$ENTREZID)
entrez_pe <- gene_maps_ga_pe$ENTREZID[!idx_na_entrez]
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

entrez_pe_out <- list()
for(i in 1:length(entrez_pe)){
  print(i)
  entrez_pe_out[[i]] <- run_g4gene(entrez_pe[i])
}
save(entrez_pe_out,file='entrez_pe_out.RData')

rm(list=ls())
library(Biostrings)
load('entrez_pe_out.RData')
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
idx_na <- sapply(entrez_pe_out,function(x)is.na(x)[1])
entrez_pe_ok <- entrez_pe_out[!idx_na]

g4_score_pe <- lapply(entrez_pe_ok,get_trans_score_per_gene)
g4_score_pe_df <- plyr::ldply(g4_score_pe)

png('fig/4A-g4-propensity-scores-within-trans-cfRNA-preclamp.png',
    width = 1200, height = 500, res = 120)
hist(g4_score_pe_df$score, 100)
dev.off()

g4_pe <- plyr::ldply(lapply(g4_score_pe, function(x){data.frame(table(x$g4))}))

png('fig/4B-g4-propensity-and-G4-counts-per-trans-cfRNA-preclamp.png',
    width = 400, height = 700, res = 120)

boxplot(log2(Freq)~Var1, data=g4_pe,
        xlab='G4 Propensity',
        ylab = 'log2(Frequency of G4 across genes in cfRNA in PE)'
)
dev.off()

