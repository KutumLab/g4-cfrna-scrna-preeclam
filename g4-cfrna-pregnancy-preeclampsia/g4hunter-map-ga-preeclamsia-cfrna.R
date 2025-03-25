# conda activate te-g4-preeclamsia
# cd /Users/rintukutum/Documents/office/ashoka/Proposals-2025-GrandChallenges/g4hunter
# R
rm(list=ls())
ga_pe_43h_genes_rf <- read.csv('../ga-preeclamsia-cfrna-43H-2022/41586_2021_4249_MOESM4_ESM.csv')
ga_pe_43h_genes_model <- read.csv('../ga-preeclamsia-cfrna-43H-2022/41586_2021_4249_MOESM4_ESM.csv')[,1]
ga_pe_43h_genes_ga_cor <- read.csv('../ga-preeclamsia-cfrna-43H-2022/ga-ga-cfrna-43H.csv')

ga_pe_43h_genes_cor <- unique(ga_pe_43h_genes_ga_cor$gene_symbol)

ga_ <- read.csv('../preeclamsia-cfrna-199-2022/sm3-ga-199-pe-de-genes-ensg.csv')[,1]
library(Biostrings)
library(org.Hs.eg.db)
library (BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#---- convert gene symbols to entrez ID
cols <- c("SYMBOL", "GENENAME", "ALIAS","ENSEMBL", "ENTREZID")
ga_pe_43h_gene_symbols <- unique(c(ga_pe_43h_genes_cor, ga_pe_43h_genes_model))
gene_maps <- select(org.Hs.eg.db, keys=ga_pe_43h_gene_symbols, columns=cols, keytype="SYMBOL")
gene_maps_ga_pe_43h_SYM_ID <- gene_maps[!duplicated(gene_maps$SYMBOL),
                                        c("SYMBOL","ENSEMBL","ENTREZID")]
idx_na_entrez <- is.na(gene_maps_ga_pe_43h_SYM_ID$ENTREZID)
entrez_43h <- gene_maps_ga_pe_43h_SYM_ID$ENTREZID[!idx_na_entrez]
source('function_G4Hunter.r')

ga_43h_rf <- gene_maps_ga_pe_43h_SYM_ID[gene_maps_ga_pe_43h_SYM_ID$SYMBOL %in% ga_pe_43h_genes_model,]
idx_na_entrez_rf <- is.na(ga_43h_rf$ENTREZID)
entrez_43h_rf <- ga_43h_rf$ENTREZID[!idx_na_entrez_rf]
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
blocks <- seq(1,length(entrez_43h_rf),100)
# library(doMC)
# registerDoMC(10)
run_block <- function(start, end){
  entrez_b <- entrez_43h_rf[start:end]
  entrez_b_out <- list()
  for(i in 1:length(entrez_b)){
    print(i)
    entrez_b_out[[i]] <- run_g4gene(entrez_b[i])
  }
  return(entrez_b_out)
}

#---b1 [done]
entrez_b1 <- entrez_43h_rf[blocks[1]:blocks[2]-1]
entrez_b1_out <- list()
for(i in 1:length(entrez_b1)){
  print(i)
  entrez_b1_out[[i]] <- run_g4gene(entrez_b1[i])
}
save(entrez_b1_out, file='entrez_b1_out.RData')
#----b 2
entrez_b2_out <- run_block(start = blocks[2], end = blocks[3]-1)
save(entrez_b2_out, file='entrez_b2_out.RData')
#----b 3 [done]
# entrez_b3_out <- run_block(start = blocks[3], end = blocks[4]-1)
save(entrez_b3_out, file='entrez_b3_out.RData')
#----b 4
entrez_b4_out <- run_block(start = blocks[4], end = blocks[5]-1)
save(entrez_b4_out, file='entrez_b4_out.RData')
#----b 5 [done]
entrez_b5_out <- run_block(start = blocks[5], end = blocks[6]-1)
save(entrez_b5_out, file='entrez_b5_out.RData')
#----b 6
entrez_b6_out <- run_block(start = blocks[6], end = blocks[7]-1)
save(entrez_b6_out, file='entrez_b6_out.RData')
#----b 7
entrez_b7_out <- run_block(start = blocks[7], end = length(entrez_43h_rf))
save(entrez_b7_out, file='entrez_b7_out.RData')
# registerDoMC(1)
#-----------
load('entrez_b2_out.RData');load('entrez_b3_out.RData');load('entrez_b4_out.RData')
load('entrez_b5_out.RData');load('entrez_b6_out.RData');load('entrez_b7_out.RData')
entrez_b_g4_rf <- c(entrez_b1_out, entrez_b2_out, entrez_b3_out, entrez_b4_out,
  entrez_b5_out, entrez_b6_out, entrez_b7_out)
names(entrez_b_g4_rf) <- entrez_43h_rf
save(entrez_b_g4_rf, file ='entrez_b_g4_rf.RData')
#-----------
rm(list=ls())
load('entrez_b_g4_rf.RData')

# G4H score >=1
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
idx_na <- sapply(entrez_b_g4_rf,function(x)is.na(x)[1])
entrez_b_g4_rf_ok <- entrez_b_g4_rf[!idx_na]

g4_score_rf <- lapply(entrez_b_g4_rf_ok,get_trans_score_per_gene)
g4_score_rf_df <- plyr::ldply(g4_score_rf)

g4_pa <- plyr::ldply(lapply(g4_score_rf, function(x){data.frame(table(x$g4))}))
