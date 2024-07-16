#!/usr/bin/Rscript
# -*- encoding: utf-8 -*-

#@File    :   dada2.R
#@Time    :   2024/07/16 15:50:34
#@Author  :   Wanjin.Hu 
#@Version :   1.0
#@Contact :   wanjin.hu@outlook.com
#@Description :

options(warn = -1)
suppressMessages(library(getopt))
command <- matrix(c(
  "help", "h", 0, "logical", "Help docs",
  "indir", "i", 1, "character", "Input raw fq data dir.",
  "refdb", "r", 1, "character", "Input 16S ref taxonomy db",
  "outdir", "o", 1, "character", "Output dir"),
  byrow = TRUE, ncol = 5
)
args <- getopt(command)
if (! is.null(args$help) || 
      is.null(args$indir) || 
      is.null(args$refdb) || 
      is.null(args$outdir)) {
  cat(paste(getopt(command, usage = TRUE), "\n"))
  q(status = 1)
}

n_cores <- parallel::detectCores()
n_used_multithread <- as.integer(n_cores / 2)

dada2_run <- function(path) {
    '
    @功能: dada2运行函数
    @参数path: 原始数据存放文件夹
    '
    if (!dir.exists(args$outdir)){
      dir.create(args$outdir)
    }
    suppressMessages(library(dada2))
    fnFs <- sort(list.files(path, pattern= "_1.fq.gz", full.names = TRUE))
    fnRs <- sort(list.files(path, pattern= "_2.fq.gz", full.names = TRUE))
    sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
    filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fq.gz"))
    filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fq.gz"))
    names(filtFs) <- sample.names
    names(filtRs) <- sample.names
    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                         truncLen=c(220,220), ## Modify!
                         maxN=0, 
                         maxEE=c(6, 6), ## Modify!
                         truncQ=2, 
                         rm.phix=TRUE,
                         compress=TRUE, 
                         multithread=n_used_multithread) ## On Windows set multithread=FALSE
    ## Learn the error rate
    errF <- learnErrors(filtFs, multithread=TRUE)
    errR <- learnErrors(filtRs, multithread=TRUE)
    ## Sample Inference
    dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
    dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
    ## Merge paired reads
    mergers <- mergePairs(dadaFs, 
                          filtFs, 
                          dadaRs, 
                          filtRs, 
                          verbose=TRUE, 
                          maxMismatch = 2, ## Modify!
                          minOverlap = 10)  ## Modify!
    ## Construct sequence table
    seqtab <- makeSequenceTable(mergers)
    ## Remove chimeras
    seqtab.nochim <- removeBimeraDenovo(seqtab, 
                                        method="consensus", 
                                        multithread=TRUE, 
                                        verbose=TRUE)
    ## Process stat
    getN <- function(x) sum(getUniques(x))
    track <- cbind(out, 
                  sapply(dadaFs, getN), 
                  sapply(dadaRs, getN), 
                  sapply(mergers, getN), 
                  rowSums(seqtab.nochim))
    colnames(track)[3:ncol(track)] <- c("denoisedF", "denoisedR", "merged", "nonchim")
    rownames(track) <- sample.names
    write.csv(track,file = paste0(args$outdir,"/SequenceProcessStat.csv"),quote=F)
    ## Assign taxonomy
    taxa <- assignTaxonomy(seqtab.nochim,
                           refFasta = args$refdb,
                           multithread=TRUE)
    suppressMessages(library(tidyverse))
    df_out <- seqtab.nochim %>% as.data.frame()
    ## 生成详细物种组成表[Kindom,Phylum,Class,Order,Family,Genus,ASV,Samples...]
    df_out_t <- t(df_out)
    df1 <- as.data.frame(df_out_t)
    df1 <- df1 %>% mutate(ASV = paste0("ASV", c(1:nrow(df1)))) %>% select(ASV, everything())
    df3 <- merge(taxa, df1, by = "row.names", all = T)
    df3 <- df3[,-1]
    df3 <- df3 %>%
      mutate(
        Kingdom = paste0("k__", Kingdom),
        Phylum = paste0("p__", Phylum),
        Class = paste0("c__", Class),
        Order = paste0("o__", Order),
        Family = paste0("f__", Family),
        Genus = paste0("g__", Genus)
    )
    write.table(df3,file = paste0(args$outdir,"/SamplesTaxonomyComposition.txt"),sep="\t",col.names=T,row.names=F,quote=F)
    ## 生成可作为picrust2输入的文件
    seqs <- colnames(df_out)
    colnames(df_out) <- paste0("ASV",c(1:ncol(df_out)))
    ids <- paste0("ASV",c(1:ncol(df_out)))
    write.csv(df_out, paste0(args$outdir,"/ASV_table.csv"), quote = FALSE)
    seqinr::write.fasta(sequences = as.list(seqs), file.out = paste0(args$outdir,"/ASV.fasta"), names = ids)
}

## Run
dada2_run(path=args$indir)
