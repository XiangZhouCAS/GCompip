#' To calculate the gene abundance in microbial communities (GAM,%).
#'
#' @param input_reads Please set the reads to include only the forward reads if the data is paired-end (PE).
#' @param result Please set the name for the result file.
#' @param threads Set the number of CPU threads, with the default value being 1.
#' @param Gene_db Please specify the directory for the target genes DIAMOND database.
#' @param USCG_db Please specify the directory for universal single-copy ribosomal genes database (e.g., 'Ribo_14.dmnd'). If you have already calculated the RPKM for these genes, you may instead specify the directory for the results (e.g., 'sample_name.USCG.hits.txt') to skip this step.
#' @param skip_fastp If you have already filtered the reads, you can set T to skip running fastp. The default is to run fastp.
#' @param min_length Set the minimum length required for filtering reads, the default is 50 bp.
#' @param run_seqkit If you have already counted the total number of reads using seqkit, you can specify the directory of the seqkit results (e.g., 'sample_name.all.reads.txt') to skip running seqkit. By default, seqkit will be executed.
#' @param filter_condition Please specify the filter_condition file. By default, the identity is set to 50 and the coverage to 80.
#' @param keep_samples By default, the temporary results will be deleted unless this parameter is setted.
#'
#' @return The gene abundance in microbial community (GAM,%)
#' @export
#'
#' @examples diamond_community(input_reads = "your/reads/data/sample_1.fastq",result = "sample_1",threads = 40,Gene_db = "your/database/target.genes.dmnd",USCG_db = "your/database/Ribo_14.dmnd",min_length = 75,filter_condition = "default")
diamond_GAM <- function(input_reads,result,threads = 1,
                        Gene_db,USCG_db,skip_fastp = F,
                        min_length = 50,run_seqkit = T,
                        filter_condition = "default",keep_samples = F){
  wd1 <- dirname(Gene_db)
  wd2 <- dirname(USCG_db)
  singleM <- paste0(wd2,"/",basename(USCG_db))
  diamond_db <- paste0(wd1,"/",basename(Gene_db))
  input_reads <- normalizePath(input_reads)
  setwd(dirname(result))
  outpath <- basename(result)
  fastp_output <- paste0(outpath,".filtered.fq.gz")
  diamond_db <- Gene_db
  singleM <- USCG_db
  min_length <- min_length
  system(sprintf("file %s > tmp.txt",singleM))
  system(sprintf("sed 's/.*://' -i tmp.txt"))
  system(sprintf("sed 's/ ASCII //g' -i tmp.txt"))
  tmp <- read.table("tmp.txt",sep = "\t")
  diamond_out <- paste0(basename(outpath),".hits.txt")
  singleM_out <- paste0(basename(outpath),".USCG.hits.txt")
  seqkit_out <- paste0(basename(outpath),".all.reads.txt")
  fastp <- sprintf("fastp -i %s -o %s --length_required %d -w %d > /dev/null 2>&1",
                   input_reads,fastp_output,min_length,threads)
  diamond <- sprintf("diamond blastx --db %s --query %s --out %s --threads %d --outfmt 6 slen stitle pident qcovhsp --max-target-seqs 1 --max-hsps 1 > /dev/null 2>&1",
                     diamond_db, fastp_output, diamond_out, threads)
  diamond_singleM <- sprintf("diamond blastx --db %s --query %s --out %s --threads %d --outfmt 6 slen stitle qcovhsp bitscore --max-target-seqs 1 --max-hsps 1 > /dev/null 2>&1",
                             singleM, fastp_output, singleM_out, threads)
  seqkit <- sprintf("seqkit stat %s > %s",
                    fastp_output,seqkit_out)
  diamond2 <- sprintf("diamond blastx --db %s --query %s --out %s --threads %d --outfmt 6 slen stitle pident qcovhsp --max-target-seqs 1 --max-hsps 1 > /dev/null 2>&1",
                      diamond_db, input_reads, diamond_out, threads)
  diamond_singleM2 <- sprintf("diamond blastx --db %s --query %s --out %s --threads %d --outfmt 6 slen stitle qcovhsp bitscore --max-target-seqs 1 --max-hsps 1 > /dev/null 2>&1",
                              singleM, input_reads, singleM_out, threads)
  seqkit2 <- sprintf("seqkit stat %s > %s",
                     input_reads,seqkit_out)
  if(skip_fastp == F){
    print("fastp is Running.")
    system(fastp)
    print("fastp is completed.")}else{
      print("Not run the fastp because you set the directory of filtered reads to skip it.")}
  if(skip_fastp == F){
    print("diamond is Running (target genes).")
    system(diamond)
    print("diamond is completed (target genes).")}else{
      print("diamond is Running (target genes).")
      system(diamond2)
      print("diamond is completed (target genes).")}
  if(skip_fastp == F){
    if(tmp$V1 != "text"){
      print("diamond is Running (USCGs).")
      system(diamond_singleM)
      print("diamond is completed (USCGs).")
    }else{print("Not count the RPKM of singleM marker genes.")}
  }else{if(tmp$V1 != "text" ){
    print("diamond is Running (USCGs).")
    system(diamond_singleM2)
    print("diamond is completed (target genes).")
  }else{print("Not count the RPKM of singleM marker genes.")}}
  if(skip_fastp == F){
    if(run_seqkit == T){
      print("seqkit is Running.")
      system(seqkit)
      print("seqkit is completed.")}else{
        print("Not run the seqkit because you set the directory of seqkit result")}
  }else{if(run_seqkit == T){
    print("seqkit is Running.")
    system(seqkit2)
    print("seqkit is completed.")
  }else{print("Not run the seqkit because you set the directory of seqkit result")}}
  if (!require(magrittr)) {
    install.packages("magrittr")
    library(magrittr)
  } else {
    library(magrittr)}
  if (!require(dplyr)) {
    install.packages("dplyr")
    library(dplyr)
  } else {
    library(dplyr)}
  if (!require(tidyr)) {
    install.packages("tidyr")
    library(tidyr)
  } else {
    library(tidyr)}
  if (!require(data.table)) {
    install.packages("data.table")
    library(data.table)
  } else {
    library(data.table)}
  if(run_seqkit == T){
    p <- read.table(seqkit_out,header = T)
  }else{p <- read.table(run_seqkit,header = T)}
  d1 <- fread(diamond_out,sep = "\t")
  colnames(d1)[c(1:4)] <- c("slen","gene","pident","qcovhsp")
  d1 <- d1%>%
    group_by_all()%>%
    count()
  d1$all_reads_nums <- p$num_seqs
  d1$all_reads_nums <- gsub(",","",d1$all_reads_nums)
  d1$all_reads_nums <- as.numeric(d1$all_reads_nums)
  d1$tmp_RPKM <- d1$n*10^9/(d1$slen*d1$all_reads_nums*3)
  d1$sample_name <- basename(outpath)
  d1$gene <- gsub(" .*","",d1$gene)
  if(filter_condition != "default"){
    f <- read.table(filter_condition,header = T,sep = "\t")
    rownames(f)[c(1:3)] <- c("gene","identity","coverage")
  }else{f <- data.frame(gene = unique(d1$gene),
                        identity = rep(50,length(unique(d1$gene))),
                        coverage = rep(80,length(unique(d1$gene))))}
  d1 <- left_join(d1,f,by = "gene")
  d1$tmp1 <- d1$pident-d1$identity
  d1$tmp2 <- d1$qcovhsp-d1$coverage
  d1 <- filter(d1,tmp1 >= 0 & tmp2 >= 0)
  d1 <- d1[,c(8,2,7)]
  d1 <- d1%>%
    group_by(gene)%>%
    mutate(RPKM = sum(tmp_RPKM))
  d1 <- unique(d1[,-3])
  if(tmp$V1 != "text"){
    d2 <- read.table(singleM_out,sep = "\t")
  }else{
    d2 <- read.table(singleM,sep = "\t")
  }
  colnames(d2)[c(1,2,3,4)] <- c("slen","sseq_id","qcovhsp","bitscore")
  d2 <- filter(d2,bitscore > 40 & qcovhsp > 80)
  d2 <- d2[,-c(3,4)]
  d2 <- d2%>%
    group_by_all()%>%
    count()
  d2$all_reads_nums <- p$num_seqs
  d2$all_reads_nums <- gsub(",","",d2$all_reads_nums)
  d2$all_reads_nums <- as.numeric(d2$all_reads_nums)
  d2$RPKM <- d2$n*10^9/(d2$slen*d2$all_reads_nums*3)
  d2$sample_name <- basename(outpath)
  d2 <- d2[,c(6,2,1,3,4,5)]
  d2$sseq_id <- gsub("-.*","",d2$sseq_id)
  d2 <- d2%>%
    group_by(sseq_id)%>%
    mutate(T_RPKM = sum(RPKM))
  d2 <- unique(d2[,c(1,2,7)])
  s_out <- paste0(basename(outpath),".USCG.hits.comb.txt")
  write.table(d2,s_out,quote = FALSE,sep = "\t",row.names = F)
  geo <- exp(mean(log(d2$T_RPKM)))
  d1$GAM <- d1$RPKM/geo*100
  d1$GAM <- ifelse(d1$GAM > 100,100,d1$GAM)
  r_out <- paste0(basename(outpath),".all.abd.txt")
  write.table(d1,r_out,quote = FALSE,sep = "\t",row.names = F)
  system(sprintf("mkdir %s",outpath))
  if(skip_fastp == F){
    system(sprintf("mv %s %s",fastp_output,outpath))
  }else{print("There is no need to move fastp result to outpath")}
  system(sprintf("mv %s %s",diamond_out,outpath))
  if(tmp$V1 != "text"){
    system(sprintf("mv %s %s",singleM_out,outpath))}else{
      print("There is no need to move singleM result to outpath")}
  if(run_seqkit == T){
    system(sprintf("mv %s %s",seqkit_out,outpath))}else{
      print("There is no need to move seqkit result to outpath")}
  system(sprintf("cat %s >> all.abd.txt",r_out))
  tmp <- read.table("all.abd.txt",sep = "\t",header = T)
  tmp <- filter(tmp,sample_name != "sample_name")
  tmp$GAM <- as.numeric(tmp$GAM)
  tmp$RPKM <- as.numeric(tmp$RPKM)
  write.table(tmp,"all.abd.txt",quote = F,sep = "\t",row.names = F)
  system(sprintf("mv %s %s",s_out,outpath))
  system(sprintf("mv %s %s",r_out,outpath))
  if(skip_fastp == F){
    system("rm \"fastp.html\" \"fastp.json\" \"tmp.txt\"")}else{
      system("rm  \"tmp.txt\"")}
  wd1 <- paste0(getwd(),"/",outpath)
  wd2 <- paste0(getwd())
  if(keep_samples == T){
    print("……")}else{if(file.exists("samples") == T){
      setwd(wd1)
      system(sprintf("mv %s ../samples",r_out))
      setwd(wd2)
      system(sprintf("rm %s -rf",outpath))}else{
        system(sprintf("mkdir samples"))
        setwd(wd1)
        system(sprintf("mv %s ../samples",r_out))
        setwd(wd2)
        system(sprintf("rm %s -rf",outpath))}}
  rpkm <- tmp%>%
    pivot_wider(id_cols = gene,
                values_from = "RPKM",
                names_from = "sample_name")
  rpkm <- as.data.frame(rpkm)
  rpkm[is.na(rpkm)] <- 0
  write.table(rpkm,"rpkm.abd.txt",
              sep = "\t",
              quote = F,row.names = F)
  com <- tmp%>%
    pivot_wider(id_cols = gene,
                values_from = "GAM",
                names_from = "sample_name")
  com <- as.data.frame(com)
  com[is.na(com)] <- 0
  com <- com[!duplicated(com[,c(1,2)]),]
  write.table(com,"GAM.abd.txt",
              sep = "\t",
              quote = F,row.names = F)
  if(!require(ggplot2)){
    install.packages("ggplot2")
    library(ggplt2)
  }else{library(ggplot2)}
  lv <- colnames(com)[-1]
  lv <- as.factor(lv)
  samplept <- com%>%
    pivot_longer(cols = !gene,values_to = "val",names_to = "sample")
  samplept$sample <- factor(samplept$sample,levels = lv)
  samplept$val <- as.numeric(samplept$val)
  hp <- ggplot(samplept,aes(sample,gene))+
    geom_tile(aes(fill = val),color = "white")+
    scale_fill_gradientn(colours = c("white","#C12554","#31357F"),
                         name = "GAM %")+
    theme(panel.background = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_text(size = 15,color = "black"),
          legend.position = "bottom",
          legend.title = element_text(vjust = 0.75),
          axis.text.y = element_text(size = 12,color = "black"),
          axis.ticks.y = element_line(color = "black"))
  w <- length(unique(samplept$sample))
  h <- length(unique(samplept$gene))
  ggsave(plot = hp,"Heatmap_GAM.pdf",width = w*2,height = h/4,limitsize = F)
  print("All Completed!")}
