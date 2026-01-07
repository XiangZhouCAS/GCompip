#' To calculate RPKM of universal single copy genes (RUSCG).
#'
#' @param input_reads Please set the directory of reads.
#' @param result Please set the directory of reads.
#' @param threads Setting the threads of CPU,default is 1.
#' @param USCG_db Please set the directory of universal single copy genes (USCGs) database.
#' @param skip_fastp If you have already filtered the reads, you can set this parameter to skip running fastp. The default is to run fastp.
#' @param min_length Set the minimum length required for filtering reads, the default is 50 bp.
#' @param run_seqkit If you have already counted the total number of reads using seqkit, you can specify the directory of the seqkit results (e.g., 'sample_name.all.reads.txt') to skip running seqkit. By default, seqkit will be executed.
#' @param keep_samples By default, the temporary results will be deleted unless this parameter is setted.
#' @return The rpkm abandance of universal single copy genes (RUSCG).
#' @export
#'
#' @examples ribo_rpkm(input_reads = "your/reads/data/sample_1.fastq",result = "sample_1",threads = 40,USCG_db = "your/database/Ribo_14.dmnd",skip_fastp = F,min_length = 75,run_seqkit = T,keep_samples = F)
uscg_rpkm <- function(input_reads,result,threads = 1,USCG_db,skip_fastp = F,
                      min_length = 50,run_seqkit = T,keep_samples = F){
  wd1 <- dirname(USCG_db)
  singleM <- paste0(wd1,"/",basename(USCG_db))
  setwd(dirname(result))
  res <- basename(result)
  singleM <- USCG_db
  singleM_out <- paste0(basename(res),".RUSCG.txt")
  fastp_output <- paste0(basename(res),".filtered.fq.gz")
  seqkit_out <- paste0(res,".all.reads.txt")
  seqkit <- sprintf("seqkit stat %s > %s",
                    fastp_output,seqkit_out)
  seqkit2 <- sprintf("seqkit stat %s > %s",
                     input_reads,seqkit_out)
  fastp <- sprintf("fastp -i %s -o %s --length_required %d -w %d",
                   input_reads,fastp_output,min_length,threads)
  diamond_singleM <- sprintf("diamond blastx --db %s --query %s --out %s --threads %d --outfmt 6 slen stitle qcovhsp bitscore --max-target-seqs 1 --max-hsps 1 > /dev/null 2>&1",
                             singleM, fastp_output, singleM_out, threads)
  diamond_singleM2 <- sprintf("diamond blastx --db %s --query %s --out %s --threads %d --outfmt 6 slen stitle qcovhsp bitscore --max-target-seqs 1 --max-hsps 1 > /dev/null 2>&1",
                              singleM, input_reads, singleM_out, threads)
  if(skip_fastp == F){
    print("fastp is Running.")
    system(fastp)
    print("fastp is completed.")}else{
      print("Not run fastp because you set the directory of filtered reads to skip it.")}
  if(skip_fastp == F){
    print("diamond is Running (USCGs).")
    system(diamond_singleM)
    print("diamond is completed (USCGs).")}else{
      print("diamond is Running (USCGs).")
      system(diamond_singleM2)
      print("diamond is completed (USCGs).")
    }
  if(skip_fastp == T){
    if(run_seqkit != T){
      print("Not run seqkit because you set the directory of seqkit result to skip it.")}else{
        print("seqkit is Running.")
        system(seqkit2)
        print("seqkit is completed.")}}else{
          print("seqkit is Running.")
          system(seqkit)
          print("seqkit is completed.")}
  if(!require(magrittr)){
    install.packages("magrittr")
    library(magrittr)
  }else{
    library(magrittr)}

  if(!require(dplyr)){
    install.packages("dplyr")
    library(dplyr)
  }else{
    library(dplyr)}

  d2 <- read.table(singleM_out,sep = "\t")
  colnames(d2)[c(1,2,3,4)] <- c("slen","sseq_id","qcovhsp","bitscore")
  d2 <- filter(d2,qcovhsp > 80 & bitscore > 40)
  d2 <- d2[,-c(3,4)]
  d2 <- d2%>%
    group_by_all()%>%
    count()
  s_out <- paste0(basename(res),".ribo.rpkm.txt")
  if(run_seqkit == T){
    p <- read.table(seqkit_out,header = T)
  }else{p <- read.table(run_seqkit,header = T)}
  d2$all_reads_nums <- p$num_seqs
  d2$all_reads_nums <- gsub(",","",d2$all_reads_nums)
  d2$all_reads_nums <- as.numeric(d2$all_reads_nums)
  d2$RPKM <- d2$n*10^9/(d2$slen*d2$all_reads_nums*3)
  d2$sample_name <- res
  d2 <- d2[,c(6,2,1,3,4,5)]
  d2$sseq_id <- gsub("-.*","",d2$sseq_id)
  d2 <- d2%>%
    group_by(sseq_id)%>%
    mutate(T_RPKM = sum(RPKM))
  d2 <- unique(d2[,c(1,2,7)])
  write.table(d2,s_out,quote = F,sep = "\t",row.names = F)
  system(sprintf("cat %s >> RUSCG.txt",s_out))
  total <- read.table("RUSCG.txt",header = T,sep = "\t")
  total <- filter(total,T_RPKM != "T_RPKM")
  write.table(total,"RUSCG.txt",quote = F,sep = "\t",row.names = F)

  system(sprintf("mkdir %s",res))
  if(run_seqkit == T){
    system(sprintf("mv %s %s",seqkit_out,res))
    system(sprintf("mv %s %s",s_out,res))
    system(sprintf("mv %s %s",singleM_out,res))}else{
      system(sprintf("mv %s %s",s_out,res))
      system(sprintf("mv %s %s",singleM_out,res))}
  if(skip_fastp == F){
    system(sprintf("mv %s %s",fastp_output,res))}else{
      print("There is no fastp result.")}
  wd1 <- paste0(getwd(),"/",res)
  wd2 <- paste0(getwd())
  if(keep_samples == T){
    print("……")}else{if(file.exists("samples") == T){
      setwd(wd1)
      system(sprintf("mv %s ../samples",s_out))
      setwd(wd2)
      system(sprintf("rm %s -rf",res))}else{
        system(sprintf("mkdir samples"))
        setwd(wd1)
        system(sprintf("mv %s ../samples",s_out))
        setwd(wd2)
        system(sprintf("rm %s -rf",res))}}
  print("All Completed!")}
