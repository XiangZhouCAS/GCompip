#' To convert RPKM to GAM of genes from geneset.
#'
#' @param input_geneset Please set the directory of geneset RPKM abundance file.
#' @param input_uscg_rpkm Please set the directory of USCGs' RPKM file, e.g.,RUSCG.txt.
#' @param output Please set the directory of output.
#'
#' @return The gene abundance in microbial community (GAM,%)
#' @export
#'
#' @examples uscg_rpkm2GAM(input_geneset = "your/data/geneset.RPKM.txt",input_uscg_rpkm = "your/database/RUSCG.txt",output = "geneset.GAM.txt")
uscg_rpkm2GAM <- function(input_geneset,input_uscg_rpkm,output){
  geneset <- input_geneset
  rpkm <- input_ucsg_rpkm
  all_data <- read.table(rpkm,header = T,sep = "\t")
  geneset <- read.table(geneset,header = T,
                        sep = "\t",check.names = F)
  if(!require(dplyr,quietly = TRUE)){
    install.packages("dplyr")
    library(dplyr)}else{
      library(dplyr)}
  if(!require(tidyr,quietly = TRUE)){
    install.packages("tidyr")
    library(tidyr)}else{
      library(tidyr)}
  all_data <- all_data%>%
    group_by(sample_name)%>%
    mutate(geomean = exp(mean(log(T_RPKM))))
  all_data <- unique(all_data[,c(1,4)])
  all_data <- all_data[match(colnames(geneset[,-1]),
                             all_data$sample_name),]
  for(i in 1:nrow(all_data)){
    ribo_rpkm <- as.numeric(all_data[i,2])
    geneset[,i+1] <- geneset[,i+1]*100/ribo_rpkm}
  df <- data.frame(lapply(geneset[,-1],
                          function(x)ifelse(x > 100, 100, x)),check.names = F)
  res <- basename(output)
  colnames(geneset)[1] <- "GeneID"
  df$GeneID <- geneset$GeneID
  df$GeneID <- geneset$GeneID
  df <- df[,c(ncol(df),1:ncol(df)-1)]
  system(sprintf("mkdir %s",res))
  setwd(dirname(result))
  system(sprintf("rm -rf %s",basename(res)))
  write.table(df,res,sep = "\t",quote = F,row.names = F)
  print("All Completed!")}
