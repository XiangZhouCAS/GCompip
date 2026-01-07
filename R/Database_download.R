#' Download database.
#'
#' @param db_filepath Download the USCG database and Hydrogenase, hydrogen metabolism terminal enzyme database.
#' @param db_type Select the type of database to download, default is "USCG", Other types include "ter", "hyd" and "all".
#' @return The database
#' @export
#'
#' @examples Database_download(db_filepath = "your/database/filepath",db_type = "USCG")
Database_download <- function(db_filepath = getwd(),db_type = "USCG"){
  if(db_type == "USCG"){
    wd1 <- getwd()
    wd2 <- db_filepath
    setwd(wd2)
    system(sprintf("wget -c https://github.com/XiangZhouCAS/GCompip/raw/refs/heads/main/database/Ribo_14.dmnd"))
    print("USCG database has been downloaded.")
    setwd(wd1)}
  if(db_type == "hyd"){
    wd1 <- getwd()
    wd2 <- db_filepath
    system(sprintf("wget -c https://github.com/XiangZhouCAS/GCompip/raw/refs/heads/main/database/hyd_id-name.script.txt"))
    system(sprintf("wget -c https://github.com/XiangZhouCAS/GCompip/raw/refs/heads/main/database/hyddb.all.dmnd"))
    print("Hydrogense database and hyd_id-name.script.txt have been downloaded.")
    setwd(wd1)}
  if(db_type == "ter"){
    wd1 <- getwd()
    wd2 <- db_filepath
    system(sprintf("wget -c https://github.com/XiangZhouCAS/GCompip/raw/refs/heads/main/database/ter.dmnd.gz"))
    print("Terminal database has been downloaded.")
    setwd(wd1)}
  if(db_type == "all"){
    wd1 <- getwd()
    wd2 <- db_filepath
    system(sprintf("wget -c https://github.com/XiangZhouCAS/GCompip/raw/refs/heads/main/database/Ribo_14.dmnd"))
    system(sprintf("wget -c https://github.com/XiangZhouCAS/GCompip/raw/refs/heads/main/database/hyd_id-name.script.txt"))
    system(sprintf("wget -c https://github.com/XiangZhouCAS/GCompip/raw/refs/heads/main/database/hyddb.all.dmnd"))
    system(sprintf("wget -c https://github.com/XiangZhouCAS/GCompip/raw/refs/heads/main/database/ter.dmnd.gz"))
    print("All databases have been downloaded.")
    setwd(wd1)}}
