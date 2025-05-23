% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uscg_rpkm2com.R
\name{uscg_rpkm2GAM}
\alias{uscg_rpkm2GAM}
\title{To convert RPKM to GAM of genes from geneset.}
\usage{
uscg_rpkm2GAM(input_geneset, input_uscg_rpkm, output)
}
\arguments{
\item{input_geneset}{Please set the directory of geneset RPKM abundance file.}

\item{input_uscg_rpkm}{Please set the directory of USCGs' RPKM file, e.g.,RUSCG.txt.}

\item{output}{Please set the directory of output.}
}
\value{
The gene abundance in microbial community (GAM,\%)
}
\description{
To convert RPKM to GAM of genes from geneset.
}
\examples{
uscg_rpkm2GAM(input_geneset = "your/data/geneset.RPKM.txt",input_uscg_rpkm = "your/database/RUSCG.txt",output = "geneset.GAM.txt")
}
