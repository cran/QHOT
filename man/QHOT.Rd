% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QHOT.R
\name{QHOT}
\alias{QHOT}
\title{QTL Hotspot Detection}
\usage{
QHOT(DataQTL, DataCrop, ScanStep, NH, NP)
}
\arguments{
\item{DataQTL}{a data-frame of values for QTL information including the trait names, which chromosomes localized, the left and right flanking marker positions of QTLs for the first to fourth columns, respectively.}

\item{DataCrop}{a data-frame of values for chromosome information consisting of the names, center positions and lengths of chromosomes for the first to third columns, respectively.}

\item{ScanStep}{a value for the length of every bin.}

\item{NH}{a value for the number of spurious hotspots in the proposed method.}

\item{NP}{a value for permutation times to calculate the threshold.}
}
\description{
This function produces both the numerical and graphical summaries of the QTL hotspot detection in the genomes that are available on the worldwide web including the flanking markers of QTLs.
Man-Hsia Yang, Dong-Hong Wu, Chen-Hung Kao. 2019. A Statistical Procedure for Genome-wide Detection of QTL Hotspots Using Public Databases with Application to Rice. G3-Genes Genom Genet DOI: 10.1534/g3.118.200922.
}
\examples{
Trait<-paste("t",sample(1:9,100,replace=TRUE,prob=c(2,rep(1,8))/10),sep="")
chr<-1
L<-sample(seq(0,90,by=10),100,replace=TRUE,prob=c(0.5,0.5,5.5,rep(0.5,7))/10)
R<-L+sample(c(0.5,1,5,10,50),100,replace=TRUE)
R[R>100]<-100
DataQTL.t<-data.frame(Trait,chr,L,R)
DataCrop.t<-data.frame(chr=1,center=75,length=100)
QHOT(DataQTL.t, DataCrop.t, ScanStep=0.5, NH=1, NP=1000)
}
\references{
Man-Hsia Yang, Dong-Hong Wu, Chen-Hung Kao. 2019. A Statistical Procedure for Genome-wide Detection of QTL Hotspots Using Public Databases with Application to Rice. G3-Genes Genom Genet DOI: 10.1534/g3.118.200922.
}
