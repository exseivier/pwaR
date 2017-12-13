% File man/parseData.Rd
\name{parseData}
\alias{parseData}
\title{Takes the arguments needed}
\description{
Takes the arguments needed to create the objects and to perform the pairwise alignment with the PWA function.
It returns a List with two objects, the first is a 'Meta' object and contains the arguments you passed (i.e. from bash),
the second is a 'orthodb' object with the names of the orthologous genes
}
\usage{
data <- parseData(args)
}
\arguments{
	\item{args}{A character vector with the arguments nedded in the folowing order. The
	first arguments should be the name of the fasta files. Then,  regulated as a couple genes list file,
	then the orthologous file, then, the next two arguments should be the position of start number and
	end number of the subset of genes you want to align, then, the next argument should specify the aligning
	method and the last argument should be the name of favourite grilfriend. Note that args character vector
	could be obtained from bash with the 'commandArgs(trailingOnly=T)' function.
	}
}
\seealso{sequences, orthodb, Meta}
\examples{
args <- c("file1.fa", "file2.fa", "file3.fa", "raac.txt", "ortho.db", "1", "2300", "global", "Veronica")
data <- parseData(args)
ortho.db <- data[[2]]
data <- data[[1]]
}

