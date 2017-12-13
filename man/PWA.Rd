% File man/PWA.Rd

\name{PWA}
\alias{PWA}
\title{Performs the pairwise alignments between the orthologous sequences and for each files pair match}
\description{
Performs the pairwise alignments between the orthologous sequences and for each files pair match.
This function performs this task using the pairwiseAlignment function of the Biostrings R package, and
returns a "matrix" with the percentage of sequence identity of every alignment comparison for every gene.
}

\usage{
pids <- PWA(ortho.db, seqs, data)
}

\arguments{
	\item{ortho.db}{'orthodb' class object. this is a database of the orthologous gene names.}
	\item{seqs}{'sequences' class object with the sequences. This object is created with loadSeqs function.}
	\item{data}{'Meta' class object with the arguments passed with the parseData function.}
}
\examples{
data <- parseData(args)
ortho.db <- data[[2]]
data <- data[[1]]
seqs <- loadSeqs(data=data, ortho.db=orthodb)
pids <- PWA(ortho.db, seqs, data)
}
\seealso{
	Meta, sequences, orthodb, parseData, loadSeqs
}
\author{
	Javier Montalvo
	\email{javier.montalvo@cinvestav.mx}
}
