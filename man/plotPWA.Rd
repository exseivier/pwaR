% File man/plotPWA.Rd

\name{plotPWA}
\alias{plotPWA}
\title{Plots the pids data in a boxplot. separating the data of pids regarding to if genes are in 'raacg' vector or not}
\description{
This function plots the percentage of sequence identity in a boxplot separating the data regarding to if genes are in 'raacg' vector or not.
The resulting plot is stored in file named PWA_boxplot.pdf
}

\usage{
plotPWA(pids, raac.file)
}

\arguments{
	\item{pids}{numeric "matrix" with the percentage of sequence identity for every orthologou genes.}
	\item{raac.file}{The name of a file containing the regulated as a couple gene names list.}
}

\examples{
data <- parseData(args)
ortho.db <- data[[2]]
data <- data[[1]]
seqs <- loadSeqs(data=data, ortho.db=ortho.db)
pids <- PWA(ortho.db, seqs, data)
plotPWA(pids, data@raac)
}

\seealso{
	Meta, sequences, orthodb, loadSeqs, PWA, parseData
}

\author{
	Javier Montalvo
	\email{javier.montalvo@cinvestav.mx}
}
