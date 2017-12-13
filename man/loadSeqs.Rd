% File man/loadSeqs.Rd

\name{loadSeqs}
\alias{loadSeqs}
\title{Load sequences from specified fasta files}
\description{
This function takes the names of fasta files and extracts the sequences and stores them in a 'sequences' object
}
\usage{
seqs <- loadSeqs(data, ortho.db)
}
\arguments{
This function takes 2 arguments. The first is a 'Meta' object and the second is a 'orthodb' object
	\item{data}{this argument should be a 'Meta' object and contains the arguments you passed with parseData}
	\item{ortho.db}{This object should be 'orthodb' class and is produced by parseData function. It contains the
	names of the orthologous genes}
}
\examples{
data <- parseData(args)
ortho.db <- data[[2]]
data <- data[[1]]
seqs <- loadSeqs(data=data, ortho.db=ortho.db)
}
\seealso{Meta, sequences, parseData, orthodb}
\author{
	Javier Montalvo
	\email{javier.montalvo@cinvestav.mx}
}

