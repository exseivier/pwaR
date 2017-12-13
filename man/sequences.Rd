% File man/sequences.Rd
\name{sequences}
\alias{sequences}
\title{Creates a 'sequences' class object used to store the sequences of every file}
\description{
Creates a 'sequences' class object used to store the sequences of every file.
Sequences of each file are stored in a 'DNAStringSet' object, and every 'DNAStringSet' object
is attached to a List object
}
\usage{
This object is created in loadSeqs function. It has two main slots.
}
\arguments{
	\item{Slot1}{'character' vector containing the names of the sequences files}
	\item{Slot2}{'list' object contains the DNAStringSet objects storing sequences set of every file}
}
\author{
Javier Montalvo. \email{javier.montalvo@cinvestav.mx}
}


