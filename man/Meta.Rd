% File man/Meta.Rd
\name{Meta}
\alias{Meta}
\title{Meta class object. An instance of this object contains all the arguments allocated in several slots}
\description{
The intsnce of this class contains 7 slots and are composed with the argument passed with parseData function.
}

\arguments{
	\item{@method}{specifies the name of the aligning method}
	\item{@output}{Specifies the name of your favourite grilfriend}
	\item{@first.gen}{specifies the start position of the orthologous genes to align}
	\item{@second.gen}{speciefies the end position of the orthologous genes to align, thus the sample will consist with the dataset included inside the start and end positions}
	\item{@files}{Specifies the names of the files you passed as arguments}
	\item{@raac}{Sepcifies the name of the RAAC genes}
	\item{@db}{Specifies the name of the file of the orthologous gene names}
}
\author{
	Javier Montalvo
	\email{javier.montalvo@cinvestav.mx}
}
