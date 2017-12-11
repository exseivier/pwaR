#	Batch pairwise alignment between orthologs of several species.
#
print("Loading Biostrings...")
library(Biostrings)

# OBJECTS
setClass("orthodb", representation(Data="data.frame"))
setClass("sequences", representation(filename="character", DNAss="list"))
setClass("Meta", representation(method="character", output="character", first.gen="numeric", second.gen="numeric", files="character", raac="character", db="character"))


# METHODS
setGeneric("parseData", function(args) standardGeneric("parseData"))
setMethod("parseData", signature("character"),
	function(args) {
		result <- list()
		data <- new("Meta",	method=args[length(args)-1], output=args[length(args)], first.gen=as.numeric(args[length(args)-3]), second.gen=as.numeric(args[length(args)-2]), files=args[1:(length(args)-6)], raac=args[length(args)-5], db=args[length(args)-4])
		database <- new("orthodb", Data=read.delim(data@db, header=T, sep="\t", row.names=1))
		result <- c(data, database)
		result
})

setGeneric("loadSeqs", function(data, ortho.db) standardGeneric("loadSeqs"))
setMethod("loadSeqs", signature("Meta", "orthodb"),
	function(data, ortho.db) {
		seqs <- list()
		filenames <- c()
		for (i in 1:length(data@files)) {
			file <- strsplit(data@files[i], "/")[[1]][length(strsplit(data@files[i], "/")[[1]])]
			DNAobject <- readDNAStringSet(data@files[i], format="fasta")
			target_genes <- as.vector(ortho.db@Data[,file])
			idx <- c()
			for (item in target_genes){
				idx <- c(idx, which(item == names(DNAobject)))
			}
			print(idx)
			DNAobject <- DNAobject[idx]
			seqs <- c(seqs, DNAobject)
			filenames <- c(filenames, file)
		}
		sequences <- new("sequences", filename=filenames, DNAss=seqs)
		names(sequences@DNAss) <- filenames
		sequences
})

setGeneric("PWA", function(orthologs, sequences, data) standardGeneric("PWA"))
setMethod("PWA", signature("orthodb", "sequences", "Meta"),
	function (orthologs, sequences, data) {
		# Creating combinations
		combination <- c()
		for (i in 1:(length(data@files)-1)) {
			a = i+1
			for (j in a:length(data@files)) {
				combination <- rbind(combination, c(i,j))
			}
		}
		print(combination)
		pids <- c()
		colns <- c()
		for (row in 1:length(combination[,1])) {
			pick.seqs1 <- combination[row, 1]
			pick.seqs2 <- combination[row, 2]
			alig <- pairwiseAlignment(seqs@DNAss[[pick.seqs1]][as.numeric(data@first.gen):as.numeric(data@second.gen)], seqs@DNAss[[pick.seqs2]][as.numeric(data@first.gen):as.numeric(data@second.gen)], type=data@method)
			pids <- cbind(pids, pid(alig))
			colns <- c(colns, paste(names(seqs@DNAss[pick.seqs1]), "_x_", names(seqs@DNAss)[pick.seqs2], sep=""))
			print(paste("Aligning ", names(seqs@DNAss)[pick.seqs1], " with ", names(seqs@DNAss)[pick.seqs2], sep=""))
		}
		rownames(pids) <- rownames(orthologs@Data)
		colnames(pids) <- colns
		pids
})
