#	Batch pairwise alignment between orthologs of several species.
#

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
		rownames(pids) <- rownames(orthologs@Data)[data@first.gen:data@second.gen]
		colnames(pids) <- colns
		pids
})

setGeneric("testDup", function(pids, raacg) standardGeneric("testDup"))
setMethod("testDup", signature("matrix", "character"),
	function(pids, raacg) {
		dup <- ifelse(rownames(pids) %in% raacg, 1, 0)
		pids <- cbind(pids, dup)
		pids
	})

setGeneric("plotPWA", function(pids, raac.file) standardGeneric("plotPWA"))
setMethod("plotPWA", signature("matrix", "character"),
	function(pids, raac.file){
		if(is.null(rownames(pids)) || is.null(colnames(pids))){
			stop("Rownames or columnames are empty!")
		}
		raacg <- readLines(raac.file)

		if(is.null(raacg)){
			stop("RAAC genes list is emnpty!")
		}
		if (colnames(pids)[length(colnames(pids))] != "dup"){
			pids <- testDup(pids, raacg)
		}

		len_col_pids <- length(pids[1,])
		pids_raac <- pids[pids[,len_col_pids] == 1,]
		pids_no_raac <- pids[pids[,len_col_pids] == 0,]
		pdf("PWA_boxplot.pdf", width=9, height=6)
		par(mfrow=c(1,2))
		boxplot(pids_raac[,-length(pids_raac[1,])], main="RAAC genes", names=colnames(pids)[-length(colnames(pids))], ylab="Percentage of sequence identity")
		boxplot(pids_no_raac[,-length(pids_no_raac[1,])], main="No-RAAC genes", names=colnames(pids)[-length(colnames(pids))], ylab="Percentage of sequence identity")
		dev.off()
	})

setGeneric("javoxon.test", function(pids, raac.file, boot, fun) standardGeneric("javoxon.test"))
setMethod("javoxon.test", signature("matrix", "character", "numeric", "character"),
	function(pids, raac.file, boot, fun) {
		if(is.null(rownames(pids)) || is.null(colnames(pids))){
			stop("Rownames or columnames are empty!")
		}
		raacg <- readLines(raac.file)

		if(is.null(raacg)){
			stop("RAAC genes list is emnpty!")
		}
		if (colnames(pids)[length(colnames(pids))] != "dup"){
			pids <- testDup(pids, raacg)
		}
		
		len_col_pids <- length(pids[1,])
		pids_raac <- pids[pids[,len_col_pids] == 1,-length(colnames(pids))]
		pids_no_raac <- pids[pids[,len_col_pids] == 0,-length(colnames(pids))]
		if(length(pids_raac[,1]) > length(pids_no_raac[,1])) {
			sample_size <- length(pids_no_raac[,1]) * 0.5
		}
		else {
			sample_size <- length(pids_raac[,1]) * 0.5
		}
		sp_raac <- c()
		sp_no_raac <- c()
		for(i in 1:boot) {
			sample_pids_raac <- pids_raac[sample((pids_raac[,1]), size=sample_size, replace=F),-length(pids[1,])]
			sample_pids_no_raac <- pids_no_raac[sample(pids_no_raac[,1], size=sample_size, replace=F),-length(pids[1,])]
			sp_raac <- rbind(sp_raac, apply(X=sample_pids_raac, MARGIN=2, fun))
			sp_no_raac <- rbind(sp_no_raac, apply(X=sample_pids_no_raac, MARGIN=2, fun))
		}
		pdf("javoxon.test.pdf", width=8, height=6)
		par(mfrow=c(1,2))
		boxplot(sp_raac, main="RAAC genes", names=colnames(pids)[-length(colnames(pids))], ylab="Percentage of sequence identity")
		boxplot(sp_no_raac, main="No RAAC genes", names=colnames(pids)[-length(colnames(pids))], ylab="Percentage of sequence identity")
		dev.off()
	})
