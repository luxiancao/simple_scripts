# Date：2020.11.23
# Author：Yusy
# Convert MCC TREE TO NEWICK FORMAT

rm(list=ls())
setwd("D:/Code/R/MCC_to_NWK/")

# needs package treeio
library(treeio)

# FUNCTION DEFINITION

# read the MCC tree file from TreeAnnotator
fname = 'H7N9_HA_cd99_GMRF_MCC_grid.trees'
doLadderize=FALSE
tol=1e-8
root0=TRUE

#lines 	<- readLines(fname)

tr_mcc <- treeio::read.beast(fname)
tr <- treeio::read.nexus(fname)
tr_r <- treeio::read.nexus(fname)

if (doLadderize) {
  tr		<- multi2di(tr)
  tr		<- ladderize(tr, right=FALSE)
  einds		<- which(tr$edge.length < tol)
  if (length(einds) > 0) {
    tr$edge.length[einds] <- tol
  }
}

# extract proper tip names
#tips 		<- unlist(apply(as.matrix(tr$tip.label), 1, strsplit, "\\["))
#tinds 	<- match(tips, taxaTbl[,1])
#tr$tip.label <- taxaTbl[tinds,2]


# extract posterior support
#tr_mcc <- treeio::read.beast('h3n2_1027_beast_mcc.tre')
#tr_nwk <- treeio::read.newick('h3n2_1027_beast_mcc.tre.nwk')
#getNodeNum(tr_mcc)
#getNodeNum(tr_nwk)

for (i in 1:length(tr_mcc@data$posterior)){
  id <- as.numeric(tr_mcc@data$node[i])
  tr_r$node.label[id] <- tr_mcc@data$posterior[i]
}
tr$node.label = na.omit(tr_r$node.label)

pp		<- as.numeric(tr$node.label)
inds		<- which(!is.finite(pp))
if (length(inds) > 0) {
  pp[inds]	<- 0
}
tr$node.label <- format(pp, digits=4)

#tr$node.label[1] <- "0.000"

trString = treeio::write.tree(tr)

# this string should end with support:root length
trString	<- gsub(";","",trString)
trString	<- paste(trString, ":0.0;", sep="")
trString	<- write.tree(tr,file='')

if (root0) {
  # set support at root to 0 for cluster picker
  tr$node.label[1] <- "0.000"
  trString	<- write.tree(tr,file='')
  
  # this string should end with support:root length
  trString	<- gsub(";","",trString)
  trString	<- paste(trString, ":0.0;", sep="")
} else {
  # remove last support value - this is not usually a good idea
  # so not doing it now
  
  trString	<- write.tree(tr)
  
  #lastEl	<- strsplit(trString, "\\)")[[1]]
  #lastEl	<- lastEl[length(lastEl)]
  #trString	<- gsub(lastEl, ";", trString, fixed=TRUE)
}

outName <- paste(fname, ".posterior.nwk", sep="")
write(trString, file=outName)

print(paste("Newick tree written to",outName))

# uncomment the line below if you actually want the tree object returned to the R-workspace
# return( tr )


