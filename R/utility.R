
#' @title Generating an Rmoleculeset from an SDF file
#' 
#' @description This function uses the ChemmineR package to read an SDF file
#' and converts it into an \code{Rmoleculeset} that can be used as input
#' for the kernel functions \code{sd2gram}, \code{sd2gramSpectrum}, ..., 
#' \code{sd2gram3Dpharma}.
#' 
#' @param sdfFileName The name of the SDF file containing the molecules.
#' 
#' @param detectArom If the molecules in the SDF file have no annotated 
#' aromatic bonds, the ChemmineR function \code{rings} is used for detecting
#' aromaticity. (Default = TRUE).
#' 
#' @param bound Detection of aromaticity can be time consuming if the
#' molecules are large. Detection is only done if the number of atoms is below
#' the given number.  (Default = 70).
#' 
#' @param type Experimental parameter to switch between to types of the
#' function.
#' 
#' @return An instance of Rmoleculeset.
#' 
#' @author Guenter Klambauer <rchemcpp@@bioinf.jku.at>
#' @export

readRmoleculeset <- function(sdfFileName, detectArom=TRUE, bound=70, type=2){
	x <- new(Rmoleculeset)
	sdf <- ChemmineR::read.SDFset(sdfFileName)
	valids <- ChemmineR::validSDF(sdf)
	if (!(all(valids))){
		warnings("Found invalid entries in the sdf file. Removing.")
		sdf <- sdf[which(valids)]
	}
	molNames <- ChemmineR::sdfid(sdf)
	
	if (type==1){
		for (i in 1:length(sdf)){
			mol <- sdf[[i]]				
			aB <- atomblock(mol)
			bB <- bondblock(mol)
			atoms <- sapply(strsplit(rownames(aB),"_"),.subset2,1)
			k <- length(atoms)
			
			if (detectArom & k <= bound){
				ringS <- ChemmineR::rings(mol,arom=TRUE)
				if (is.null(ringS$RINGS)){
					aromaticIdx <- rep(FALSE,k)
				} else {
					aromaticIdx <-  (1:k) %in% as.integer(unlist(sapply(
											ringS$RINGS[which(ringS$AROMATIC)],
											function(x){
												sapply(strsplit(x,"_"),.subset2,2)
											})))
				}
				
				aIdx <- which(aromaticIdx)
				
				bB[which((bB[,1] %in% aIdx) & (bB[,2] %in% aIdx)),3] <- 4	
				
			}
			
			#browser()
			
			m <- new(Rmolecule)
			sapply(atoms,function(kk)  m$addAtom(kk))
			bB[,1:2] <- bB[,1:2]-1
			#bB[which(bB[,3]==4),3] <- 3 
			
			#browser()
			apply(bB,1, function(x) m$linkAtoms(x[1],x[2],x[3]))
			
			
			x$addMoleculeCopy(m)
			
		}
	} else {
		if (detectArom){
			
			for (i in 1:length(sdf)){
				mol <- sdf[[i]]			
				aB <- atomblock(mol)
				bB <- bondblock(mol)
				mode(bB) <- "integer"
				atoms <- sapply(strsplit(rownames(aB),"_"),.subset2,1)
				k <- length(atoms)
				
				if (k <=bound){
					ringS <- ChemmineR::rings(mol,arom=TRUE)
					if (is.null(ringS$RINGS)){
						aromaticIdx <- rep(FALSE,k)
					} else {
						aromaticIdx <-  (1:k) %in% as.integer(unlist(sapply(
												ringS$RINGS[which(ringS$AROMATIC)],
												function(x){
													sapply(strsplit(x,"_"),.subset2,2)
												})))
					}
					
					aIdx <- which(aromaticIdx)
					bB[which((bB[,1] %in% aIdx) & (bB[,2] %in% aIdx)),3] <- 4
					mol@bondblock <- bB
				}
				sdf[[i]] <- mol
			}
			
		}
		
		# make 7 columns for bondblocks
		bb <- ChemmineR::bondblock(sdf)
		bbrep <- lapply(bb, function(b){
					if (ncol(b)>7){
						b <- b[,1:7,drop=FALSE]
					} else if (ncol(b) <7){
						b <- cbind(b,matrix(0,nrow(b),7-ncol(b)))
					} 
					return(as.matrix(b))
					
				})
		ChemmineR::bondblock(sdf) <- bbrep
		
		sdf2strWrite(sdf,"tmp.sdf",cid=TRUE)
		x$addSD("tmp.sdf",TRUE)
		
	}
	
	
	ll <- list(x,sdf,molNames)
	names(ll) <- c("Rmoleculeset","SDFset","moleculeNames")
	return(ll)
	#return(x)
}


#' @title Converting SDFset objects to Rmoleculesets
#' 
#' @param SDFset The SDFset object of ChemmineR. 
#' @param detectArom Flag that determines, whether aromatic structures should 
#' be detected. (Default = TRUE).
#' @param bound Detection of aromaticity can be time consuming if the
#' molecules are large. Detection is only done if the number of atoms is below
#' the given number.  (Default = 70).
#' 
#' @return An instance of Rmoleculeset
#'  
#' @noRd 
#' 
#' @author Guenter Klambauer

SDFsetToRmoleculeset <- function(SDFset, detectArom=TRUE, bound=70,type=2){
	x <- new(Rmoleculeset)
	
	valids <- ChemmineR::validSDF(SDFset)
	if (!(all(valids))){
		warnings("Found invalid entries in the sdf file. Removing.")
		SDFset <- SDFset[which(valids)]
	}
	molNames <- ChemmineR::sdfid(SDFset)
	
	if (type==1){
		
		for (i in 1:length(SDFset)){
			mol <- SDFset[[i]]			
			aB <- atomblock(mol)
			bB <- bondblock(mol)
			mode(bB) <- "integer"
			atoms <- sapply(strsplit(rownames(aB),"_"),.subset2,1)
			k <- length(atoms)
			
			if (detectArom & k<=bound){
				ringS <- ChemmineR::rings(mol,arom=TRUE)
				if (is.null(ringS$RINGS)){
					aromaticIdx <- rep(FALSE,k)
				} else {
					aromaticIdx <-  (1:k) %in% as.integer(unlist(sapply(
											ringS$RINGS[which(ringS$AROMATIC)],
											function(x){
												sapply(strsplit(x,"_"),.subset2,2)
											})))
				}
				
				aIdx <- which(aromaticIdx)
				bB[which((bB[,1] %in% aIdx) & (bB[,2] %in% aIdx)),3] <- 4	
			}
			
			#browser()
			
			m <- new(Rmolecule)
			sapply(atoms,function(kk)  m$addAtom(kk))
			bB[,1:2] <- bB[,1:2]-1
			#bB[which(bB[,3]==4),3] <- 3 
			
			apply(bB,1, function(x) m$linkAtoms(as.integer(x[1]),
								as.integer(x[2]),as.integer(x[3])))
			
			
			x$addMoleculeCopy(m)
			
		}
	} else {
		
		if (detectArom){
			
			for (i in 1:length(SDFset)){
				mol <- SDFset[[i]]			
				aB <- atomblock(mol)
				bB <- bondblock(mol)
				mode(bB) <- "integer"
				atoms <- sapply(strsplit(rownames(aB),"_"),.subset2,1)
				k <- length(atoms)
				
				if (k <= bound){
					ringS <- ChemmineR::rings(mol,arom=TRUE)
					if (is.null(ringS$RINGS)){
						aromaticIdx <- rep(FALSE,k)
					} else {
						aromaticIdx <-  (1:k) %in% as.integer(unlist(sapply(
												ringS$RINGS[which(ringS$AROMATIC)],
												function(x){
													sapply(strsplit(x,"_"),.subset2,2)
												})))
					}
					
					aIdx <- which(aromaticIdx)
					bB[which((bB[,1] %in% aIdx) & (bB[,2] %in% aIdx)),3] <- 4
					mol@bondblock <- bB
				}
				SDFset[[i]] <- mol
			}
			
		}
		
		# make 7 columns for bondblocks
		bb <- ChemmineR::bondblock(SDFset)
		bbrep <- lapply(bb, function(b){
					if (ncol(b)>7){
						b <- b[,1:7,drop=FALSE]
					} else if (ncol(b) <7){
						b <- cbind(b,matrix(0,nrow(b),7-ncol(b)))
					} 
					return(as.matrix(b))
					
				})
		ChemmineR::bondblock(SDFset) <- bbrep
		
		
		sdf2strWrite(SDFset,"tmp.sdf",cid=TRUE)
		x$addSD("tmp.sdf",TRUE)
		
	}
	
	return(list(x,SDFset,molNames))
}



#' @title createRMolecule
#' 
#' @description Creates an \"Rmolecule\" from an atom-vector and a bond-matrix 
#' 
#' @usage createRMolecule(atoms, bonds)
#' 
#' @param atoms A vector containing the symbol names of all atoms in the molecule
#' @param bonds A matrix with the same number of rows and columns as the atoms-vector
#' containing the type of bonds between the atoms
#' @return an instance of "molecule"
#' @author Michael Mahr <rchemcpp@@bioinf.jku.at>
#' 
#' @examples 
#' m <- createRMolecule(c("C","C"),matrix(c(0,3,3,0),nrow=2))
#' @export


createRMolecule = function(atoms, bonds)
{
	if(!is.vector(atoms)) stop("atoms must be vector")
	if(!is.matrix(bonds)) stop("bonds must be matrix")
	
	bonds = apply(bonds,1,as.integer)
	
	if (length(atoms) != nrow(bonds)) stop("matrix must have as many rows as vector")
	if (length(atoms) != ncol(bonds)) stop("matrix must have as many cols as vector has rows")
	
	if (! isSymmetric(bonds) ) stop("matrix must be symmetric")
	if (! all(diag(bonds) == 0)) stop("cannot connect atoms with themselves")	
	
	if ( (max(bonds) > 3) || (min(bonds) < 0 ) ) stop("Bond types can only range from 1 to 3") 
	
	m = new(Rchemcpp::Rmolecule)
	
	for (k in atoms)
	{
		m$addAtom(k)
	}
	
	for (i in 1:length(atoms)){
		for (j in 1:length(atoms)){
			if ((i < j) && (bonds[i,j] != 0)){
				m$linkAtoms(i-1,j-1,bonds[i,j]);
			}
		}
	}
	
	return (m)
}



#.makeValidSDF <- function(inFile,outFile){
#	library(ChemmineR)
#	x <- ChemmineR::read.SDFset(inFile)
#	idx <- ChemmineR::validSDF(x)
#	x <- x[idx]
#	
#	# make 6 columns for bondblocks
#	bb <- ChemmineR::bondblock(x)
#	bbrep <- lapply(bb, function(b){
#				if (ncol(b)>6){
#					b <- b[,1:6,drop=FALSE]
#				} else if (ncol(b) <6){
#					b <- cbind(b,matrix(0,nrow(b),6-ncol(b)))
#				} 
#				return(as.matrix(b))
#				
#			})
#	ChemmineR::bondblock(x) <- bbrep
#	ChemmineR::write.SDF(x,outFile,cid=TRUE)
#}
## Convert SDF to SDFstr Class for Export to File
## Function allows to customize output via optional arguments 

sdf2strWrite <- function(sdf, file, cid=NULL, ...) {
	## Checks
	if(class(sdf)!="SDFset") stop("Function expects molecule object of class SDF as input!")	
	if(cid==TRUE) {	sdflist <- lapply(cid(sdf), function(x) sdf2str(sdf=sdf[[x]], cid=x, ...)) }	
	if(cid==FALSE) { sdflist <- lapply(cid(sdf), function(x) sdf2str(sdf=sdf[[x]], ...)) } 
	cat(unlist(sdflist), sep="\n", file=file)
	
}

sdf2str <- function (sdf, file, head, ab, bb, db, cid=NULL, sig=FALSE, ...){
	
	if(missing(head)) {
		head <- as.character(sdf[[1]])
		if(sig==TRUE) head[2] <- paste("ChemmineR-", format(Sys.time(), "%m%d%y%H%M"), "XD", sep="")	
		if(length(cid)==1) head[1] <- cid
	}
	
	## Atom block
	if(missing(ab)) {
		ab <- sdf[[2]]
		ab <- cbind(Indent="", format(ab[,1:3], width=9, justify="right"), 
				A=format(gsub("_.*", "", rownames(ab)), width=1, justify="left"), 
				Space="", format(ab[,-c(1:3)], width=2, justify="right"))
		ab <- sapply(seq(along=ab[,1]), function(x) paste(ab[x, ], collapse=" "))
	}
	
	## Bond block
	if(missing(bb)) {
		bb <- sdf[[3]]
		bb <- cbind(Indent="", format(bb, width=3, justify="right"))
		bb <- sapply(seq(along=bb[,1]), function(x) paste(bb[x, ], collapse=""))
	}
	
	## Data block
	if(missing(db)) {
		db <- sdf[[4]]
		if(length(db)>0) {
			dbnames <- paste("> <", names(db), ">", sep="")
			dbvalues <- as.character(db)
			db <- as.vector(rbind(dbnames, dbvalues, ""))
		} else {
			db <- NULL
		}
	}
	
	## Assemble in character vector
	sdfstrvec <- c(head, ab, bb, "M  END", db, "$$$$")
	return(sdfstrvec)
	
}



##########################
## (7) Plotting Methods ##
##########################
## Plot single CMP Structure
plotStruc <- function(sdf, atomcex=1.2, atomnum=FALSE, no_print_atoms=c("C"), noHbonds=TRUE, bondspacer=0.12, colbonds=NULL, bondcol="red", ...) {
	toplot <- list(atomblock=cbind(atomblock(sdf)[,c(1:2)], as.matrix(bonds(sdf, type="bonds")[,-1])), bondblock=cbind(as.matrix(as.data.frame(bondblock(sdf))[,1:3]), bondcol=1))
	## Add bond color
	toplot[[2]][, "bondcol"] <- toplot[[2]][,"bondcol"] + as.numeric((toplot[[2]][,"C1"] %in% colbonds) & (toplot[[2]][,"C2"] %in% colbonds))
	## Create empty plot with proper dimensions
	plot(toplot[[1]], type="n", axes=F, xlab="", ylab="", ...)
	## Remove C-hydrogens including their bonds 
	if(noHbonds==TRUE) {
		nonbonded <- !1:length(toplot[[1]][,1]) %in% sort(unique(as.vector(toplot[[2]][,1:2])))
		nonbonded <- as.data.frame(toplot[[1]])[nonbonded,]
		CHbondindex <- sapply(seq(toplot[[2]][,1]), function(x) paste(sort(gsub("_.*", "", rownames(toplot[[1]]))[toplot[[2]][x,1:2]]), collapse="") == "CH")
		toplot[[1]] <- toplot[[1]][sort(unique(as.numeric(toplot[[2]][!CHbondindex,1:2]))), ]
		toplot[[2]] <- as.matrix(as.data.frame(toplot[[2]])[!CHbondindex,]) 
		toplot[[1]] <- as.matrix(rbind(toplot[[1]], nonbonded))
	}	
	## Plot bonds
	z <- toplot[[2]][, "bondcol"]; z[z==2] <- bondcol # Stores bond coloring data
	for(i in seq(along=toplot[[2]][,1])) {
		x <- toplot[[1]][gsub("*.*_", "", rownames(toplot[[1]])) %in% toplot[[2]][i,1:2],1]
		y <- toplot[[1]][gsub("*.*_", "", rownames(toplot[[1]])) %in% toplot[[2]][i,1:2],2]
		## Plot single bonds
		if(toplot[[2]][i,3]==1) {
			lines(x=x, y=y, lty=1, lwd=3, col=z[i]) 
		} 
		## Plot double bonds
		if(toplot[[2]][i,3]==2) {
			rslope <- (atan(diff(y)/diff(x))*180/pi)/90
			lines(x=x-rslope*bondspacer, y=y+(1-abs(rslope))*bondspacer, lty=1, lwd=3, col=z[i])
			lines(x=x+rslope*bondspacer, y=y-(1-abs(rslope))*bondspacer, lty=1, lwd=3, col=z[i])
		}
		## Plot triple bonds
		if(toplot[[2]][i,3]==3) {
			rslope <- (atan(diff(y)/diff(x))*180/pi)/90
			bondspacer <- bondspacer * 2
			lines(x=x, y=y, lty=1, lwd=3, col=z[i]) 
			lines(x=x-rslope*bondspacer, y=y+(1-abs(rslope))*bondspacer, lty=1, lwd=3, col=z[i])
			lines(x=x+rslope*bondspacer, y=y-(1-abs(rslope))*bondspacer, lty=1, lwd=3, col=z[i])
			
		}
		
		## Plot aromatic or any other bonds
		if(toplot[[2]][i,3] > 3 | toplot[[2]][i,3] < 1 ) {
			lines(x=x, y=y, lty=2, lwd=3, col=z[i]) 
		}
	}
	## Exclude certain atoms from being printed
	exclude <- paste("(^", no_print_atoms, "_)", sep="", collapse="|")
	labelMA <- toplot[[1]][!grepl(exclude, rownames(toplot[[1]])), , drop=FALSE] # Added July 31, 2012: 'drop=FALSE'
	## Add charges 
	charge <- c("0"="", "3"="3+", "2"="2+", "1"="+", "-1"="-", "-2"="2-", "-3"="3-")
	charge <- charge[as.character(labelMA[,"charge"])]
	## Add hydrogens to non-charged/non-C atoms according to valence rules (some SD files require this)
	Nhydrogens <- c("0"="", "1"="H", "2"="H2", "3"="H3", "4"="H4", "5"="H5", "6"="H6", "7"="H7", "8"="H8") 
	hydrogens <- (labelMA[, "Nbondrule"] + labelMA[, "charge"]) - labelMA[,"Nbondcount"]
	hydrogens[labelMA[,"charge"]!=0] <- 0; hydrogens[hydrogens < 0] <- 0
	hydrogens <- Nhydrogens[as.character(hydrogens)]
	## Plot data
	if(is.vector(labelMA)) labelMA <- matrix(labelMA, 1, 2, byrow=TRUE, dimnames=list(rownames(toplot[[1]])[!grepl(exclude, rownames(toplot[[1]]))], c("C1", "C2")))
	if(is.matrix(labelMA) & length(labelMA[,1])>=1) {
		atomcol <- gsub("_.*", "", rownames(labelMA)); atomcol[!grepl("N|C|O|H", atomcol)] <- "any"; mycol <- c(C="black", H="black", N="blue", O="red", any="green"); atomcol <- mycol[atomcol]
		
		## Overplot nodes to display atom labels
		points(x=labelMA[,1], y=labelMA[,2], col="white", pch=16, cex=2.8)
		## Plot atom labels
		if(atomnum==TRUE) {
			text(x=labelMA[,1], y=labelMA[,2], paste(gsub("_", "", rownames(labelMA)), hydrogens, charge, sep=""), cex=atomcex, col=atomcol) 
		} else {
			text(x=labelMA[,1], y=labelMA[,2], paste(gsub("_.*", "", rownames(labelMA)), hydrogens, charge, sep=""), cex=atomcex, col=atomcol)
		}
	}
}
## Usage:
# plotStruc(sdf=sdfset[[2]], atomcex=1.2, atomnum=F, no_print_atoms=c("C"), noHbonds=TRUE, bondspacer=0.08)
# par(mfrow=c(2,3)); for(i in 1:6) plotStruc(sdf=sdfset[[i]], atomcex=1.8, atomnum=F, no_print_atoms=c("C"), noHbonds=TRUE, bondspacer=0.08)

## Plot method for single SDF object
setMethod(f="plot", signature="SDF", definition=function(x, print=TRUE, ...) { 
			plotStruc(sdf=x, ...)
			if(print==TRUE) { return(x) } 
		}
)

## Plot method for multiple SDF objects in SDFset
setMethod(f="plot", signature="SDFset",
		definition=function(x, griddim, print_cid=cid(x), print=TRUE, ...) {
			if(missing(griddim)) {
				mydim <- ceiling(sqrt(length(x)))
				griddim <- c(mydim, mydim)
			}
			par(mfrow=griddim)
			for(i in 1:length(x)) { plotStruc(sdf=x[[i]], main=print_cid[i], ...) }
			if(print==TRUE) { return(SDFset2SDF(x)) }
		})

