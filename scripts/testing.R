library(Rchemcpp)
library(apcluster)

sdfolder <- system.file("sample_data",package="Rchemcpp")
sdf <- list.files(sdfolder,full.names=TRUE,pattern="small")
K <- sd2gram(sdf,returnNormalized=TRUE)

depth <- 10
K <- sd2gramSpectrum(sdf,kernelType="spectrum",depthMax=as.integer(depth),returnNormalized=TRUE)
K <- sd2gramSpectrum(sdf,kernelType="tanimoto",depthMax=as.integer(depth),returnNormalized=TRUE)
K <- sd2gramSpectrum(sdf,kernelType="minmaxTanimoto",depthMax=as.integer(depth),returnNormalized=TRUE)
K <- sd2gramSpectrum(sdf,kernelType="marginalized",depthMax=as.integer(depth),returnNormalized=TRUE)
K1 <- sd2gramSpectrum(sdf,kernelType="lambda",depthMax=as.integer(5),lambdaKernelLambda=4,returnNormalized=TRUE)
K2 <- sd2gramSpectrum(sdf,kernelType="lambda",depthMax=as.integer(depth),lambdaKernelLambda=4,returnNormalized=TRUE)
all(K1==K2)

K1 <- sd2gramSubtree(sdf,kernelType="sizebased",depthMax=as.integer(10),branchKernelUntilN=TRUE)
K2 <- sd2gramSubtree(sdf,kernelType="branchingbased",depthMax=as.integer(10),branchKernelUntilN=TRUE)
all(K1==K2)


K <- sd2gram3Dspectrum(sdf,returnNormalized=TRUE)

K <- sd2gram3Dpharma(sdf,returnNormalized=TRUE)



r <- apcluster(K)
plot(r,K)


#getMoleculeNamesFromSDF(sdf)
getMoleculePropertyFromSDF(sdf,property="Activity")


#http://www.clab.kwansei.ac.jp/mining/datasets/PAKDD2000/okd.htm
#[Debnath 91] A. K. Debnath, R. L. Lopez de Compadre, G. Debnath, A. J. Shusterman & C. Hansch: Structure-Activity Relationship of Mutagenic Aromatic and Heteroaromatic Nitro Compounds. Correlation with Molecular Orbital Energies and Hydrophobicity, Journal of Medicinal Chemistry, Vol.34, pp.786-797 (1991).


#perl one-liner to remove non-ASCII Characters: perl -i.bk -pe 's/[^[:ascii:]]//g;' file

## svm prediction ################################################################
library(kernlab)
library(Rchemcpp)

sdfolder <- system.file("sample_data",package="Rchemcpp")
sdf <- list.files(sdfolder,full.names=TRUE,pattern="mutag.sdf")
K <- sd2gramSpectrum(sdf,kernelType="lambda",depthMax=as.integer(15),lambdaKernelLambda=4,returnNormalized=TRUE)
y <- as.numeric(getMoleculePropertyFromSDF(sdf,"Activity"))
idxRm <- which(y==-99)
K <- K[-idxRm,-idxRm]
y <- y[-idxRm]


prediction <- vector("numeric",length(y))

for (idx in 1:length(y)){
	trainK <- as.kernelMatrix(K)[-idx,-idx]
	model <- ksvm(trainK,y=y[-idx],kernel="matrix",type="nu-svr")
	testK <- as.kernelMatrix(K[idx,-idx,drop=FALSE][ ,SVindex(model),drop=FALSE])
	prediction[idx] <- as.numeric(predict(model, testK))
}

plot(prediction,y)
cvError <- y - prediction


################################################################################
library(Rchemcpp)
sdfolder <- system.file("sample_data",package="Rchemcpp")
sdf <- list.files(sdfolder,full.names=TRUE,pattern="small")

sdf <- read.SDFset(sdf)
sdf2 <- sdf[1]

#source("/home/klambaue/Dropbox/WorkspaceSequencing/RchemcppRep/Rchemcpp/R/sd2gramSpectrum.R")
#source("/home/klambaue/Dropbox/WorkspaceSequencing/RchemcppRep/Rchemcpp/R/utility.R")
#source("/home/klambaue/Dropbox/WorkspaceSequencing/RchemcppRep/Rchemcpp/R/utility.R")
#source("/home/klambaue/Dropbox/WorkspaceSequencing/RchemcppRep/Rchemcpp/R/sd2gramSubtree.R")


# input: sdf filenames
library(Rchemcpp)
sdfolder <- system.file("sample_data",package="Rchemcpp")


sdf <- list.files(sdfolder,full.names=TRUE,pattern="small")
#sdf <- "~/Desktop/test.sdf"
sd2gram(sdf) #works
sd2gramSubtree(sdf) #works
sd2gramSpectrum(sdf) #segfault


sdf <- list.files(sdfolder,full.names=TRUE,pattern="tiny")
sd2gram3Dpharma(sdf) #works
sd2gram3Dspectrum(sdf) #works


sdf <- list.files(sdfolder,full.names=TRUE,pattern="tiny")
sdf2 <- list.files(sdfolder,full.names=TRUE,pattern="small")
sd2gram(sdf,sdf2) #setComparisonSet does not work!!
sd2gramSpectrum(sdf,sdf2) #segfault
sd2gramSubtree(sdf,sdf2)  #works


sdf <- list.files(sdfolder,full.names=TRUE,pattern="tiny")
sdf2 <- list.files(sdfolder,full.names=TRUE,pattern="small")
sd2gram3Dpharma(sdf,sdf2) #works 
sd2gram3Dspectrum(sdf,sdf2) #segfault


# input: SDFsets
library(Rchemcpp)
sdfolder <- system.file("sample_data",package="Rchemcpp")
sdf <- list.files(sdfolder,full.names=TRUE,pattern="small")
#sdf <- "/media/Daten_/BioinfData/Chemoinformatics/Mohr/ci900367j_si_001/compoundWCoordinates.sdf"
sdf <- read.SDFset(sdf)

sd2gram(sdf) #works
sd2gramSpectrum(sdf) #segfault
sd2gramSubtree(sdf) #works
sd2gram3Dpharma(sdf) #works
sd2gram3Dspectrum(sdf) #segfault

sdf2 <- read.SDFset(list.files(sdfolder,full.names=TRUE,pattern="tiny"))

sd2gram(sdf,sdf2) #setComparisonSet does not work!!gives back 5x5!!
sd2gramSpectrum(sdf,sdf2) #segfault
sd2gramSubtree(sdf,sdf2)  #works
sd2gram3Dpharma(sdf[1:2],sdf2[1:3]) #works
sd2gram3Dspectrum(sdf[1:5],sdf2[1:4]) #segfault



# input: Rmoleculesets
library(Rchemcpp)
sdfolder <- system.file("sample_data",package="Rchemcpp")
sdf <- list.files(sdfolder,full.names=TRUE,pattern="small")
sdf <- Rchemcpp:::readRmoleculeset(sdf)
sdf <- sdf[[1]]

sd2gram(sdf) #works
sd2gramSpectrum(sdf) #segfault!!
sd2gramSubtree(sdf) #works
sd2gram3Dpharma(sdf) #works
sd2gram3Dspectrum(sdf) #segfault

sdf2 <- Rchemcpp:::readRmoleculeset(list.files(sdfolder,full.names=TRUE,pattern="tiny"))
sdf2 <- sdf2[[1]]

sd2gram(sdf,sdf2) #setComparisonSet does not work!! gives back 5x5!!
sd2gramSpectrum(sdf,sdf2) #segfault
sd2gramSubtree(sdf,sdf2)  #works
sd2gram3Dpharma(sdf,sdf2) #works but gives back 5x5!!
sd2gram3Dspectrum(sdf,sdf2) #segfault

###
library(Rchemcpp)
sdfolder <- system.file("sample_data",package="Rchemcpp")
sdf <- list.files(sdfolder,full.names=TRUE,pattern="small")
sdf <- read.SDFset(sdf)
sd2gramSpectrum(sdf[1:5])
sd2gramSpectrum(sdf[6:7])
sd2gramSpectrum(sdf[7:8])
sd2gramSpectrum(sdf[8:9])
sd2gramSpectrum(sdf[10:12])
for (i in 1:100){
	cat("iteration: ",i,"\n")
	j <- round(runif(n=1,min=1,max=30))
	k <- round(runif(n=1,min=1,max=30))
	if (j!=k)
		sd2gramSpectrum(sdf[c(j,k)])
	
}

##############################################################################

library(Rchemcpp)
sdfolder <- system.file("sample_data",package="Rchemcpp")
sdf <- list.files(sdfolder,full.names=TRUE,pattern="small")
K <- sd2gram(sdf,returnNormalized=TRUE,
		moleculeNameProperty="Compound Name")
sdf2 <- read.SDFset(sdf)
K2 <- sd2gram(sdf2,detectArom=FALSE)


sdfolder <- system.file("sample_data",package="Rchemcpp")
sdf <- list.files(sdfolder,full.names=TRUE,pattern="mutag.sdf")
K <- sd2gram(sdf,returnNormalized=TRUE,
		moleculeNameProperty="Compound Name")
sdf2 <- read.SDFset(sdf)
K2 <- sd2gram(sdf2,detectArom=FALSE)


