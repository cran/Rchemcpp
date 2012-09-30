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
