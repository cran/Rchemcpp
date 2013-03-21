
#' @title sd2gram3Dspectrum - Similarity of molecules by fast 
#' approximations of the pharmacophore kernel
#' 
#' @description This tool implements the six discrete approximations of the pharmacophore
#' kernel presented in "The pharmacophore kernel for virtual screening 
#' with support vector machines" (\cite{Mahe, 2006}).
#' 
#' @param sdf File containing the molecules. Must be in MDL file format
#' (MOL and SDF files). For more information on the file format see 
#' http://en.wikipedia.org/wiki/Chemical_table_file. Default = "missing".
#' @param sdf2 A second file containing molecules. Must also be in SDF.
#' If specified the molecules of the first file will be compared with the 
#' molecules of this second file. Default = "missing".
#' @param chargesFileName A character with the name of the file containing
#' the atom charges. Default = missing.
#' @param chargesFileName2 A character with the name of the file containing
#' the atom charges. Default = missing.
#' @param kernelType Type of kernel to be used. Possible choices are 
#' 3-points spectrum kernel ("3Pspectrum"), 
#' 3-points binary kernel ("3Pbinary"),
#' 3-points Tanimoto kernel ("3Ptanimoto"),
#' 2-points spectrum kernel ("2Pspectrum"),
#' 2-points binary kernel ("2Pbinary"),
#' 2-points Tanimoto kernel ("2Ptanimoto"). Default = "3Pspectrum".		
#' @param depthMax The maximal length of the molecular fragments. Default = 3.
#' @param nBins number of bins used to discretize the inter-atomic lengths. 
#' An adequate value for the number of bins is between 20 and 30. 
#' Default = 20.
#' @param distMin minimum distance for inter-atomic distance range. 
#' Default = 0.
#' @param distMax maximum distance in angstrom for inter-atomic distance range. 
#' Default = 20. 
#' @param chargesThreshold specifies a threshold above which partial charges
#' are considered as positive/negative.
#' By default this threshold is zero, and every positive (resp. negative) 
#' partial charge is seen as a positive (resp. negative) charge. 
#' However, it might be interesting to consider a threshold of 0.2 for
#' example, in which case only partial charges greater than 0.2 
#' (resp. smaller than -0.2) would be seen as positive (resp. negative). 
#' Default = 0. 
#' @param flagRemoveH A logical that indicates whether H-atoms should be 
#' removed or not.
#' @param morganOrder The order of the DeMorgan Indices to be used. If set to
#' zero no DeMorgan Indices are used. The higher the order the more different
#' types of atoms exist and consequently the more dissimilar will be the molecules.
#' @param silentMode Whether or not the program should print progress reports
#' to the standart output.
#' @param returnNormalized A logical specifying whether a normalized kernel
#' matrix should be returned. Default = TRUE.
#' @param detectArom Whether aromatic rings should be detected and aromatic
#' bonds should a special bond type. (Default = TRUE).
#' @examples 
#' sdfolder <- system.file("sample_data",package="Rchemcpp")
#' sdf <- list.files(sdfolder,full.names=TRUE,pattern="small")
#' K <- sd2gram3Dspectrum(sdf)
#' @return A numeric matrix containing the similarity values between the
#' molecules.
#' @author Michael Mahr <rchemcpp@@bioinf.jku.at>
#' c++ function written by Jean-Luc Perret and Pierre Mahe
#' @references (Mahe, 2006) --  P. Mahe, L. Ralaivola, V. Stoven, and J.-P. Vert.
#' The pharmacophore kernel for virtual screening
#' with support vector machines. Technical Report, HAL:ccsd-00020066, Ecole des
#' Mines de Paris, March 2006.
#' 
#'  
#' 
#' @export


sd2gram3Dspectrum = function(sdf, sdf2, 
		chargesFileName = "", chargesFileName2 = "", 
		kernelType = c("3Pspectrum", "3Pbinary", "3Ptanimoto", 
				"2Pspectrum", "2Pbinary", "2Ptanimoto" ),
		depthMax = as.integer(3), nBins = as.integer(20), distMin = 0, 
		distMax = 20, flagRemoveH = FALSE, morganOrder = as.integer(0), 
		chargesThreshold = 0, silentMode = FALSE, returnNormalized = TRUE, detectArom=TRUE)
{
	
	if(!is.character(chargesFileName)) stop("chargesFileName must be string")
	if(!is.character(chargesFileName2)) stop("chargesFileName2 must be string")
	
	if(!is.character(kernelType)) stop("kernelType must be a string")
	
	if(!is.integer(depthMax)) stop("depthMax must be integer")
	if(!is.logical(flagRemoveH)) stop("flagRemoveH must be logical")
	if(!is.numeric(morganOrder)) stop("morganOrder must be integer")
	morganOrder <- as.integer(morganOrder)
	if(!is.numeric(chargesThreshold)) stop("chargesThreshold must be numeric")
	
	if(!is.numeric(nBins)) stop("nBins must be integer")
	nBins <- as.integer(nBins)
	if(!is.numeric(distMin)) stop("distMin must be numeric")
	if(!is.numeric(distMax)) stop("distMax must be numeric")
	
	if(!is.logical(silentMode)) stop("silentMode must be logical")
	if(!is.logical(returnNormalized)) stop("returnNormalized must be logical")
	
	if((missing(sdf2)) && (chargesFileName2 != "")) print("chargesFilename2 is useless") 
	
	if((missing(sdf2)) && (chargesFileName == "") && (chargesFileName2 != ""))
		stop("both chargesFileNames have to be specified") 
	if((missing(sdf2)) && (chargesFileName != "") && (chargesFileName2 == "")) 
		stop("both chargesFileNames have to be specified") 
	
	
	if (missing(kernelType)) {kernelType = kernelType[1]}
	kernelType = match.arg(kernelType)
	
	#for passing
	kernelTypeIndex = as.integer(which(c("3Pspectrum", "3Pbinary", 
							"3Ptanimoto", "2Pspectrum", "2Pbinary", 
							"2Ptanimoto" ) == kernelType)-1)
	
	
	if(inherits(sdf,"SDFset")){
		aSet <- SDFsetToRmoleculeset(sdf,detectArom=detectArom)[[1]]
		molnames <- ChemmineR::sdfid(sdf)
		molnames2 <- molnames
		if (!missing(sdf2)){
			if (inherits(sdf2,"SDFset")){
				aSet2 <- SDFsetToRmoleculeset(sdf2,detectArom=detectArom)[[1]]
				molnames2 <- ChemmineR::sdfid(sdf2)
			} else 
				stop("Input must be existing SDF files or \"SDFset\" objects.")	
		}
		
	} else if (inherits(sdf,"Rcpp_Rmoleculeset")) {
		aSet <- sdf
		molnames <- NULL
		molnames2 <- NULL
		if (!missing(sdf2)){
			if (inherits(sdf2,"Rcpp_Rmoleculeset")){
				aSet2 <- sdf2
			} else 
				stop("Input must be existing SDF files or \"SDFset\" objects.")
		}
		
	} else if (is.character(sdf) & file.exists(sdf)) {
		
		TT <- try({
					aSet <- new(Rmoleculeset)					
					aSet$addSD(sdf,TRUE)
					molnames <- getMoleculeNamesFromSDF(sdf)
					molnames2 <- molnames
					
					if (!missing(sdf2)){
						if (is.character(sdf2) & file.exists(sdf2)){
							aSet2 <- new(Rmoleculeset)
							aSet$addSD(sdf2,TRUE)
							molnames2 <- getMoleculeNamesFromSDF(sdf2)		
						} else 
							stop("Input must be existing SDF files or \"SDFset\" objects.")	
					} 
				})
		if (inherits(TT,"try-error")){
			aSetList <- readRmoleculeset(sdf,detectArom=detectArom)
			aSet <- aSetList[[1]]
			molnames <- aSetList[[3]]
			molnames2 <- aSetList[[3]]
			if (!missing(sdf2)){
				if (is.character(sdf2) & file.exists(sdf2)){
					aSetList2  <-  readRmoleculeset(sdf2,detectArom=detectArom)
					aSet2 <- aSetList2[[1]]
					molnames2 <- aSetList2[[3]]
				} else 
					stop("Input must be existing SDF files or \"SDFset\" objects.")	
			} 
		}
		
	} else {
		stop("Input must be existing SDF files or \"SDFset\" objects.")
	}
	
	
	
	# if sflag = 1 --> SEPARATE TEST SET
	if( !missing(sdf2)){
		
	
		if( !silentMode ){
			print(paste("*** sd file added ; number of molecules in the set 1 =  ",
							aSet$numMolecules() ));
			print(paste("*** sd file added ; number of molecules in the set 2 =  ",
							aSet2$numMolecules() ));
		}
		
		
		# read partial charges if specified on command line
		if(chargesFileName != ""){
			if( !silentMode ){
				print("reading partial charges");
			}
			aSet$readPartialCharges(chargesFileName);
			aSet2$readPartialCharges(chargesFileName2);
		}
		
		
		# remove H when specified on command line
		if( flagRemoveH == TRUE ){
			if( !silentMode ){
				print("removing hydrogens");
			}
			aSet$hideHydrogens();
			aSet2$hideHydrogens();
		}
		
		# compute Morgan labels
		if( !silentMode ){
			print(paste("setting morgan labels ", morganOrder));
		}
		aSet$setMorganLabels( morganOrder );
		aSet2$setMorganLabels( morganOrder );
		
		
		# introduce partial charges in Morgan labels if specified on command line
		if(chargesFileName != ""){
			if( !silentMode ){
				print("setting partial charges");
			}
			aSet$setMorganChargesLabels(chargesThreshold);
			aSet2$setMorganChargesLabels(chargesThreshold);
		}
		
		
		# set kashima probabilities if kernelType = marginalized
		#if(kernelType == "marginalized"){
		#  aSet$setKashimaKernelParam( kernelParam, convgce, skipSkeleton );
		#  aSet2$setKashimaKernelParam( kernelParam, convgce, skipSkeleton );
		#}
		
		
		# initialize gram matrices
		aSet$initializeSelfKernel( 0.0);
		aSet2$initializeSelfKernel(0.0);
		aSet$setComparisonSetCopy( aSet2 );
		aSet$initializeGram( 0.0 );
		if( !silentMode){
			print("#### initialization Gram OK");
		}
		
		
		# 3 - compute the gram matrix
		# ----------------------------
		gramSpectrum3D_test( aSet, depthMax, kernelTypeIndex, nBins, distMin,
				distMax, silentMode);
		if( !silentMode ){
			print("gramComputeSpectrum (test) OK");
		}
		
		if (kernelType == "3Pspectrum"){
			aSet$normalizeGram();
		}
		else if (kernelType == "3Pbinary"){
			aSet$normalizeGram();
		}
		else if (kernelType == "3Ptanimoto"){
			aSet$normalizeTanimoto();
		}
		else if (kernelType == "2Pspectrum"){
			aSet$normalizeGram();
		}
		else if (kernelType == "2Pbinary"){
			aSet$normalizeGram();
		}
		else if (kernelType == "2Ptanimoto"){
			aSet$normalizeTanimoto();
		}
		
		
		if( !silentMode ){
			print("normalize gram (test) OK");
		}
		
		if (returnNormalized == FALSE)
		{	
			K <- do.call(rbind,aSet$getGram() )
			
		}
		else
		{		
			K <- do.call(rbind,aSet$getGramNormal() )
			
		}
		
		#aSet$writeSelfKernelList( outputDir + baseName + "_test", silentMode );
		#aSet$getComparisonSet$writeSelfKernelList( outputDir + baseName + "_train", silentMode ); !!!
		
	}else{ # --> SELF KERNEL

		
		# read partial charges if specified on command line
		if(chargesFileName != ""){
			if( !silentMode ){
				print("reading partial charges");
			}
			aSet$readPartialCharges(chargesFileName);
		}
		
		
		# remove H when specified on command line
		if( flagRemoveH == TRUE ){
			if( !silentMode ){
				print("removing hydrogens");
			}
			aSet$hideHydrogens();
		}
		
		# compute Morgan labels
		if( !silentMode ){
			print(paste("setting morgan labels ", morganOrder));
		}
		aSet$setMorganLabels( morganOrder );
		
		
		# introduce partial charges in Morgan labels if specified on command line
		if(chargesFileName != ""){
			if( !silentMode ){
				print("setting partial charges");
			}
			aSet$setMorganChargesLabels(chargesThreshold);
		}
		
		
		# set kashima probabilities if kernelType = marginalized
		#if(kernelType == "marginalized"){
		#  aSet$setKashimaKernelParam( kernelParam, convgce, skipSkeleton );
		#}
		
		
		# initialize gram matrices
		aSet$setComparisonSetSelf();
		aSet$initializeGram( 0.0 );
		aSet$initializeSelfKernel( 0.0);
		if( !silentMode)
		{
			print("#### initialization Gram OK");
		}
		
		
		# 3 - compute the gram matrix
		# ---------------------------
		gramSpectrum3D_self( aSet, depthMax, kernelTypeIndex, nBins, distMin, distMax, silentMode);
		if( !silentMode ){
			print("gramComputeSpectrum (self) OK");
		}
		
		#normalize Gram
		if (kernelType == "3Pspectrum"){
			aSet$normalizeGram();
		}
		else if (kernelType == "3Pbinary"){
			aSet$normalizeGram();
		}
		else if (kernelType == "3Ptanimoto"){
			aSet$normalizeTanimoto();
		}
		else if (kernelType == "2Pspectrum"){
			aSet$normalizeGram();
		}
		else if (kernelType == "2Pbinary"){
			aSet$normalizeGram();
		}
		else if (kernelType == "2Ptanimoto"){
			aSet$normalizeTanimoto();
		}
		
		
		if( !silentMode ){
			print("normalize gram (self) OK");
		}
		
		
		if (returnNormalized == FALSE)
		{	
			K <- do.call(rbind,aSet$getGram() )
			
		}
		else
		{		
			K <- do.call(rbind,aSet$getGramNormal() )
			
		}
		
		
		#aSet$writeSelfKernelList( outputDir + baseName + "_train", false);
		
	}  
	
	
	
	rownames(K) <- molnames
	colnames(K) <- molnames2
	
	
	return ( K )
	
	
}

