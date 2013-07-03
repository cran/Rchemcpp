#' @title sd2gramSpectrum - Similarity of molecules by walk-based graph kernels
#' 
#' @description This function computes several walk-based graph kernel functions
#' based on finite length walks and a fast implementation for input sd file(s).
#'
#' @usage sd2gramSpectrum(sdFileName, sdFileName2 = "", 
#'		kernelType = c("spectrum", "tanimoto", "minmaxTanimoto","marginalized","lambda"), 
#'		margKernelEndProbability = 0.1, lambdaKernelLambda = 1.0, 
#'		depthMax = as.integer(0), onlyDepthMax = FALSE , flagRemoveH = FALSE, 
#'		morganOrder = as.integer(0), fileType = c("sd","genericsd","kcf"), 
#'		silentMode = FALSE, returnNormalized = FALSE, moleculeNameProperty = "",
#'		moleculeNameProperty2 = "")
#' 
#' @param sdFileName File containing the molecules. Must be in MDL file format
#' (MOL and SDF files). For more information on the file format see 
#' http://en.wikipedia.org/wiki/Chemical_table_file. Default = "missing"
#' @param sdFileName2 A second file containing molecules. Must also be in SDF.
#' If specified the molecules of the first file will be compared with the 
#' molecules of this second file. Default = "missing".
#' @param kernelType Sets which kernel is to be used. Options are "spectrum (Spectrum 
#' kernel) , "tanimoto" (Tanimoto kernel), "minmaxTanimoto" (MinMax Tanimoto kernel),
#' "marginalized (Marginalized kernel approximation) and 
#' "lambda" (LambdaK kernel). Default = "spectrum".
#' @param margKernelEndProbability The ending probability for the marginalized
#' kernel. Default = 0.1. 
#' @param lambdaKernelLambda The lambda parameter of the LambdaK kernel. 
#' Default = 1.0.
#' @param depthMax The maximal length of the molecular fragments.
#' @param onlyDepthMax Whether fragments up to the given length should be 
#' used or only fragments of the given length.
#' @param flagRemoveH A logical that indicates whether H-atoms should be 
#' removed or not.
#' @param morganOrder The order of the DeMorgan Indices to be used. If set to
#' zero no DeMorgan Indices are used. The higher the order the more different
#' types of atoms exist and consequently the more dissimilar will be the molecules.
#' @param fileType Which filetype was submitted.
#' @param silentMode Whether or not the program should print progress reports
#' to the standart output.
#' @param returnNormalized A logical specifying whether a normalized kernel
#' matrix should be returned. Default = TRUE.
#' @param moleculeNameProperty A string which specifies the name of the property
#' of the molecules in the sdFile from which the row names (and column names)
#' are read from. Default = "".
#' @param moleculeNameProperty2 A string which specifies the name of the property
#' of the molecules in the sdFile 2 from which the column names are read from.
#' Default = "".
#' @examples 
#' sdfolder <- system.file("sample_data",package="Rchemcpp")
#' sdf <- list.files(sdfolder,full.names=TRUE,pattern="small")
#' K <- sd2gramSpectrum(sdf, moleculeNameProperty="Compound Name")
#' @return A numeric matrix containing the similarity values between the
#' molecules.
#' @author Michael Mahr <rchemcpp@@bioinf.jku.at>
#' c++ function written by Jean-Luc Perret and Pierre Mahe
#' 
#' @examples 
#' sdfolder <- system.file("sample_data",package="Rchemcpp")
#' sdf <- list.files(sdfolder,full.names=TRUE,pattern="tiny")
#' moleculeNames <- sd2gramSpectrum(sdf)
#' 
#' @export


sd2gramSpectrum = function(sdFileName, sdFileName2 = "", 
		kernelType = c("spectrum", "tanimoto", "minmaxTanimoto","marginalized","lambda"), 
		margKernelEndProbability = 0.1, lambdaKernelLambda = 1.0, 
		depthMax = as.integer(0), onlyDepthMax = FALSE , flagRemoveH = FALSE, 
		morganOrder = as.integer(0), fileType = c("sd","genericsd","kcf"), 
		silentMode = FALSE, returnNormalized = FALSE, moleculeNameProperty = "",
		moleculeNameProperty2 = "")
{
	margKernelConvgce = 10000;
	
	margKernelSkipSkeleton = FALSE;
	
	if(!is.character(sdFileName)) stop("sdFileName must be string")
	if(!is.character(sdFileName2)) stop("sdFileName2 must be string")
	
	if(!is.character(kernelType)) stop("kernelType must be a string")
	if(!is.numeric(margKernelEndProbability)) stop("margKernelEndProbability must be numeric")
	if(!is.numeric(lambdaKernelLambda)) stop("lambdaKernelLambda must be numeric")
	
	if(!is.numeric(depthMax)) stop("depthMax must be integer")
	depthMax <- as.integer(depthMax)
	if(!is.logical(onlyDepthMax)) stop("onlyDepthMax must be logical")
	if(!is.logical(flagRemoveH)) stop("flagRemoveH must be logical")
	if(!is.numeric(morganOrder)) stop("morganOrder must be integer")
	morganOrder <- as.integer(morganOrder)
	
	if(!is.character(fileType)) stop("fileType must be string")
	if(!is.logical(silentMode)) stop("silentMode must be logical")
	if(!is.logical(returnNormalized)) stop("returnNormalized must be logical")
	
	
	if (missing(fileType)) {fileType = fileType[1]}
	fileType = match.arg(fileType)
	
	if (missing(kernelType)) {kernelType = kernelType[1]}
	kernelType = match.arg(kernelType)
	
	
	#For passing...
	if (kernelType == "marginalized"){
		kernelParam = margKernelEndProbability;
	}else if (kernelType == "lambda"){
		kernelParam = lambdaKernelLambda;
	}else{
		kernelParam = 0.0;
	}
	kernelTypeIndex = as.integer(which(c("spectrum", "tanimoto", "minmaxTanimoto","marginalized","lambda") == kernelType)-1)
	
	
	aSet = new (Rchemcpp::Rmoleculeset);
	aSet2 = new (Rchemcpp::Rmoleculeset);
	
	
	# if sflag = 1 --> SEPARATE TEST SET
	if( sdFileName2 != "" ){
		
		# 1 - data initialization
		# -----------------------
		# read the set of molecules
		if( fileType == "sd" ){
			aSet$addSD( sdFileName, FALSE );
			aSet2$addSD( sdFileName2, FALSE );  
		}else if( fileType == "genericsd" ){
			aSet$addSD( sdFileName, TRUE );
			aSet2$addSD( sdFileName2, TRUE );   
		}else if( fileType == "kcf" ){
			aSet$addKCF( sdFileName );
			aSet2$addKCF( sdFileName2 );
		}
		if( !silentMode ){
			print(paste("*** sd file added ; number of molecules in the set 1 =  ", aSet$numMolecules() ));
			print(paste("*** sd file added ; number of molecules in the set 2 =  ", aSet2$numMolecules() ));
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
		
		# set kashima probabilities if kernelType = marginalized
		if(kernelType == "marginalized"){
			aSet$setKashimaKernelParam( kernelParam, margKernelConvgce, margKernelSkipSkeleton );
			aSet2$setKashimaKernelParam( kernelParam, margKernelConvgce, margKernelSkipSkeleton );
		}
		
		# initialize gram matrices
		aSet2$initializeSelfKernel( 0.0);  
		aSet$initializeSelfKernel( 0.0);   
		aSet$setComparisonSetCopy( aSet2 );
		aSet$initializeGram( 0.0 );
		
		
		if( !silentMode){
			print("#### initialization Gram OK");
		}
		
		
		# 3 - compute the gram matrix
		# ----------------------------
		gramSpectrum_test( aSet, depthMax, kernelTypeIndex, kernelParam, onlyDepthMax, silentMode);
		if( !silentMode ){
			print("gramComputeSpectrum (test) OK");
		}
		
		# normalize gram
		if (kernelType == "tanimoto")
		{
			aSet$normalizeTanimoto();
		}
		else if (kernelType == "minmaxTanimoto") 
		{
			aSet$normalizeTanimotoMinMax();
		}
		else
		{
			aSet$normalizeGram();
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
		
		# 1 - data initialization
		# -----------------------
		# read the set of molecules  
		if( fileType == "sd" ){
			aSet$addSD( sdFileName, FALSE );
		}else if( fileType == "genericsd" ){
			aSet$addSD( sdFileName, TRUE );
		}else if( fileType == "kcf" ){
			aSet$addKCF( sdFileName );
		}
		if( !silentMode ){
			print (paste("*** sd file added ; number of molecules =  ", aSet$numMolecules() ));
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
		
		# set kashima probabilities if kernelType = marginalized
		if(kernelType == "marginalized"){
			aSet$setKashimaKernelParam( kernelParam, margKernelConvgce, margKernelSkipSkeleton );
		}
		
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
		gramSpectrum_self( aSet, depthMax, kernelTypeIndex, kernelParam, onlyDepthMax, silentMode);
		if( !silentMode ){
			print("gramComputeSpectrum (self) OK");
		}
		
		# normalize gram
		if (kernelType == "tanimoto")
		{
			aSet$normalizeTanimoto();
		}
		else if (kernelType == "minmaxTanimoto") 
		{
			aSet$normalizeTanimotoMinMax();
		}
		else
		{
			aSet$normalizeGram();
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
	
	
	#name the molecules
	if (moleculeNameProperty != "")
	{
		molnames = c()
		for (i in 0:(aSet$numMolecules() -1))
		{
			mol = aSet$getMolByIndex(i);
			
			if( moleculeNameProperty %in% mol$listStringDescriptors() )
			{
				molnames = c(molnames, mol$getStringDescriptorValue(moleculeNameProperty))
			}
			else
			{
				molnames = c(molnames, i)
			}	
		}
		rownames(K) <- molnames

		if (sdFileName2==""){
			colnames(K) <- molnames
		} else {
			molnames2 = c()
			for (i in 0:(aSet2$numMolecules() -1))
			{
				mol2 = aSet2$getMolByIndex(i);
			
				if( moleculeNameProperty2 %in% mol2$listStringDescriptors() )
				{
					molnames2 = c(molnames2, mol2$getStringDescriptorValue(moleculeNameProperty2))
				}
				else
				{
					molnames2 = c(molnames2, i)
				}	
			}
			colnames(K) <- molnames2
		}
			
	}
	
	return ( K )
	
}



