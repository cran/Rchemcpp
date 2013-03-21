#' @title sd2gramSubtree - Similarity of molecules by several graph kernels
#' based on the count of common subtrees
#'
#' @description This tools computes several graph kernels based on the detection of common 
#' subtrees: the so-called
#' tree-pattern graph kernels, originally introduced in (\cite{Ramon, 2003}), 
#' and revisited in (\cite{Mahe, 2006}). 
#' 
#' 
#' @param sdf File containing the molecules. Must be in MDL file format
#' (MOL and SDF files). For more information on the file format see 
#' http://en.wikipedia.org/wiki/Chemical_table_file. Default = "missing"
#' @param sdf2 A second file containing molecules. Must also be in SDF.
#' If specified the molecules of the first file will be compared with the 
#' molecules of this second file. Default = "missing".
#' @param kernelType Determines whether subtrees of the molecule are penalized
#' size-based or branching-based. Default = "sizebased".
#' @param branchKernelUntilN Logical whether tree patterns of until N should be 
#' considered. Default = FALSE.
#' @param lambda Weighted contribution of tree-patterns depending on their sizes
#' Default = 1.
#' @param depthMax tree-patterns of depth. Default = 3. 
#' @param flagRemoveH A logical that indicates whether H-atoms should be 
#' removed or not.
#' @param filterTottering A logical that indicates whether tottering paths
#' should be removed. Default = FALSE.
#' @param morganOrder The order of the DeMorgan Indices to be used. If set to
#' zero no DeMorgan Indices are used. The higher the order the more different
#' types of atoms exist and consequently the more dissimilar will be the molecules.
#' @param silentMode Whether the program should print progress reports
#' to the standart output.
#' @param returnNormalized A logical specifying whether a normalized kernel
#' matrix should be returned. Default = TRUE.
#' @param detectArom Whether aromatic rings should be detected and aromatic
#' bonds should a special bond type. (Default = TRUE).
#' @examples 
#' sdfolder <- system.file("sample_data",package="Rchemcpp")
#' sdf <- list.files(sdfolder,full.names=TRUE,pattern="small")
#' K <- sd2gram(sdf)
#' @return A numeric matrix containing the similarity values between the
#' molecules.
#' @author Michael Mahr <rchemcpp@@bioinf.jku.at>
#' c++ function written by Jean-Luc Perret and Pierre Mahe
#' @references 
#' (Mahe, 2006) -- P. Mahe and J.-P. Vert. Graph kernels based on tree patterns 
#' for molecules. Technical Report, HAL:ccsd-00095488, Ecoles des Mines de Paris,
#' September 2006.
#' (Ramon, 2003) -- J. Ramon and T. Gaertner. Expressivity versus efficiency of 
#' graph kernels. In T. Washio and L. De Raedt, editors, 
#' Proceedings of the First International Workshop on Mining Graphs, 
#' Trees and Sequences, pages 65-74, 2003.
#' 
#' @export


sd2gramSubtree = function(sdf, sdf2, 
		kernelType = c("sizebased","branchingbased") , 
		branchKernelUntilN = FALSE, lambda = 1,  depthMax = as.integer(3), 
		flagRemoveH = FALSE, filterTottering = FALSE, morganOrder = as.integer(0),
		silentMode = FALSE, 
		returnNormalized = TRUE,detectArom=TRUE)
{
	
	branchFlag = FALSE;
	
		
	if(!is.logical(branchKernelUntilN)) stop("branchKernelUntilN must be logical")
	if(!is.numeric(lambda)) stop("lambda must be numeric")
	if(!is.numeric(depthMax)) stop("depthMax must be integer")
	depthMax <- as.integer(depthMax)
	if(!is.logical(flagRemoveH)) stop("flagRemoveH must be logical")
	
	if(!is.logical(filterTottering)) stop("filterTottering must be logical")
	if(!is.numeric(morganOrder)) stop("morganOrder must be integer")
	morganOrder <- as.integer(morganOrder)
	
	if(!is.logical(silentMode)) stop("silentMode must be logical")
	if(!is.logical(returnNormalized)) stop("returnNormalized must be logical")
	
	if (missing(kernelType)) {kernelType = kernelType[1]}
	kernelType = match.arg(kernelType)
	
	
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
	
	
	if(kernelType == 1)
		branchFlag = TRUE;
	
	
	# JP: For the fast version we need a static variable to store the list of subsets of p indeices among n choices
	# The static variable tuples is initialized once with the initialize_tuples( int nmax) function
	# After initialization, tuples[n][p] is a vector of vector<int> that contains all different p-tuples of {1,..,n}.
	initialize_tuples(4);
	
	if( !missing(sdf2) ){
	
		# remove H when specified on command line
		if( flagRemoveH == TRUE ){
			
			if( !silentMode ){
				print( "removing hydrogens" );
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
		
		# no-totters transformation
		if(filterTottering){
			aSet$noTottersTransform();
			aSet2$noTottersTransform();
		}
		
		# initialize gram matrices
		aSet$initializeSelfKernel( 0.0);
		aSet2$initializeSelfKernel( 0.0);
		aSet$setComparisonSetCopy( aSet2 );
		aSet$initializeGram( 0.0 );
		if( !silentMode)
			print("initialization Gram OK");  
		
		
		# Initialize the indexes
		#aSet$initialize_extended(); #moved to function
		#aSet2$initialize_extended(); #moved to function
		
		
		gramSubtree_test( aSet, lambda, depthMax, filterTottering, branchFlag, branchKernelUntilN, silentMode);
		
		
		if( !silentMode ){
			print("compute gram OK");
		}
		
		# normalize gram
		aSet$normalizeGram();
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
		#aSet2$writeSelfKernelList( outputDir + baseName + "_train", silentMode );
		
		
	}else{ #  SELF KERNEL
			
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
		
		# no-totters transformation
		if(filterTottering){
			aSet$noTottersTransform();
		}
		
		# initialize gram matrices
		aSet$setComparisonSetSelf();
		aSet$initializeGram( 0.0 );
		aSet$initializeSelfKernel( 0.0 );
		if( !silentMode){
			print("initialization Gram OK");
		}
		
		
		#aSet$initialize_extended(); #moved to function
		
		gramSubtree_self(aSet, lambda, depthMax, filterTottering, branchFlag, branchKernelUntilN, silentMode);
		
		if( !silentMode )
			print("compute gram OK");
		
		# normalize gram
		aSet$normalizeGram();
		if( !silentMode )
			print("normalize gram (self) OK");
		# write gram matrices and self kernel values
		
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




