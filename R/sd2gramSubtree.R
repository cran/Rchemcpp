#' @title sd2gramSubtree - Similarity of molecules by several graph kernels
#' based on the count of common subtrees
#'
#' @description This tools computes several graph kernels based on the detection of common 
#' subtrees: the so-called
#' tree-pattern graph kernels, originally introduced in (\cite{Ramon, 2003}), 
#' and revisited in (\cite{Mahe, 2006}). 
#' 
#' @usage sd2gramSubtree(sdFileName, sdFileName2 = "", 
#'		kernelType = c("sizebased","branchingbased") , 
#'		branchKernelUntilN = FALSE, lambda = 1,  depthMax = as.integer(3), 
#'		flagRemoveH = FALSE, filterTottering = FALSE, morganOrder = as.integer(0), 
#'		fileType = c("sd", "genericsd","kcf"), silentMode = FALSE, 
#'		returnNormalized = FALSE, moleculeNameProperty = "",
#'		moleculeNameProperty2 = "")
#' 
#' 
#' @param sdFileName File containing the molecules. Must be in MDL file format
#' (MOL and SDF files). For more information on the file format see 
#' http://en.wikipedia.org/wiki/Chemical_table_file. Default = "missing"
#' @param sdFileName2 A second file containing molecules. Must also be in SDF.
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
#' @param fileType Which filetype was submitted.
#' @param silentMode Whether the program should print progress reports
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
#' K <- sd2gram(sdf, moleculeNameProperty="Compound Name")
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
#' @examples 
#' sdfolder <- system.file("sample_data",package="Rchemcpp")
#' sdf <- list.files(sdfolder,full.names=TRUE,pattern="tiny")
#' moleculeNames <- sd2gramSubtree(sdf)
#' 
#' @export


sd2gramSubtree = function(sdFileName, sdFileName2 = "", 
		kernelType = c("sizebased","branchingbased") , 
		branchKernelUntilN = FALSE, lambda = 1,  depthMax = as.integer(3), 
		flagRemoveH = FALSE, filterTottering = FALSE, morganOrder = as.integer(0), 
		fileType = c("sd", "genericsd","kcf"), silentMode = FALSE, 
		returnNormalized = FALSE, moleculeNameProperty = "",
		moleculeNameProperty2 = "")
{
	
	branchFlag = FALSE;
	
	if(!is.character(sdFileName)) stop("sdFileName must be string")
	if(!is.character(sdFileName2)) stop("sdFileName2 must be string")
	
	if(!is.logical(branchKernelUntilN)) stop("branchKernelUntilN must be logical")
	if(!is.numeric(lambda)) stop("lambda must be numeric")
	if(!is.numeric(depthMax)) stop("depthMax must be integer")
	depthMax <- as.integer(depthMax)
	if(!is.logical(flagRemoveH)) stop("flagRemoveH must be logical")
	
	if(!is.logical(filterTottering)) stop("filterTottering must be logical")
	if(!is.numeric(morganOrder)) stop("morganOrder must be integer")
	morganOrder <- as.integer(morganOrder)
	
	if(!is.character(fileType)) stop("fileType must be string")
	if(!is.logical(silentMode)) stop("silentMode must be logical")
	if(!is.logical(returnNormalized)) stop("returnNormalized must be logical")
	
	
	
	if (missing(fileType)) {fileType = fileType[1]}
	fileType = match.arg(fileType)
	
	if (missing(kernelType)) {kernelType = kernelType[1]}
	kernelType = match.arg(kernelType)
	
	
	
	aSet = new (Rchemcpp::Rmoleculeset);
	aSet2 = new (Rchemcpp::Rmoleculeset);
	
	
	if(kernelType == 1)
		branchFlag = TRUE;
	
	
	# JP: For the fast version we need a static variable to store the list of subsets of p indeices among n choices
	# The static variable tuples is initialized once with the initialize_tuples( int nmax) function
	# After initialization, tuples[n][p] is a vector of vector<int> that contains all different p-tuples of {1,..,n}.
	initialize_tuples(4);
	
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
			print(paste("*** sd file added ; number of molecules in the set 1 =  ", aSet$numMolecules() ) );
			print(paste("*** sd file added ; number of molecules in the set 2 =  ", aSet2$numMolecules() ) );
		}
		
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
		if( !silentMode )
			print(paste("*** sd file added ; number of molecules in the set 1 =  ", aSet$numMolecules() ) );
		
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




