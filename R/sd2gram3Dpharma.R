
#' @title sd2gram3Dpharma - Similarity of molecules by the exact pharmacophore
#' kernel.
#' 
#' This tool implements the (exact version of) pharmacophore kernel for 3D 
#' structures of molecules (\cite{Mahe, 2006}).
#' 
#' @param sdFileName File containing the molecules. Must be in MDL file format
#' (MOL and SDF files). For more information on the file format see 
#' http://en.wikipedia.org/wiki/Chemical_table_file. Default = "missing".
#' @param sdFileName2 A second file containing molecules. Must also be in SDF.
#' If specified the molecules of the first file will be compared with the 
#' molecules of this second file. Default = "missing".
#' @param chargesFileName A character with the name of the file containing
#' the atom charges. Default = missing.
#' @param chargesFileName2 A character with the name of the file containing
#' the atom charges. Default = missing.
#' @param edgeKernelType Options to specify the kernel function comparing 
#' distances between atoms. Choices are "RBF" or "triangular". Default = "RBF".
#' @param edgeKernelParameter Specifies the parameter associated to these
#' kernels. Either the bandwith of the RBF kernel or the cut-off of the
#' triangular kernel. Default = 1.
#' @param atomKernelMatrix A string that sets the similarity measure between
#' atoms that should be used. Dfault = "missing". 
#' @param flagRemoveH A logical that indicates whether H-atoms should be 
#' removed or not.
#' @param morganOrder The order of the DeMorgan Indices to be used. If set to
#' zero no DeMorgan Indices are used. The higher the order the more different
#' types of atoms exist and consequently the more dissimilar will be the molecules.
#' @param morganChargesThreshold specifies a threshold above which partial 
#' Morgan charges are considered as positive/negative.
#' @param fileType Which filetype was submitted.
#' @param silentMode Whether or not the program should print progress reports
#' to the standart output.
#' @param returnNormalized A logical specifying whether a normalized kernel
#' matrix should be returned. Default = TRUE.
#' @examples 
#' sdfolder <- system.file("sample_data",package="Rchemcpp")
#' sdf <- list.files(sdfolder,full.names=TRUE,pattern="small")
#' K <- sd2gram(sdf)
#' @return A numeric matrix containing the similarity values between the
#' molecules.
#' @author Michael Mahr <rchemcpp@@bioinf.jku.at>
#' c++ function written by Jean-Luc Perret and Pierre Mahe
#' @references (Mahe, 2006) --  P. Mahe, L. Ralaivola, V. Stoven, and J.-P. Vert.
#' The pharmacophore kernel for virtual screening
#' with support vector machines. Technical Report, HAL:ccsd-00020066, Ecole des
#' Mines de Paris, March 2006.
#' @export


sd2gram3Dpharma = function(sdFileName, sdFileName2 = "", chargesFileName = "", 
		chargesFileName2 = "",  edgeKernelType = c("RBF", "triangular"), 
		edgeKernelParameter = 1, atomKernelMatrix = "", flagRemoveH = FALSE, 
		morganOrder = as.integer(0), morganChargesThreshold = 0, 
		fileType = c("sd","genericsd","kcf"), silentMode = FALSE, 
		returnNormalized = FALSE)
{
	
	if(!is.character(sdFileName)) stop("sdFileName must be string")
	if(!is.character(sdFileName2)) stop("sdFileName2 must be string")
	if(!is.character(chargesFileName)) stop("chargesFileName must be string")
	if(!is.character(chargesFileName2)) stop("chargesFileName2 must be string")
	
	if(!is.character(edgeKernelType)) stop("edgeKernelType must be a string")
	if(!is.numeric(edgeKernelParameter)) stop("edgeKernelParameter must be numeric")
	if(!is.character(atomKernelMatrix)) stop("atomKernelMatrix must be string")
	
	if(!is.logical(flagRemoveH)) stop("flagRemoveH must be logical")
	if(!is.numeric(morganOrder)) stop("morganOrder must be integer")
	morganOrder <- as.integer(morganOrder)
	if(!is.numeric(morganChargesThreshold)) stop("morganchargesThreshold must be numeric")
	
	if(!is.character(fileType)) stop("fileType must be string")
	if(!is.logical(silentMode)) stop("silentMode must be logical")
	if(!is.logical(returnNormalized)) stop("returnNormalized must be logical")
	
	if (missing(fileType)) {fileType = fileType[1]}
	fileType = match.arg(fileType)
	
	if (missing(edgeKernelType)) {edgeKernelType = edgeKernelType[1]}
	edgeKernelType = match.arg(edgeKernelType)
	
	#--
	
	if((sdFileName2 == "") && (chargesFileName2 != ""))
		print("chargesFilename2 is useless") 
	
	if((sdFileName2 != "") && (chargesFileName == "") && (chargesFileName2 != "")) 
		stop("both chargesFileNames have to be specified") 
	if((sdFileName2 != "") && (chargesFileName != "") && (chargesFileName2 == ""))
		stop("both chargesFileNames have to be specified") 
	
	
	
	if( atomKernelMatrix != "" ){
		loadGramAtoms( atomKernelMatrix );	# atomKernelMatrix is set by the
	}
	
	
	if( sdFileName2 != "" ){ # i.e., test set version.
		
		aSet = new (Rchemcpp::Rmoleculeset);
		aSet2 = new (Rchemcpp::Rmoleculeset);
		
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
		
		# read partial charges if specified on command line
		if(chargesFileName != ""){
			if( !silentMode ){
				print("reading partial charges");
			}
			aSet$readPartialCharges(chargesFileName);
			aSet2$readPartialCharges(chargesFileName2);
		}
		
		# remove hydrogens
		if( flagRemoveH == TRUE ){
			if( !silentMode ){
				print("removing hydrogens");
			}
			aSet$hideHydrogens();
			aSet2$hideHydrogens();
		}
		
		# set morgan labels
		if( !silentMode )
			print(paste("setting morgan labels ", morganOrder));
		aSet$setMorganLabels( morganOrder );
		aSet2$setMorganLabels( morganOrder );
		
		# introduce partial charges in Morgan labels if specified on command line
		if(chargesFileName != ""){
			if( !silentMode ){
				print("setting partial charges");
			}
			aSet$setMorganChargesLabels(morganChargesThreshold);
			aSet2$setMorganChargesLabels(morganChargesThreshold);
		}
		
		# compute gram matrix
		aSet$setComparisonSetCopy( aSet2 );
		
		if( atomKernelMatrix != "" ){
			if( !silentMode ){
				print("using external atomKernel");
			}
			
			if (edgeKernelType == "RBF"){
				# --> RBF edge kernel
				aSet$gramCompute3D (edgeKernelParameter, silentMode, TRUE, TRUE) #external,rbf
				
			}else if (edgeKernelType == "triangular"){
				# --> triangular edge kernel
				aSet$gramCompute3D (edgeKernelParameter, silentMode, TRUE, FALSE) #external,triangular
			}
		}
		else{
			if (edgeKernelType == "RBF"){
				# --> RBF edge kernel
				aSet$gramCompute3D (edgeKernelParameter, silentMode, FALSE, TRUE) #morgan,rbf
				
			}else if (edgeKernelType == "triangular"){
				# --> triangular edge kernel
				aSet$gramCompute3D (edgeKernelParameter, silentMode, FALSE, FALSE) #morgan,triangular
				
			}
		}
		
		
		if( !silentMode ){
			print("gram matrix computation OK");
		}
		
		# write the gram matrix to a file (raw and normalized matrices)
		#aSet$writeGramMatrix( outputDir + baseName, false, false, silentMode);
		#aSet$writeGramMatrix( outputDir + baseName, true, false, silentMode);
		# --> write self kernel values
		
		if (returnNormalized == FALSE)
		{	
			K <- do.call(rbind,aSet$getGram() )
			
		}
		else
		{		
			K <- do.call(rbind,aSet$getGramNormal() )
			
		}
		
		
		
		xx <- try(molNames1 <- getMoleculeNamesFromSDF(sdFileName))
		if (inherits(xx,"try-error") | length(molNames1)!=nrow(K)){
			molNames1 <- paste("Mol",1:nrow(K),sep="")
		}
		
		if (sdFileName2==""){
			molNames2 <- molNames1
		} else {
			yy <- try(molNames2 <- getMoleculeNamesFromSDF(sdFileName2))
			if (inherits(yy,"try-error") | length(molNames2)!=ncol(K)){
				molNames2 <- paste("Mol",(nrow(K)+1):(nrow(K)+ncol(K)))
			}
		}
		
		
		rownames(K) <- molNames1
		colnames(K) <- molNames2
		
		return ( K )
		
		
		#aSet$writeSelfKernelList( outputDir + baseName + "_test", silentMode );
		#aSet2$writeSelfKernelList( outputDir + baseName + "_train", silentMode );
		# and exit
		#aSet$deleteAll();   should be done by destructor
		#aSet2$deleteAll();  should be done by destructor
		
		
	}
	else{ # self-set version
		
		# create a new moleculeSet
		aSet = new (Rchemcpp::Rmoleculeset);
		
		if( fileType == "sd" ){
			aSet$addSD( sdFileName, FALSE );
		}else if( fileType == "genericsd" ){
			aSet$addSD( sdFileName, TRUE );
		}else if( fileType == "kcf" ){
			aSet$addKCF( sdFileName );
		}
		
		# read partial charges if specified on command line
		if(chargesFileName != ""){
			if( !silentMode )
				print("reading partial charges");
			aSet$readPartialCharges(chargesFileName);
		}
		
		# remove hydrogens
		if( flagRemoveH == TRUE ) {
			if( !silentMode ){
				print("removing Hydrogens");
			}
			aSet$hideHydrogens();
		}
		# set morgan labels
		if( !silentMode ){
			print(paste("setting morgan labels ", morganOrder));
		}
		aSet$setMorganLabels(morganOrder);
		if( !silentMode ){
			print("setting morgan labels OK");
		}
		# introduce partial charges in Morgan labels if specified on command line
		if(chargesFileName != ""){
			if( !silentMode ){
				print("setting partial charges");
			}
			aSet$setMorganChargesLabels(morganChargesThreshold);
		}   
		
		# compute gram matrix
		aSet$setComparisonSetSelf();
		
		if( atomKernelMatrix != "" ){
			if( !silentMode ) {
				print("using external atomKernel");
			}
			
			if (edgeKernelType == "RBF"){
				# --> RBF edge kernel
				aSet$gramCompute3D (edgeKernelParameter, silentMode, TRUE, TRUE) #external,rbf
				
			}else if (edgeKernelType == "triangular"){
				# --> triangular edge kernel
				aSet$gramCompute3D (edgeKernelParameter, silentMode, TRUE, FALSE) #external,triangular
			}
		}
		else{
			if (edgeKernelType == "RBF"){
				# --> RBF edge kernel
				aSet$gramCompute3D (edgeKernelParameter, silentMode, FALSE, TRUE) #morgan,rbf
				
			}else if (edgeKernelType == "triangular"){
				# --> triangular edge kernel
				aSet$gramCompute3D (edgeKernelParameter, silentMode, FALSE, FALSE) #morgan,triangular
				
			}
			
		}
		
		if( !silentMode )
		{
			print("gram matrix computation OK");
		}
		
		
		if (returnNormalized == FALSE)
		{	
			K <- do.call(rbind,aSet$getGram() )
			
		}
		else
		{		
			K <- do.call(rbind,aSet$getGramNormal() )
			
		}
		
		
		
		xx <- try(molNames1 <- getMoleculeNamesFromSDF(sdFileName))
		if (inherits(xx,"try-error") | length(molNames1)!=nrow(K)){
			molNames1 <- paste("Mol",1:nrow(K),sep="")
		}
		
		if (sdFileName2==""){
			molNames2 <- molNames1
		} else {
			yy <- try(molNames2 <- getMoleculeNamesFromSDF(sdFileName2))
			if (inherits(yy,"try-error") | length(molNames2)!=ncol(K)){
				molNames2 <- paste("Mol",(nrow(K)+1):(nrow(K)+ncol(K)))
			}
		}
		
		
		rownames(K) <- molNames1
		colnames(K) <- molNames2
		
		return ( K )
		
		
		# write the gram matrix to a file (raw and normalized matrices)
		#aSet$writeGramMatrix( outputDir + baseName, false, false, silentMode);
		#aSet$writeGramMatrix( outputDir + baseName, true, false, silentMode);
		# --> write self kernel values
		#aSet$writeSelfKernelList( outputDir + baseName + "_train", false);      
		# and exit
		#aSet$deleteAll();  should be done by destructor
	}
	
	
	
}


