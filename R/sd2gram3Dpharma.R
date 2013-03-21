
#' @title sd2gram3Dpharma - Similarity of molecules by the exact pharmacophore
#' kernel.
#' 
#' @description This tool implements the (exact version of) pharmacophore kernel for 3D 
#' structures of molecules (\cite{Mahe, 2006}).
#'  
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
#' @param silentMode Whether or not the program should print progress reports
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
#' @references (Mahe, 2006) --  P. Mahe, L. Ralaivola, V. Stoven, and J.-P. Vert.
#' The pharmacophore kernel for virtual screening
#' with support vector machines. Technical Report, HAL:ccsd-00020066, Ecole des
#' Mines de Paris, March 2006.
#' 
#' @export


sd2gram3Dpharma = function(sdf, sdf2, chargesFileName = "", 
		chargesFileName2 = "",  edgeKernelType = c("RBF", "triangular"), 
		edgeKernelParameter = 1, atomKernelMatrix = "", flagRemoveH = FALSE, 
		morganOrder = as.integer(0), morganChargesThreshold = 0, silentMode = FALSE, 
		returnNormalized = TRUE, detectArom=TRUE)
{
	
	if(!is.character(chargesFileName)) stop("chargesFileName must be string")
	if(!is.character(chargesFileName2)) stop("chargesFileName2 must be string")
	
	if(!is.character(edgeKernelType)) stop("edgeKernelType must be a string")
	if(!is.numeric(edgeKernelParameter)) stop("edgeKernelParameter must be numeric")
	if(!is.character(atomKernelMatrix)) stop("atomKernelMatrix must be string")
	
	if(!is.logical(flagRemoveH)) stop("flagRemoveH must be logical")
	if(!is.numeric(morganOrder)) stop("morganOrder must be integer")
	morganOrder <- as.integer(morganOrder)
	if(!is.numeric(morganChargesThreshold)) stop("morganchargesThreshold must be numeric")
	
	if(!is.logical(silentMode)) stop("silentMode must be logical")
	if(!is.logical(returnNormalized)) stop("returnNormalized must be logical")
	

	if (missing(edgeKernelType)) {edgeKernelType = edgeKernelType[1]}
	edgeKernelType = match.arg(edgeKernelType)
	
	#--
	
	if((missing(sdf2)) && (chargesFileName2 != ""))
		print("chargesFilename2 is useless") 
	
	if((missing(sdf2)) && (chargesFileName == "") && (chargesFileName2 != "")) 
		stop("both chargesFileNames have to be specified") 
	if((missing(sdf2)) && (chargesFileName != "") && (chargesFileName2 == ""))
		stop("both chargesFileNames have to be specified") 
	if( atomKernelMatrix != "" ){
		loadGramAtoms( atomKernelMatrix );	# atomKernelMatrix is set by the
	}
	
	
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
	
	
	if(!missing(sdf2)){ # i.e., test set version.
		
		
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
		
		
		#aSet$writeSelfKernelList( outputDir + baseName + "_test", silentMode );
		#aSet2$writeSelfKernelList( outputDir + baseName + "_train", silentMode );
		# and exit
		#aSet$deleteAll();   should be done by destructor
		#aSet2$deleteAll();  should be done by destructor
		
	}
	else{ # self-set version
		
		
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
		
		
		# write the gram matrix to a file (raw and normalized matrices)
		#aSet$writeGramMatrix( outputDir + baseName, false, false, silentMode);
		#aSet$writeGramMatrix( outputDir + baseName, true, false, silentMode);
		# --> write self kernel values
		#aSet$writeSelfKernelList( outputDir + baseName + "_train", false);      
		# and exit
		#aSet$deleteAll();  should be done by destructor
	}
	
	
	rownames(K) <- molnames
	colnames(K) <- molnames2
	
	
	return ( K )
	
}


