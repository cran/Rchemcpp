#' @title sd2gram - Similarity of molecules by the marginalized kernel and 
#' proposed extensions.
#' 
#' This tools compute the marginalized kernel (\cite{Kashima, 2004})
#' and its proposed extensions (\cite{Mahe, 2005)}.
#' 
#' @param sdFileName File containing the molecules. Must be in MDL file format
#' (MOL and SDF files). For more information on the file format see 
#' http://en.wikipedia.org/wiki/Chemical_table_file. Default = "missing".
#' @param sdFileName2 A second file containing molecules. Must also be in SDF.
#' If specified the molecules of the first file will be compared with the 
#' molecules of this second file. Default = "missing".
#' @param stopP ... . Default = 0.1.
#' @param filterTottering A logical specifying whether tottering paths should
#' be removed. Default = FALSE.
#' @param converg A numeric value specifying when convergence is reached. The
#' algorithm stops when the kernel value does not change by more
#' than 1/c, where c is the value specified by the converg option. 
#' Default = 1000.
#' @param atomKernelMatrix A string that sets the similarity measure between
#' atoms that should be used. Dfault = "missing". 
#' @param flagRemoveH A logical that indicates whether H-atoms should be 
#' removed or not.
#' @param morganOrder The order of the DeMorgan Indices to be used. If set to
#' zero no DeMorgan Indices are used. The higher the order the more different
#' types of atoms exist and consequently the more dissimilar will be the molecules.
#' @param fileType Which filetype was submitted.
#' @param silentMode Whether or not the program should print progress reports
#' to the standart output. Default = FALSE.
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
#' @references
#' (Kashima, 2004) -- H. Kashima, K. Tsuda, and A. Inokuchi. Kernels for graphs. 
#' In B. Schoelkopf, K. Tsuda, and J.P.
#' Vert, editors, Kernel Methods in Computational Biology, pages 155-170.
#' MIT Press, 2004.
#' 
#' (Mahe, 2005) -- P. Mahe, N. Ueda, T. Akutsu, J.-L. Perret, and J.-P. Vert. 
#' Graph kernels for molecular structure-
#' activity relationship analysis with support vector machines. 
#' J Chem Inf Model, 45(4):939-51, 2005.
#' @export


sd2gram = function(sdFileName, sdFileName2 = "", stopP = 0.1, 
		filterTottering = FALSE, converg = as.integer(1000),  
		atomKernelMatrix = "", flagRemoveH = FALSE, morganOrder = as.integer(0),
		fileType = c("sd","genericsd","kcf"), silentMode = FALSE, 
		returnNormalized = TRUE)
{
	
	nbThreadsWanted = as.integer(1)
	fromN = as.integer(1)
	
	if(!is.numeric(stopP)) stop("stopP must be numeric")
	if(!is.logical(filterTottering)) stop("filterTottering must be logical")
	if(!is.numeric(converg)) stop("converg must be integer")
	converg <- as.integer(converg)
	if(!is.character(sdFileName)) stop("sdFileName must be string")
	if(!is.character(sdFileName2)) stop("sdFileName2 must be string")
	if(!is.character(atomKernelMatrix)) stop("atomKernelMatrix must be string")
	if(!is.logical(flagRemoveH)) stop("flagRemoveH must be logical")
	
	#if(!is.numeric(nbThreadsWanted)) stop("nbThreadsWanted must be integer")
	#nbThreadsWanted <- as.integer(nbThreadsWanted)
	if(!is.numeric(morganOrder)) stop("morganOrder must be integer")
	morganOrder <- as.integer(morganOrder)
	#if(!is.numeric(fromN)) stop("fromN must be integer")
	#fromN <- as.integer(fromN)
	
	if(!is.character(fileType)) stop("fileType must be string")
	if(!is.logical(silentMode)) stop("silentMode must be logical")
	if(!is.logical(returnNormalized)) stop("returnNormalized must be logical")
	
	if (missing(fileType)) {fileType = fileType[1]}
	fileType = match.arg(fileType)
	
	
	aSet = new (Rchemcpp::Rmoleculeset);
	aSet2 = new (Rchemcpp::Rmoleculeset);
	
	if( atomKernelMatrix != "" ){
		loadGramAtoms( atomKernelMatrix );	# atomKernelMatrix is set by the
		# -k argument on the command line
	}
	
	
	if( !silentMode ){
		print("reading file");
	}
	
	if( fileType == "sd" ){
		aSet$addSD( sdFileName, FALSE);
		if( flagRemoveH == TRUE ) {
			if( !silentMode ){
				print("removing Hydrogens");
			}
			aSet$hideHydrogens();
		}
		
	}else if( fileType == "genericsd" ){
		aSet$addSD( sdFileName, TRUE);
		
	}else if( fileType == "kcf" ){
		aSet$addKCF( sdFileName );
		
		if( flagRemoveH == TRUE ) {
			if( !silentMode ){
				print("removing Hydrogens");
			}
			aSet$hideHydrogens();
		}
	}
	
	if( !silentMode ){
		print("reading file done");
	}
	
	if(sdFileName2 == ""){
		# if the gram matrix is to be computed with all mol files in a single directory
		
		if( !silentMode ){
			print("setting morgan labels");
		}
		aSet$setMorganLabels(morganOrder);
		
		
		if( atomKernelMatrix != "" ){
			
			if( !silentMode ){
				print("using external atomKernel");
			}
			
			
			aSet$gramCompute(stopP, converg, fromN, nbThreadsWanted, silentMode, filterTottering, TRUE); #external			
			
		}else{
			if( !silentMode ){
				print("using moleculeKernel Kashima");
			}
			
			aSet$gramCompute(stopP, converg, fromN, nbThreadsWanted, silentMode, filterTottering, FALSE); #morgan			
		}
		
		
		# write the gram matrix to a file
		# raw matrix and normalised matrix
		
		
		#aSet$deleteAll(); #Should be done by destructor!!!
		
		
	}else{
		if( !silentMode ){
			print("test set mode");
		}
		
		
		if( !silentMode ){
			print("reading second file");
		}
		
		if( fileType == "sd" ){
			aSet2$addSD( sdFileName2, FALSE);
			if( flagRemoveH == TRUE ) {
				if( !silentMode ){
					print("removing Hydrogens");
				}
				aSet2$hideHydrogens();
			}
			
		}else if( fileType == "genericsd" ){
			aSet2$addSD( sdFileName2, TRUE );
			
		}else if( fileType == "kcf" ){
			aSet2$addKCF( sdFileName2 );
			
			if( flagRemoveH == TRUE ) {
				if( !silentMode ){
					print("removing Hydrogens");
				}
				aSet2$hideHydrogens();
			}
		}
		if( !silentMode ){
			print("reading second file done");
		}
		
		
		# if the gram matrix is to be computed between two datasets of mol files
		# in two separate directories
		#deleted
		
		
		if( !silentMode ){
			print ("setting morgan labels");
		}
		
		aSet$setMorganLabels( morganOrder );
		aSet2$setMorganLabels( morganOrder );
		
		if( atomKernelMatrix != "" ){
			
			if( !silentMode ){
				print("using external atomKernel");
			}
			
			aSet$gramCompute(stopP, converg, fromN, nbThreadsWanted, silentMode, filterTottering, TRUE); #external			
			
		}else{
			
			if( !silentMode ){
				print("using moleculeKernel Kashima");
			}
			
			aSet$gramCompute(stopP, converg, fromN, nbThreadsWanted, silentMode, filterTottering, FALSE); #morgan			
			
		}
		
		
		#aSet$deleteAll(); #Should be done by destructor!!!
		#aSet2$deleteAll(); #Should be done by destructor!!!
	}
	
	if( !silentMode ){
		print("end");
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
	
	
}



