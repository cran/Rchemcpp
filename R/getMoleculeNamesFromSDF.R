#' @title getMoleculeNamesFromSDF - a helper function
#' 
#' @description This function helps to extract a certain property from an SDF file. Usually
#' the molecule class, like "active/non-active" or a property of the molecule,
#' like "biological activity", is also stored in the SDF file. These values
#' often serve as targets for a prediction task. This function is a small
#' wrapper that extracts the information.
#' 
#' @usage 
#' getMoleculeNamesFromSDF(sdfile)
#' 
#' @param sdfile A character containing the name of the SDF file.
#' @return A character vector with one name per molecule.
#' @author Guenter Klambauer <rchemcpp@@bioinf.jku.at>
#' 
#' @examples 
#' sdfolder <- system.file("sample_data",package="Rchemcpp")
#' sdf <- list.files(sdfolder,full.names=TRUE,pattern="small")
#' moleculeNames <- getMoleculeNamesFromSDF(sdf)
#' 
#' @export


getMoleculeNamesFromSDF <- function(sdfile){
	rl <- readLines(sdfile)
	mnames <- rl[grep("V[2|3]000",rl)-3]
	if (length(unique(mnames))!=length(mnames)){
		message("Molecule names not unique. Adding ID number.")
		mnames <- paste(mnames,1:length(mnames),sep="_")
	}
	#check if number of molecules is ok
	if (length(grep("\\$\\$\\$\\$",rl))!=length(mnames)){
		stop("Inconsistency in SDF/Mol file. Could not read names.")
	}
	
	return(mnames)
}



#' @title getMoleculePropertyFromSDF - a helper function
#' 
#' @description This function helps to extract a certain property from an SDF 
#' file. Usually
#' the molecule class, like "active/non-active" or a property of the molecule,
#' like "biological activity", is also stored in the SDF file. These values
#' often serve as targets for a prediction task. This function is a small
#' wrapper that extracts the information.
#' 
#' @usage 
#' getMoleculePropertyFromSDF(sdfile,property)
#' 
#' 
#' @param sdfile A character containing the name of the SDF file.
#' @param property The name of the slot in the SDF. 
#' @return A character vector with one value per molecule.
#' @author Guenter Klambauer <rchemcpp@@bioinf.jku.at>
#'
#'  @examples 
#' sdfolder <- system.file("sample_data",package="Rchemcpp")
#' sdf <- list.files(sdfolder,full.names=TRUE,pattern="small")
#' moleculeNames <- getMoleculePropertyFromSDF(sdf,"Activity")
#' 
#' @export

getMoleculePropertyFromSDF <- function(sdfile,property){
	if (missing(sdfile)){
		stop("You must submit an SDF file containing the molecules.")
	} 
	if (missing(property)){
		stop("You must submit the name of the property to be extracted.")
	}
	rl <- readLines(sdfile)
	property <- rl[grep(paste("<",property,">",sep=""),rl)+1]
	#check if number of molecules is ok
	if (length(grep("\\$\\$\\$\\$",rl))!=length(property)){
		stop("Inconsistency in SDF/Mol file. Could not read properties.")
	}
	
	return(property)
}

