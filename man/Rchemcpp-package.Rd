\name{Rchemcpp-package}
\alias{Rchemcpp-package}
\alias{Rchemcpp}
\docType{package}
\title{Rchemcpp provides tools for comparing chemical compounds}
\description{
Compares sets of chemical compounds given as SD/SDF/MOL- or KCF-files
and returns pairwise similarities as a matrix (gram matrix).
It uses the compiled-in c++ library "chemcpp" to emulate the five chemcpp tools
"sd2gram", "sd2gram3Dspectrum", "sd2gramSubtree", "sd2gram3Dpharma" and 
"sd2gramSpectrum". The tools are made accessible as R functions.
}
\details{
\tabular{ll}{
	Package: \tab Rchemcpp\cr
	Type: \tab Package\cr
	Version: \tab 1.0.2\cr
	Date: \tab 2012-10-15\cr
	License: \tab GPL2.1\cr
	}
}
\author{
Michael Mahr and Guenter Klambauer
}
\references{
(Kashima, 2004) -- H. Kashima, K. Tsuda, and A. Inokuchi. Kernels for graphs. 
In B. Schoelkopf, K. Tsuda, and J.P.
Vert, editors, Kernel Methods in Computational Biology, pages 155-170.
 MIT Press, 2004.

(Mahe, 2005) -- P. Mahe, N. Ueda, T. Akutsu, J.-L. Perret, and J.-P. Vert. 
Graph kernels for molecular structure-
activity relationship analysis with support vector machines. 
J Chem Inf Model, 45(4):939-51, 2005.

(Ralaivola, 2005) -- L. Ralaivola, S. J. Swamidass, H. Saigo, and P. Baldi. G
raph kernels for chemical informatics.
Neural Netw., 18(8):1093-1110, Sep 2005.

(Gaertner, 2003) -- T. Gaertner, P. Flach, and S. Wrobel. On graph kernels: 
hardness results and efficient alternatives.
In B. Schoelkopf and M. Warmuth, editors, Proceedings of the Sixteenth Annual 
Conference on Computational Learning Theory and the Seventh Annual Workshop on 
Kernel Machines, volume 2777 of Lecture Notes in Computer Science,
 pages 129-143, Heidelberg, 2003.

(Mahe, 2006a) -- P. Mahe and J.-P. Vert. Graph kernels based on tree patterns 
for molecules. Technical Report, HAL:ccsd-00095488, Ecoles des Mines de Paris,
September 2006.

(Mahe, 2006b) --  P. Mahe, L. Ralaivola, V. Stoven, and J.-P. Vert.
The pharmacophore kernel for virtual screening
with support vector machines. Technical Report, HAL:ccsd-00020066, Ecole des
Mines de Paris, March 2006.

(Leslie, 2002) -- C. Leslie, E. Eskin, and W.S. Noble. The spectrum kernel: 
a string kernel for SVM protein clas-
sification. In Russ B. Altman, A. Keith Dunker, Lawrence Hunter, 
Kevin Lauerdale, and Teri E.
Klein, editors, Proceedings of the Pacific Symposium on Biocomputing 2002, 
pages 564-575. World Scientific, 2002.

(Ramon, 2003) -- J. Ramon and T. Gaertner. Expressivity versus efficiency of 
graph kernels. In T. Washio and L. De Raedt, editors, 
Proceedings of the First International Workshop on Mining Graphs, 
Trees and Sequences, pages 65-74, 2003.
}
\keyword{ package }
\seealso{
\code{\link{sd2gram}}
\code{\link{sd2gram3Dpharma}}
\code{\link{sd2gramSpectrum}}
\code{\link{sd2gram3Dspectrum}}
\code{\link{sd2gramSubtree}}
}
\examples{
sdfolder <- system.file("sample_data",package="Rchemcpp")

sdf_small <- list.files(sdfolder,full.names=TRUE,pattern="small")
K1 <- sd2gram(sdf_small)
K2 <- sd2gramSpectrum(sdf_small)
K3 <- sd2gramSubtree(sdf_small)

sdf_tiny <- list.files(sdfolder,full.names=TRUE,pattern="tiny")
K4 <- sd2gram3Dspectrum(sdf_tiny)
K5 <- sd2gram3Dpharma(sdf_tiny)
}

