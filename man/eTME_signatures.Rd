\name{eTME_signatures}
\alias{eTME_signatures}
\docType{data}
\title{Empirically-defined Immune Signature Genes in RCC Bulk Tumor
}
\description{Based on DisHet analysis of 35 RCC trio RNA-Seq data, we defined immune-specific genes with empirical evidence, named eTME, for empirically-defined immune signature of tumor. Using eTME, we refined previously published Immunome signatures. We also assigned other eTME genes to specific immune cell types using the BLUEPRINT project data. These two sets of refined gene signatures were consolidated and documented in the DisHet R package as the "eTME" immune cell gene signatures.
}
\usage{data("eTME_signatures")}
\format{
  A List with 26 items/vectors.
  \describe{
    \item{'M2 macrophages'}{:   a vector of genes that are abundantly expressed in M2 macrophages}
    \item{'M1 macrophages'}{:   a vector of genes that are abundantly expressed in M1 macrophages}
    \item{Macrophages}{:   a vector of genes that are abundantly expressed in Macrophages}
    \item{Monocytes}{:   a vector of genes that are abundantly expressed in Monocytes}
    \item{Thymocytes}{:   a vector of genes that are abundantly expressed in Thymocytes}
    \item{'Megakaryocyte cells'}{:   a vector of genes that are abundantly expressed in Megakaryocyte cells}
    \item{'B cells'}{:   a vector of genes that are abundantly expressed in B cells}
    \item{'CD8 T cells'}{:   a vector of genes that are abundantly expressed in CD8 T cells}
    \item{'T cells'}{:   a vector of genes that are abundantly expressed in T cells}
    \item{'Common lymphoid progenitors'}{:   a vector of genes that are abundantly expressed in Common lymphoid progenitors}
    \item{'Dendritic cells'}{:   a vector of genes that are abundantly expressed in Dendritic cells}
    \item{'CD56dim NK cells'}{:   a vector of genes that are abundantly expressed in CD56dim NK cells}
    \item{'CD56bright NK cells'}{:   a vector of genes that are abundantly expressed in CD56bright NK cells}
    \item{'NK cells'}{:   a vector of genes that are abundantly expressed in NK cells}
    \item{'Endothelial cells'}{:   a vector of genes that are abundantly expressed in Endothelial cells}
    \item{Erythroblasts}{:   a vector of genes that are abundantly expressed in Erythroblasts}
    \item{Eosinophils}{:   a vector of genes that are abundantly expressed in Eosinophils}
    \item{Neutrophils}{:   a vector of genes that are abundantly expressed in Neutrophils}
    \item{'Treg cells'}{:   a vector of genes that are abundantly expressed in Treg cells}
    \item{'Th1 cells'}{:   a vector of genes that are abundantly expressed in Th1 cells}
    \item{'Th2 cells'}{:   a vector of genes that are abundantly expressed in Th2 cells}
    \item{'Th cells'}{:   a vector of genes that are abundantly expressed in Th cells}
    \item{aDCs}{:   a vector of genes that are abundantly expressed in aDCs}
    \item{iDCs}{:   a vector of genes that are abundantly expressed in iDCs}
    \item{pDCs}{:   a vector of genes that are abundantly expressed in pDCs}
    \item{'Mast cells'}{:   a vector of genes that are abundantly expressed in Mast cells}
  }
}
%% \details{
%%  ~~ If necessary, more details than the __description__ above ~~
%% }
%% \source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%% Manuscript "An Empirical Approach Leveraging Tumorgrafts to Dissect the Tumor Microenvironment in Renal Cell Carcinoma Identifies Missing Link to Prognostic Inflammatory Factors" (recently submitted to Nature Genetics). 
%% }
%%\references{
%%  ~~ possibly secondary sources and usages ~~
%%}
\examples{
data(eTME_signatures)
eTME_signatures$Macrophages
eTME_signatures$'T cells'
}
\keyword{datasets}
