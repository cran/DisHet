\name{StromaExp}
\alias{StromaExp}
\title{Stroma (Immune) Component Gene Expression Estimation}
\description{
  This function estimate the stroma component gene expression profiles of all patients, using the proportion estimates obtained from function \code{DisHet}. The estimates are stored in a p-by-k matrix, where p is the number of genes and k is the number of patients.The order of genes and the order of patients are the same as in the input bulk sample expression matrix.  
}
\usage{
StromaExp(exp_T,exp_N,exp_G, rho)
}
\arguments{
  \item{exp_T}{
    Gene expression in bulk RNA-seq samples. The rows correspond to different genes. The columns correspond to different patients.
  }
  \item{exp_N}{
    Gene expression in the corresponding normal samples. The rows list the same set of genes as in exp_G. The columns correspond to patients matched with exp_T.
  }
  \item{exp_G}{
    Gene expression in the corresponding tumor samples. The rows list the same set of genes as in exp_G. The columns correspond to patients matched with exp_T.
  }
  \item{rho}{
    Output from function \code{DisHet}: the patient-specific proportion estimates corresponding to tumor, normal, stroma components in order.  
  }
}
\examples{
  load(system.file("example/example_data.RData",package="DisHet"))
  exp_T <- exp_T[1:200,]
  exp_N <- exp_N[1:200,]
  exp_G <- exp_G[1:200,]
  
  rho <- DisHet(exp_T, exp_N, exp_G, save=FALSE, n_cycle=500, mean_last=50)
  S <- StromaExp(exp_T,exp_N,exp_G, rho)
}