\name{DisHet}
\alias{DisHet}
\title{Heterogeneity Dissection}
\description{This function performs dissection of bulk sample gene expression using matched normal and tumorgraft RNA-seq data. It outputs the final proportion estiamtes of the three components for all patients. 
  
  The patient-specific dissection proportion estimates are saved in a 3-by-k matrix named "rho", where k is the number of patients. The 3 rows of "rho" matrix correspond to the tumor, normal, stroma components in order. That is, the proportion estimate of tumor component for patient i is stored in rho[1,i]; the normal component proportion estimate of this patient is stored in rho[2,i];and stroma component proportion in rho[3,i]. 
}
\usage{
DisHet(exp_T,exp_N,exp_G, save=TRUE, MCMC_folder, 
      n_cycle=10000, save_last=500, mean_last=200, dirichlet_c=1, S_c=1, rho_small=1e-2, 
      initial_rho_S=0.02,initial_rho_G=0.96,initial_rho_N=0.02)
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
  \item{save}{
    When save==TRUE, as in default, all component proportion estimates during MCMC iterations can be saved into a user-specified directory using the "MCMC_folder" argument. 
  }
  \item{MCMC_folder}{
    Directory for saving the estimated mixture proportion matrix updates during MCMC iterations. The default setting is to create a "DisHet" folder under the current working directory. 
  }
  \item{n_cycle}{
    Number of MCMC iterations(chain length). The default value is 10,000. 
  }
  \item{save_last}{
    Save the rho matrix updates for the last "save_last" Number of MCMC iterations. The default value is 500. 
  }
  \item{mean_last}{
    Calculate the final proportion estiamte matrix using the last "mean_last" number of MCMC iterations. The default value is 200.  
  }
  \item{dirichlet_c}{
    Stride scale in sampling rho. Larger value leads to smaller steps in sampling rho. The default value is 1. 
  }
  \item{S_c}{
    Stride scale in sampling Sij. Larger value leads to larger steps in sampling Sij. The default value is 1. 
  }
  \item{rho_small}{
    The smallest rho updates allowed during MCMC. The default is 1e-2. This threshold is set to help improve numerical stability of the algorithm. 
  }
  \item{initial_rho_S}{
    Initial value of the proportion estimate for the stroma component. The default value is 0.02. 
  }
  \item{initial_rho_G}{
    Initial value of the proportion estimate for the tumor component. The default value is 0.96. 
  }
  \item{initial_rho_N}{
    Initial value of the proportion estimate for the normal component. The default value is 0.02. 
  }
}
\details{
  Un-logged expression values should be used in exp_N/T/G matrices, and their rows and columns must match each other corresponding to the same set of genes and patients. 
  
  The values specified for "initial_rho_S", "initial_rho_G", and "initial_rho_S" all have to be positive. If the three proportion initials are not summing to 1, normalization is performed automatically to force the sum to be 1. 
}
%% \references{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%% Manuscript "An Empirical Approach Leveraging Tumorgrafts to Dissect the Tumor Microenvironment in Renal Cell Carcinoma Identifies Missing Link to Prognostic Inflammatory Factors" (recently submitted to Nature Genetics). 
%% }
\examples{
  load(system.file("example/example_data.RData",package="DisHet"))
  exp_T <- exp_T[1:200,]
  exp_N <- exp_N[1:200,]
  exp_G <- exp_G[1:200,]
  
  rho <- DisHet(exp_T, exp_N, exp_G, save=FALSE, n_cycle=200, mean_last=50)
}

