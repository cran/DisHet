DisHet <- function(exp_T,exp_N,exp_G, save=TRUE, MCMC_folder=NULL, 
      n_cycle=10000, save_last=500, mean_last=200, dirichlet_c=1, S_c=1, rho_small=1e-2, 
      initial_rho_S=0.02,initial_rho_G=0.96,initial_rho_N=0.02)
{
  if(is.null(MCMC_folder)){MCMC_folder<- paste(getwd(),"/DisHet",sep="")}
  if(save==TRUE & !dir.exists(MCMC_folder)){dir.create(MCMC_folder)}
  
  rho=matrix(0,nrow=3,ncol=dim(exp_G)[2],dimnames=list(c("G","N","S"),colnames(exp_G)))
  rho_final=rho
  if(initial_rho_S>0 & initial_rho_G>0 & initial_rho_N>0){iniSum <- initial_rho_S+initial_rho_G+initial_rho_N}
  else{print("Initialization Error: initial values have to be positive.");stop();}
  if(iniSum!=1){
    print("Initialization Warning: sum of the initial proportions is not 1. Normalization of the specified initials is performed.");
    initial_rho_S <- initial_rho_S/iniSum; 
    initial_rho_G <- initial_rho_G/iniSum;
    initial_rho_N <- initial_rho_N/iniSum;}
  rho["S",]=initial_rho_S # initial state
  rho["G",]=initial_rho_G
  rho["N",]=initial_rho_N
  
  I=dim(exp_G)[2] ## number of patients
  J=dim(exp_G)[1] ## number of genes
  S_small=log(quantile(exp_T,0.00001)) # smallest possible Sij values
  S=apply(log(exp_T),1,mean)
  mu=mean(S)
  tao2=var(S)
  delta2=apply(log(exp_T),1,var)/100
  
  # Estimating an Inverse Gamma distribution, A. Llera, C. F. Beckmann
  alpha=mean(delta2)^2/var(delta2)+2
  beta=mean(delta2)*(alpha-1)
  alphas=c("G"=0.978,"N"=0.011,"S"=0.011)*300
  
  ############  MCMC  ##################
  
  log_exp_T=log(exp_T)
  t_exp_G=t(exp_G)
  t_exp_N=t(exp_N)
  
  for (cycle in 1:n_cycle)
  {
    ##########  update S  ############################
    S_new=rnorm(length(S),S,sqrt(delta2)*S_c)
    S_new[S_new<S_small]=S_small
    
    X=(log_exp_T-log(t(t_exp_G*rho["G",]+t_exp_N*rho["N",])+
                       exp(S) %*% rho["S",,drop=F]))^2
    X_new=(log_exp_T-log(t(t_exp_G*rho["G",]+t_exp_N*rho["N",])+
                           exp(S_new) %*% rho["S",,drop=F]))^2
    
    ratio=exp(-rowSums(X_new-X)/2/delta2-((S_new-mu)^2-(S-mu)^2)/2/tao2)
    ratio[ratio>1]=1
    cat(paste("  Overall acceptance rate for S:",round(mean(ratio),digits=2),"\n"))
    
    keep=runif(J)<ratio
    S[keep]=S_new[keep]
    X[keep,]=X_new[keep,]
    
    ##########  update delta2  ########################
    
    scale=rowSums(X)/2+beta
    for (j in 1:J) {delta2[j]=1/rgamma(1, shape=alpha+I/2,scale=1/scale[j])}
    
    
    ##########  update rho  ###########################
    rho_new=rho # propose new rho 
    d_rho_new=d_rho=rep(0,I)
    
    for (i in 1:I) 
    {
      rho_new[,i]=rdirichlet(1,rho[,i]*dirichlet_c)[1,]
      rho_new[rho_new[,i]<rho_small,i]=rho_small
      rho_new[,i]=rho_new[,i]/sum(rho_new[,i])
      d_rho_new[i]=ddirichlet(rho_new[,i],rho[,i]*dirichlet_c) # given rho, density of rho_new
      d_rho[i]=ddirichlet(rho[,i],rho_new[,i]*dirichlet_c) # given rho_new, density of rho
    }
    
    X_new=(log_exp_T-log(t(t_exp_G*rho_new["G",]+t_exp_N*rho_new["N",])+
                           exp(S) %*% rho_new["S",,drop=F]))^2
    
    # calculate probabilities
    ratio=exp(colSums(-(X_new-X)/delta2)/2+ # ratio of posterior probability
                log((rho_new["G",]/rho["G",]))*(alphas["G"]-1)+
                log((rho_new["N",]/rho["N",]))*(alphas["N"]-1)+
                log((rho_new["S",]/rho["S",]))*(alphas["S"]-1)+
                log(d_rho/d_rho_new)) # ratio of jumping probability
    ratio[ratio>1]=1
    cat(paste("  Overall acceptance rate for rho:",round(mean(ratio),digits=2),"\n"))
    
    # accept
    keep=runif(I)<=ratio
    rho[,keep]=rho_new[,keep]
    
    ############  save rho estimates  ##########
    cat(paste("Finished iteration",cycle,"\n\n\n"))
    if(save==TRUE & (cycle>(n_cycle-save_last))){
    save_file=paste(MCMC_folder,"/",cycle,".RData",sep="")
    save(rho,file=save_file)  }
    
    if(cycle>(n_cycle-mean_last)){ rho_final = rho_final + rho }
  }#end of MCMC
  
  rho_final = rho_final/mean_last
  return(rho_final)
}
