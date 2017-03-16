StromaExp <- function(exp_T,exp_N,exp_G, rho)
{
  exp_S=t((t(exp_T)-t(exp_G)*rho["G",]-t(exp_N)*rho["N",])/rho["S",])
  S_small=log(quantile(exp_T,0.00001))
  exp_S[exp_S<exp(S_small)]=exp(S_small)
  return(exp_S)
}
