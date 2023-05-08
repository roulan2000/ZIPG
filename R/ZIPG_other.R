#' ZIPG EM algorithm to optimize log-likelihood
#'
#' @param W Count data
#' @param M Sequencing depth, ZIPG use log(M) as offset by default
#' @param X Formula of covariates of differential abundance
#' @param X_star Formula of covariates of differential varibility
#' @param A (No use, set A=1)
#' @param d The number of covariates of differential abundance
#' @param d_star The number of covariates of differential varibility
#' @param parms initialization parameters
#' @param optim_method Default : 'BFGS'
#' @param fix_index Index of parameters that are fixed in optimization
#' @param tol convergence criterion
#' @param maxit maximum number of iteration
#' @noRd

ZIPG_optim_EM <- function(W,M,X,X_star,
                        A=1,d,d_star,
                        parms,optim_method ='BFGS',
                        fix_index=NULL,tol = 1e-5,maxit = 200)
{
  grad_ZIPG <- function(par,fix_par,fix_index=NULL) {
    if(!is.null(fix_index))
    {
      full_index = c(1:(length(par)+length(fix_par)))
      opt_index = setdiff(full_index,fix_index)
      parms = rep(0,length(full_index))
      parms[fix_index] = fix_par
      parms[opt_index] = par
    }else{
      parms = par
    }
    length_parms = length(parms)
    n = length(M)
    beta = matrix(parms[1:(A*(d+1))],d+1,A)
    beta_star = matrix(0,d_star+1,A)
    beta_star[1,] = parms[(A*(d+1)+1) : (A*(d+2))]
    if(d_star!=0){
      beta_star[2:(d_star+1),]  =  parms[(A*(d+2)+1) : (A*(d+2)+d_star)]}
    gamma_0 = t(as.matrix(parms[(A*(d+2)+d_star+1) : (A*(d+3)+d_star)]))
    p_0 = 1/(1+exp(-gamma_0))
    lambda= matrix(0,n,A)
    theta= matrix(0,n,A)

    nb0 = matrix(0,n,A)

    grad_beta = matrix(0,d+1,A)

    for(k in 1:A)
    {
      lambda[,k] = exp(X %*% beta[,k])*M
      theta[,k] = exp(X_star %*% beta_star[,k])
      nb0[,k] = exp(gamma_0[k])*(1+lambda[,k]*theta[,k])^(1/theta[,k])
    }

    phi_lambda =
      ifelse(W!=0,W- (1+W*theta)/(1+lambda*theta)*lambda,
             -lambda/(1+nb0)/(1+lambda*theta))


    phi_theta =(-digamma(W + 1/theta)/theta + digamma(1/theta)/theta+
                  log(1+lambda*theta)/theta + (W-lambda)/(1+lambda*theta))*
      ifelse(W!=0, 1, 1/(1+nb0))


    phi_gamma_0= -matrix(rep(p_0,each = n),n,A) +
      ifelse(W!=0, 0, 1/(1+1/nb0))


    for(k in 1:A){
      phi_beta = phi_lambda[,k]*X
      grad_beta[,k] = colSums(phi_beta)
    }

    phi_beta_star = t(phi_theta) %*% X_star
    grad_beta_star_0 = phi_beta_star[,1]
    grad_beta_star = colSums(phi_beta_star)[-1]

    grad_gamma_0 = colSums(phi_gamma_0)

    grad = c(c(grad_beta),grad_beta_star_0,
             grad_beta_star,grad_gamma_0)
    if(!is.null(fix_index))
    {
      grad = grad[-fix_index]

    }
    return(grad)
  }


  Estep_zp <- function(parms,init = FALSE)
  {


    length_parms = length(parms)
    n = length(M)

    beta = matrix(parms[1:(A*(d+1))],d+1,A)
    beta_star = matrix(0,d_star+1,A)
    beta_star[1,] = parms[(A*(d+1)+1) : (A*(d+2))]
    if(d_star!=0){
      beta_star[2:(d_star+1),]  =  parms[(A*(d+2)+1) : (A*(d+2)+d_star)]}
    gamma_0 = t(as.matrix(parms[(A*(d+2)+d_star+1) : (A*(d+3)+d_star)]))
    p_0 = 1/(1+exp(-gamma_0))
    if(init == TRUE){
      p_0_col = matrix(rep(p_0,n),n,A,byrow = TRUE)
      return(p_0_col)
    }else{

      lambda= matrix(0,n,A)
      theta= matrix(0,n,A)

      zp = matrix(0,n,A)
      nb_prob = matrix(0,n,A)
      zp0 = matrix(0,n,A)
      for(k in 1:A)
      {

        lambda[,k] = exp(X %*% beta[,k])*M
        theta[,k] = exp(X_star %*% beta_star[,k])

        nb_prob[,k] = stats::dnbinom(W[,k],
                              mu = lambda[,k],size = 1/theta[,k],
                              log = TRUE)

        zp0[,k] = p_0[,k]/(p_0[,k]+exp(log(1-p_0[,k])+nb_prob[,k]))
        zp[,k] = ifelse(W[,k]==0,zp0,0)
        zp[which(is.na(zp[,k])),k]= 0
      }
      return(zp)
    }
  }
  Mstep_zp <- function(parms,fix_index,zp){
    likelihood_ZIPG_EM <- function(par,fix_par,fix_index=NULL,zp) {
      if(!is.null(fix_index))
      {
        full_index = c(1:(length(par)+length(fix_par)))
        opt_index = setdiff(full_index,fix_index)
        parms = rep(0,length(full_index))
        parms[fix_index] = fix_par
        parms[opt_index] = par
      }else{
        parms = par
      }

      length_parms = length(parms)
      n = length(M)



      beta = matrix(parms[1:(A*(d+1))],d+1,A)
      beta_star = matrix(0,d_star+1,A)
      beta_star[1,] = parms[(A*(d+1)+1) : (A*(d+2))]
      if(d_star!=0){
        beta_star[2:(d_star+1),]  =  parms[(A*(d+2)+1) : (A*(d+2)+d_star)]}
      gamma_0 = t(as.matrix(parms[(A*(d+2)+d_star+1) : (A*(d+3)+d_star)]))
      p_0 = 1/(1+exp(-gamma_0))

      lambda= matrix(0,n,A)
      theta= matrix(0,n,A)

      logli_per_W = matrix(0,n,A)
      nb_prob = matrix(0,n,A)
      zp1= matrix(0,n,A)
      zp0 = matrix(0,n,A)
      for(k in 1:A)
      {

        lambda[,k] = exp(X %*% beta[,k])*M
        theta[,k] = exp(X_star %*% beta_star[,k])

        nb_prob[,k] = stats::dnbinom(W[,k],
                              mu = lambda[,k],size = 1/theta[,k],
                              log = TRUE)


        logli_per_W[,k] = zp[,k]*log(p_0[,k])+(1-zp[,k])*
          (log(1-p_0[,k])+nb_prob[,k])
      }
      return(sum(sum(logli_per_W)))

    }

    grad_ZIPG_EM <- function(par,fix_par,fix_index=NULL,zp) {



      if(!is.null(fix_index))
      {
        full_index = c(1:(length(par)+length(fix_par)))
        opt_index = setdiff(full_index,fix_index)
        parms = rep(0,length(full_index))
        parms[fix_index] = fix_par
        parms[opt_index] = par
      }else{
        parms = par
      }
      length_parms = length(parms)
      n = length(M)
      beta = matrix(parms[1:(A*(d+1))],d+1,A)
      beta_star = matrix(0,d_star+1,A)
      beta_star[1,] = parms[(A*(d+1)+1) : (A*(d+2))]
      if(d_star!=0){
        beta_star[2:(d_star+1),]  =  parms[(A*(d+2)+1) : (A*(d+2)+d_star)]}
      gamma_0 = t(as.matrix(parms[(A*(d+2)+d_star+1) : (A*(d+3)+d_star)]))
      gamma_0_col =matrix( rep(gamma_0,n),n,A)
      p_0 = 1/(1+exp(-gamma_0))

      lambda= matrix(0,n,A)
      theta= matrix(0,n,A)

      nb0 = matrix(0,n,A)

      grad_beta = matrix(0,d+1,A)


      for(k in 1:A)
      {
        lambda[,k] = exp(X %*% beta[,k])*M
        theta[,k] = exp(X_star %*% beta_star[,k])

      }

      phi_lambda =(W- (1+W*theta)/(1+lambda*theta)*lambda)*(1-zp)


      phi_theta =(-digamma(W + 1/theta)/theta + digamma(1/theta)/theta+
                    log(1+lambda*theta)/theta + (W-lambda)/(1+lambda*theta))*
        (1-zp)


      phi_gamma_0= zp-1/(1+exp(-gamma_0_col))


      for(k in 1:A){
        phi_beta = phi_lambda[,k]*X
        grad_beta[,k] = colSums(phi_beta)
      }

      phi_beta_star = t(phi_theta) %*% X_star
      grad_beta_star_0 = phi_beta_star[,1]
      grad_beta_star = colSums(phi_beta_star)[-1]

      grad_gamma_0 = colSums(phi_gamma_0)

      grad = c(c(grad_beta),grad_beta_star_0,
               grad_beta_star,grad_gamma_0)
      if(!is.null(fix_index))
      {
        grad = grad[-fix_index]

      }
      return(grad)
    }

    if(A==1){

      if(is.null(fix_index))
      {
        par = parms
        fix_par = NULL
      } else{
        par = parms[-fix_index]
        fix_par = parms[fix_index]
      }

      if(optim_method == 'L-BFGS-B'){
        fit <- tryCatch(stats::optim(fn = likelihood_ZIPG_EM, gr = grad_ZIPG_EM,
                              par = par,
                              method = "L-BFGS-B", hessian = TRUE,
                              control = list(maxit = 10000,trace = 0,
                                             fnscale = -1,factr = 1e-6),
                              fix_index = fix_index,fix_par = fix_par,zp=zp),
                        error = function(e){
                          warning("optim ERROR :",conditionMessage(e),"\n")
                          return(conditionMessage(e))
                        })
      } else {
        fit <- tryCatch(stats::optim(fn = likelihood_ZIPG_EM, gr = grad_ZIPG_EM,
                              par = par,
                              method = "BFGS", hessian = TRUE,
                              control = list(maxit = 10000,trace = 0,
                                             fnscale = -1,reltol = 1e-6),
                              fix_index = fix_index,fix_par = fix_par,zp=zp),
                        error = function(e){
                          warning("optim ERROR :",conditionMessage(e),"\n")
                          return(conditionMessage(e))
                        })

        if(length(fit)>=2){
          fit$gardient = grad_ZIPG(fit$par,fix_index = fix_index,fix_par = fix_par)}
      }
      return(fit)
    }else{
      length_parms = length(parms)
      parms_index = c(1:length_parms)
      if(d_star!=0){
        beta_star_index  =  c((A*(d+2)+1) : (A*(d+2)+d_star))}
      taxa_index = matrix(0,d+3,A)
      for(t in 1:A){
        taxa_index[,t] = c(
          ((t-1)*(d+1)+1):((t-1)*(d+1)+d+1),
          A*(d+1)+t,
          A*(d+2)+d_star+t
        )
      }
      parms_fit_single = parms
      if(1){
      for(t in 1:A){
        fix_index = parms_index[-taxa_index[,t]]
        par = parms_fit_single[-fix_index]
        fix_par = parms_fit_single[fix_index]
        fit <- tryCatch(stats::optim(fn = likelihood_ZIPG_EM, gr = grad_ZIPG_EM,
                              par = par,
                              method = "BFGS", hessian = TRUE,
                              control = list(maxit = 10000,trace = 0,
                                             fnscale = -1,reltol = 1e-6),
                              fix_index = fix_index,fix_par = fix_par,zp=zp),
                        error = function(e){
                          warning("optim ERROR :",conditionMessage(e),"\n")
                          return(conditionMessage(e))
                        })

        if(length(fit)>=2){

          parms_fit_single[taxa_index[,t]] = fit$par
        }
      }

        if(0){
          fix_index = parms_index[-beta_star_index]
          par = parms_fit_single[-fix_index]
          fix_par = parms_fit_single[fix_index]
          fit <- tryCatch(stats::optim(fn = likelihood_ZIPG_EM, gr = grad_ZIPG_EM,
                                par = par,
                                method = "BFGS", hessian = TRUE,
                                control = list(maxit = 10000,trace = 0,
                                               fnscale = -1,reltol = 1e-6),
                                fix_index = fix_index,fix_par = fix_par,zp=zp),
                          error = function(e){
                            warning("optim ERROR :",conditionMessage(e),"\n")
                            return(conditionMessage(e))
                          })

          if(length(fit)>=2){
            fit$gardient = grad_ZIPG(fit$par,fix_index = fix_index,fix_par = fix_par)
            parms_fit_single[beta_star_index] = fit$par
            fit$hessian = diag(c(fit$hessian),length_parms)
          }
        }
      }
      fit_final <- tryCatch(stats::optim(fn = likelihood_ZIPG_EM, gr = grad_ZIPG_EM,
                            par = parms_fit_single,
                            method = "BFGS", hessian = TRUE,
                            control = list(maxit = 1,trace = 0,
                                           fnscale = -1,reltol = 1e-6),
                            fix_index = NULL,fix_par = NULL,zp=zp),
                      error = function(e){
                        warning("optim ERROR :",conditionMessage(e),"\n")
                        return(conditionMessage(e))
                      })


      return(fit_final)
    }
  }

  likelihood_ZIPG_EM_zp <- function(parms,zp) {
    length_parms = length(parms)
    n = length(M)



    beta = matrix(parms[1:(A*(d+1))],d+1,A)
    beta_star = matrix(0,d_star+1,A)
    beta_star[1,] = parms[(A*(d+1)+1) : (A*(d+2))]
    if(d_star!=0){
      beta_star[2:(d_star+1),]  =  parms[(A*(d+2)+1) : (A*(d+2)+d_star)]}
    gamma_0 = t(as.matrix(parms[(A*(d+2)+d_star+1) : (A*(d+3)+d_star)]))
    p_0 = 1/(1+exp(-gamma_0))

    lambda= matrix(0,n,A)
    theta= matrix(0,n,A)

    logli_per_W = matrix(0,n,A)
    nb_prob = matrix(0,n,A)
    zp1= matrix(0,n,A)
    zp0 = matrix(0,n,A)
    for(k in 1:A)
    {

      lambda[,k] = exp(X %*% beta[,k])*M
      theta[,k] = exp(X_star %*% beta_star[,k])

      nb_prob[,k] = stats::dnbinom(W[,k],
                            mu = lambda[,k],size = 1/theta[,k],
                            log = TRUE)


      zp0 <- log(p_0[,k] + exp(log(1 - p_0[,k]) +nb_prob[,k]))
      zp1 <- log(1 - p_0[,k]) + nb_prob[,k]
      logli_per_W[,k] = zp[,k]*log(p_0[,k])+(1-zp[,k])*
        (log(1-p_0[,k])+nb_prob[,k])
    }
    return(sum(sum(logli_per_W)))

  }

  parms0 = parms
  zp0 = Estep_zp(parms0,FALSE)
  Q0 = likelihood_ZIPG_EM_zp(parms0,zp0)
  p_zp0 = colMeans(zp0,na.rm = TRUE)
  gamma_zp0 = log(p_zp0/(1-p_zp0))

  nIter = 0
  while(TRUE){
    nIter = nIter + 1
    parms_it = parms
    if(d_star!=0 && A >3){
      beta_star_index  =  c((A*(d+2)+1) : (A*(d+2)+d_star))
      parms_it[beta_star_index] = parms0[beta_star_index]
    }

    fit1 = Mstep_zp(parms_it,fix_index,zp0)
    if(length(fit1)==1){
      return(fit1)
    }
    Q1 = fit1$value
    parms1 = fit1$par
    zp1 = Estep_zp(parms1)
    if(abs(Q1-Q0)/abs(Q0) < tol | nIter>=maxit) break
    Q0 = Q1
    parms0 =parms1
    zp0 = zp1
    p_zp1 = colMeans(zp0,na.rm = TRUE)
    gamma_zp1 = log(p_zp0/(1-p_zp0))
  }
  fit1$EM_nIter = nIter
  fit1$EM_converge = !(nIter==maxit)

  fit1$init = parms
  return(fit1)
}

#' Log-likelihood for ZIPG
#'
#' @param W Count data
#' @param M Sequencing depth, ZIPG use log(M) as offset by default
#' @param X Formula of covariates of differential abundance
#' @param X_star Formula of covariates of differential varibility
#' @param A (No use, set A=1)
#' @param d The number of covariates of differential abundance
#' @param d_star The number of covariates of differential varibility
#' @param parms Vector of parameters including beta,beta_star and gamma.
#' @param sep (No use, set FALSE)Return log-likelihodd of each sample separately.
#'
#' @return The value of log-likelihood of observed data under ZIPG model.
#' @noRd
ZIPG_logli<- function(W,M,X,X_star,
                    A=1,d,d_star,
                    parms,sep=FALSE)
{
  W = as.matrix(W)
  length_parms = length(parms)
  n = length(M)
  beta = matrix(parms[1:(A*(d+1))],d+1,A)
  beta_star = matrix(0,d_star+1,A)
  beta_star[1,] = parms[(A*(d+1)+1) : (A*(d+2))]
  if(d_star!=0){
    beta_star[2:(d_star+1),]  =  parms[(A*(d+2)+1) : (A*(d+2)+d_star)]}
  gamma_0 = t(as.matrix(parms[(A*(d+2)+d_star+1) : (A*(d+3)+d_star)]))
  p_0 = 1/(1+exp(-gamma_0))

  lambda= matrix(0,n,A)
  theta= matrix(0,n,A)

  logli_per_W = matrix(0,n,A)
  logli_per_A = rep(0,A)
  nb_prob = matrix(0,n,A)
  zinb_1 = matrix(0,n,A)
  zinb_0 = matrix(0,n,A)
  for(k in 1:A)
  {

    lambda[,k] = exp(X %*% beta[,k])*M
    theta[,k] = exp(X_star %*% beta_star[,k])

    nb_prob[,k] = stats::dnbinom(W[,k],
                          mu = lambda[,k],size = 1/theta[,k],
                          log = TRUE)

    zinb0 <- log(p_0[,k] + exp( log(1 - p_0[,k]) + nb_prob[,k] ))
    zinb1 <- log(1 - p_0[,k]) + nb_prob[,k]
    logli_per_W[,k] = ifelse(W[,k]==0,zinb0,zinb1)
    logli_per_A[k] = sum(logli_per_W[,k])
  }
  if(sep == TRUE)
    return(logli_per_A)
  else
    return(sum(logli_per_A))

}

#' Log-likelihood for ZIPG given latent variable
#'
#' @param W Count data
#' @param M Sequencing depth, ZIPG use log(M) as offset by default
#' @param X Formula of covariates of differential abundance
#' @param X_star Formula of covariates of differential varibility
#' @param A (No use, set A=1)
#' @param d The number of covariates of differential abundance
#' @param d_star The number of covariates of differential varibility
#' @param parms Vector of parameters including beta,beta_star and gamma.
#' @param zp Latent variable of zero-inflation
#'
#' @return The value of log-likelihood of complete data under ZIPG model.
#' @noRd
ZIPG_logli_EM <- function(W,M,X,X_star,
                        A,d,d_star,
                        parms,zp) {
  W = as.matrix(W)
  length_parms = length(parms)
  n = length(M)

  beta = matrix(parms[1:(A*(d+1))],d+1,A)
  beta_star = matrix(0,d_star+1,A)
  beta_star[1,] = parms[(A*(d+1)+1) : (A*(d+2))]

  if(d_star!=0){
    beta_star[2:(d_star+1),]  =  parms[(A*(d+2)+1) : (A*(d+2)+d_star)]}
  gamma_0 = t(as.matrix(parms[(A*(d+2)+d_star+1) : (A*(d+3)+d_star)]))
  p_0 = 1/(1+exp(-gamma_0))

  lambda= matrix(0,n,A)
  theta= matrix(0,n,A)

  logli_per_W = matrix(0,n,A)
  nb_prob = matrix(0,n,A)
  zp1= matrix(0,n,A)
  zp0 = matrix(0,n,A)
  for(k in 1:A)
  {

    lambda[,k] = exp(X %*% beta[,k])*M
    theta[,k] = exp(X_star %*% beta_star[,k])

    nb_prob[,k] = stats::dnbinom(W[,k],
                          mu = lambda[,k],size = 1/theta[,k],
                          log = TRUE)


    zp0 <- log(p_0[,k] + exp(log(1 - p_0[,k]) +nb_prob[,k]))
    zp1 <- log(1 - p_0[,k]) + nb_prob[,k]
    logli_per_W[,k] = zp[,k]*log(p_0[,k])+(1-zp[,k])*
      (log(1-p_0[,k])+nb_prob[,k])
  }
  return(sum(sum(logli_per_W)))

}


#' Initialize ZIPG from pscl estimation
#'
#' @param W Count data
#' @param M Sequencing depth, ZIPG use log(M) as offset by default
#' @param X Formula of covariates of differential abundance
#' @param X_star Formula of covariates of differential varibility
#' @param A (No use, set A=1)
#' @param d The number of covariates of differential abundance
#' @param d_star The number of covariates of differential varibility
#' @param return_model Whether return the full model
#'
#' @return A list of pscl result.
#' @noRd
ZIPG_init_pscl <- function(W,M,X,X_star,
                         A=1,d,d_star,
                         return_model = TRUE){
  length_parms = A*(d+3)+d_star
  beta_init = matrix(0,d+1,A)
  beta_star0_init = rep(0,A)
  beta_star_init = rep(0,d_star)
  gamma0_init = rep(0,A)

  beta_init_pval = matrix(0,d+1,A)
  beta_star0_init_pval = rep(0,A)
  beta_star_init_pval = rep(0,d_star)
  gamma0_init_pval = rep(0,A)

  beta_init_SE = matrix(0,d+1,A)
  beta_star0_init_SE = rep(0,A)
  beta_star_init_SE = rep(0,d_star)
  gamma0_init_SE = rep(0,A)

  logli = rep(0,A)
  pscl = list()
  for(k in 1:A){
    if(d==0){
      pscl[[k]] <- pscl::zeroinfl(W[,k]~offset(log(M))|1, dist = "negbin")
    }else{
      pscl[[k]] <- pscl::zeroinfl(W[,k]~offset(log(M))+X[,-1]|1, dist = "negbin")
    }
    beta_init[,k] =pscl[[k]][["coefficients"]][["count"]]
    beta_star0_init[k]=-log(pscl[[k]]$theta)
    gamma0_init[k] = pscl[[k]][["coefficients"]][["zero"]]

    ls_pscl = pval_getpval(pscl[[k]])

    beta_init_pval[,k] = ls_pscl$pval[1:(d+1)]
    beta_star0_init_pval[k] = ls_pscl$pval[d+2]
    gamma0_init_pval[k] = ls_pscl$pval[d+3]

    beta_init_SE[,k] = ls_pscl$SE[1:(d+1)]
    beta_star0_init_SE[k] = ls_pscl$SE[d+2]
    gamma0_init_SE[k] = ls_pscl$SE[d+3]

  }

  pscl_init=c(c(beta_init),beta_star0_init,beta_star_init,gamma0_init)
  pscl_init_pval = c(c(beta_init_pval),beta_star0_init_pval,
                     beta_star_init_pval,gamma0_init_pval)
  pscl_init_SE=c(c(beta_init_SE),beta_star0_init_SE,
                 beta_star_init_SE,gamma0_init_SE)

  li_pscl = ZIPG_logli(W,M,X,X_star,A,d,d_star,pscl_init,sep=TRUE)
  if(!return_model)
  {
    pscl = NULL
  }
  return(list(
    init_model = pscl[[1]],
    init_parms = pscl_init,
    init_SE = pscl_init_SE,
    init_pval = pscl_init_pval,
    init_logli = li_pscl
  ))
}

#' Get Wald test p-value
#'
#' @param res Output from ZIPG_optim_EM
#' @param test_index Index of which parameters equal to zero under H0.
#'
#' @return A list of Wald test result
#' @noRd
pwald <- function(res,test_index)
{
  solve_hessian<- function(hessian){
    svd_res = svd(hessian)
    smalls = which(svd_res$d<1e-5)
    if(length(smalls)>=1) warning('Small eigen values ',svd_res$d[smalls])
    svd_res$d[smalls] = 1e-5
    return(svd_res$v %*% diag(1/svd_res$d) %*% t(svd_res$u))
  }
  if(length(res)==1){
    warning('optim error')
    return(list(SE=NA,
                twald = NA,
                pval = NA))
  }

  n= length(res$par)
  cov_mat = tryCatch(solve(-res$hessian),
                     error = function(e){
                       return("error")
                     })
  if(length(cov_mat)==1){
    cov_mat = tryCatch(solve_hessian(-res$hessian),
                       error = function(e){
                         return("error")
                       })
  }
  if(length(cov_mat)==1){
    warning('solve hessian error')
    return(list(SE=rep(NA,n),
                twald = rep(NA,n),
                pval = rep(NA,n)))

  }else{

    twald=c(1:n)
    pval=c(1:n)
    SE = c(1:n)



    if(is.null(test_index))
    {
      for(i in 1:n){
        par = as.matrix(res$par[i])
        cov = as.matrix(cov_mat[i,i])
        SE[i] = sqrt(abs(cov))
        if(cov==0 | is.na(cov)){cov=NA}
        cov = abs(solve(cov))
        twald[i] = t(par) %*% cov %*% par
        pval[i] = stats::pchisq(twald[i],df=1,lower.tail = FALSE)
      }
    }else{
      par = as.matrix(res$par[test_index])
      cov = as.matrix(cov_mat[test_index,test_index])
      SE = sqrt(abs(cov))
      cov = solve(cov)
      twald = t(par) %*% cov %*% par
      pval = stats::pchisq(twald,df=length(test_index),lower.tail = FALSE)
    }

    return(list(SE=SE,
                twald = twald,
                pval = pval))
  }
}

#' pscl function to get p-value from pscl result
#'
#' @param object pscl result
#'
#' @return pvalue
#' @noRd
pval_getpval<- function(object){
  object$residuals <- stats::residuals(object, type = "pearson")
  kc <- length(object$coefficients$count)
  kz <- length(object$coefficients$zero)
  se <- sqrt(diag(object$vcov))
  coef <- c(object$coefficients$count, object$coefficients$zero)
  if (object$dist == "negbin") {
    coef <- c(coef[1:kc], `Log(theta)` = log(object$theta),
              coef[(kc + 1):(kc + kz)])
    se <- c(se[1:kc], object$SE.logtheta, se[(kc + 1):(kc +
                                                         kz)])
    kc <- kc + 1
  }
  zstat <- coef/se
  pval <- 2 * stats::pnorm(-abs(zstat))
  return(list(SE = se,
              pval = pval))
}

#' Run EM algorithm to fit ZIPG model.
#'
#' @param W Count data
#' @param M Sequencing depth, ZIPG use log(M) as offset by default
#' @param X Formula of covariates of differential abundance
#' @param X_star Formula of covariates of differential varibility
#' @param optim_method Default : 'BFGS'
#' @param init initialization parameters
#' @param return_model whether return full complete imfomation for fitted model
#' @param ref_init Reference
#'
#' @return A list of EM optimize result
#' @noRd
ZIPG_main_EM <- function(W,M,X,X_star,
                         optim_method = 'BFGS',init=NULL,
                         return_model = TRUE,
                         ref_init = NULL){
  test_index=NULL
  doNBZIMM_init = FALSE
  A =1
  W = as.matrix(W)
  X = as.matrix(X)
  X_star = as.matrix(X_star)
  d = dim(X)[2]-1
  d_star = dim(X_star)[2]-1
  length_parms = A*(d+3)+d_star
  if(is.null(init)){
    pscl_init=tryCatch(ZIPG_init_pscl(W,M,X,X_star,A,d,d_star,return_model),
                       error = function(e){
                         warning("pscl ERROR :",conditionMessage(e))
                         return(list(
                           init_model = NULL,
                           init_parms = rep(0,length_parms),
                           init_SE = rep(NA,length_parms),
                           init_pval = rep(NA,length_parms),
                           init_logli = NA
                         ))
                       })

    parms_pscl = pscl_init$init_parms

    pz=colSums(W==0)/length(M)
    gamma_init = log(pz/(1-pz))
    parms_pscl[(A*(d+2)+d_star+1) : (A*(d+3)+d_star)] = gamma_init

    res1 = ZIPG_optim_EM(W,M,X,X_star,
                         A,d,d_star,
                         parms_pscl,
                         optim_method ,NULL)
    if(length(res1)==1 || sum(is.na(pscl_init$init_pval))>=1){
      warning('pscl init optim error')
      warning(parms_pscl)
      if(!is.null(ref_init)){
        res1 = ZIPG_optim_EM(W,M,X,X_star,
                             A,d,d_star,
                             ref_init,
                             optim_method ,NULL)
      }
    }
    if(length(res1)==1){
      res1 = list(
        par = rep(NA,length_parms),
        value = NA,
        hessian = NA
      )
      t_ls1 = list(SE=rep(NA,length_parms),
                   twald = rep(NA,length_parms),
                   pval = rep(NA,length_parms))
      logli1 = rep(NA,A)
    } else {
      logli1 = ZIPG_logli(W,M,X,X_star,A,d,d_star,res1$par,TRUE)
      t_ls1 = pwald(res1,test_index)
      LR1 = logli1-pscl_init$init_logli
      LR1 = 2*c(LR1,sum(LR1))
      LRT_pval1 = stats::pchisq(LR1,df = d_star,lower.tail = FALSE)
    }


    r_df = list(pscl_init = list(init = pscl_init,
                                 res = res1,
                                 wald_test = t_ls1,
                                 logli = logli1))
  }else{
    parms=init
    res1 = ZIPG_optim_EM(W,M,X,X_star,
                         A,d,d_star,
                         parms,
                         optim_method ,NULL)
    if(length(res1)==1){
      warning('pscl_init optim error')
      warning(parms)
      res1 = list(
        par = rep(NA,length_parms),
        value = NA,
        hessian = NA
      )
      t_ls1 = list(SE=rep(NA,length_parms),
                   twald = rep(NA,length_parms),
                   pval = rep(NA,length_parms))
      logli1 = rep(NA,A)
      LRT1 = list(
        LR = rep(NA,A+1),
        pval = rep(NA,A+1)
      )
    } else {
      t_ls1 = pwald(res1,test_index)
      logli1 = ZIPG_logli(W,M,X,X_star,A,d,d_star,res1$par,TRUE)
    }
    r_df = list(pscl_init = list(init = parms,
                                 res = res1,
                                 wald_test = t_ls1,
                                 logli = logli1))

  }
  return(r_df$pscl_init)

}

