#' ZIPG EM algorithm to optimize log-likelihood
#'
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
    if(init == T){
      p_0_col = matrix(rep(p_0,n),n,A,byrow = T)
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

        nb_prob[,k] = dnbinom(W[,k],
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

        nb_prob[,k] = dnbinom(W[,k],
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
        fit <- tryCatch(optim(fn = likelihood_ZIPG_EM, gr = grad_ZIPG_EM,
                              par = par,
                              method = "L-BFGS-B", hessian = T,
                              control = list(maxit = 10000,trace = 0,
                                             fnscale = -1,factr = 1e-6),
                              fix_index = fix_index,fix_par = fix_par,zp=zp),
                        error = function(e){
                          cat("optim ERROR :",conditionMessage(e),"\n")
                          return(conditionMessage(e))
                        })
      } else {
        fit <- tryCatch(optim(fn = likelihood_ZIPG_EM, gr = grad_ZIPG_EM,
                              par = par,
                              method = "BFGS", hessian = T,
                              control = list(maxit = 10000,trace = 0,
                                             fnscale = -1,reltol = 1e-6),
                              fix_index = fix_index,fix_par = fix_par,zp=zp),
                        error = function(e){
                          cat("optim ERROR :",conditionMessage(e),"\n")
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
        fit <- tryCatch(optim(fn = likelihood_ZIPG_EM, gr = grad_ZIPG_EM,
                              par = par,
                              method = "BFGS", hessian = T,
                              control = list(maxit = 10000,trace = 0,
                                             fnscale = -1,reltol = 1e-6),
                              fix_index = fix_index,fix_par = fix_par,zp=zp),
                        error = function(e){
                          cat("optim ERROR :",conditionMessage(e),"\n")
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
          fit <- tryCatch(optim(fn = likelihood_ZIPG_EM, gr = grad_ZIPG_EM,
                                par = par,
                                method = "BFGS", hessian = T,
                                control = list(maxit = 10000,trace = 0,
                                               fnscale = -1,reltol = 1e-6),
                                fix_index = fix_index,fix_par = fix_par,zp=zp),
                          error = function(e){
                            cat("optim ERROR :",conditionMessage(e),"\n")
                            return(conditionMessage(e))
                          })

          if(length(fit)>=2){
            fit$gardient = grad_ZIPG(fit$par,fix_index = fix_index,fix_par = fix_par)
            parms_fit_single[beta_star_index] = fit$par
            fit$hessian = diag(c(fit$hessian),length_parms)
          }
        }
      }
      fit_final <- tryCatch(optim(fn = likelihood_ZIPG_EM, gr = grad_ZIPG_EM,
                            par = parms_fit_single,
                            method = "BFGS", hessian = T,
                            control = list(maxit = 1,trace = 0,
                                           fnscale = -1,reltol = 1e-6),
                            fix_index = NULL,fix_par = NULL,zp=zp),
                      error = function(e){
                        cat("optim ERROR :",conditionMessage(e),"\n")
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

      nb_prob[,k] = dnbinom(W[,k],
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
  zp0 = Estep_zp(parms0,F)
  Q0 = likelihood_ZIPG_EM_zp(parms0,zp0)
  p_zp0 = colMeans(zp0,na.rm = T)
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
    p_zp1 = colMeans(zp0,na.rm = T)
    gamma_zp1 = log(p_zp0/(1-p_zp0))
  }
  fit1$EM_nIter = nIter
  fit1$EM_converge = !(nIter==maxit)

  fit1$init = parms
  return(fit1)
}

#' Log-likelihood for ZIPG
#'
#' @param W
#' @param M
#' @param X
#' @param X_star
#' @param A
#' @param d
#' @param d_star
#' @param parms
#' @param sep
#'
#' @return
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

    nb_prob[,k] = dnbinom(W[,k],
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
#' @param W
#' @param M
#' @param X
#' @param X_star
#' @param A
#' @param d
#' @param d_star
#' @param parms
#' @param zp
#'
#' @return
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

    nb_prob[,k] = dnbinom(W[,k],
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
#' @param W
#' @param M
#' @param X
#' @param X_star
#' @param A
#' @param d
#' @param d_star
#' @param return_model
#'
#' @return
ZIPG_init_pscl <- function(W,M,X,X_star,
                         A=1,d,d_star,
                         return_model = T){
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
#' @param res
#' @param test_index
#'
#' @return
pwald <- function(res,test_index)
{
  solve_hessian<- function(hessian){
    svd_res = svd(hessian)
    smalls = which(svd_res$d<1e-5)
    if(length(smalls)>=1) cat(svd_res$d[smalls])
    svd_res$d[smalls] = 1e-5
    return(svd_res$v %*% diag(1/svd_res$d) %*% t(svd_res$u))
  }
  if(length(res)==1){
    print('optim error')
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
    print('solve hessian error')
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
        pval[i] = pchisq(twald[i],df=1,lower.tail = F)
      }
    }else{
      par = as.matrix(res$par[test_index])
      cov = as.matrix(cov_mat[test_index,test_index])
      SE = sqrt(abs(cov))
      cov = solve(cov)
      twald = t(par) %*% cov %*% par
      pval = pchisq(twald,df=length(test_index),lower.tail = F)
    }

    return(list(SE=SE,
                twald = twald,
                pval = pval))
  }
}


#' Simulate W from ZIPG model
#'
#' @param M Sequencing depth
#' @param X Covariates matrix with intercept, n * (d+1)
#' @param X_star Covariates matrix with intercept, n * (d_star+1)
#' @param A no use, reserved for multi-taxa
#' @param d number of covariates in X
#' @param d_star number of covariates in X_star
#' @param parms model paraneters, input c(beta,beta*,gamma)
#' @param N repetition times
#' @param zi whether generate zero-inflated distribution
#' @param returnU whether return fluctuation factor U
#'
#' @return A list of W generated from ZIPG model with input parameter
#' @export
#'
#' @examples
#' data(Dietary)
#' dat = Dietary
#' sim_M = sample(dat$M,100,replace = TRUE)
#' sim_pre = rep(sample(rep(c(0,1),each = 10)),each = 5)
#' sim_PC1_mean = rep(rnorm(20,mean = 0,sd = 1),each = 5)
#' sim_PC1_error = rnorm(100,0,0.1)
#' sim_PC1 = sim_PC1_mean + sim_PC1_error
#' X = as.matrix(cbind(1,data.frame(X1 = sim_pre,X2 = sim_PC1)))
#' parms = c(-4.23,1,0.45,0.6,1,0,0) #p = 0.5
#' W_sim <- ZIPG_simulate(M = sim_M,X=X,X_star=X,d=2,d_star=2,parms = parms,N=100)
#' hist(W_sim$W_list[[1]])
#' ZIPG_res <- ZIPG_main(data = data.frame(X1 = sim_pre,X2 = sim_PC1),
#' X = ~X1+X2, X_star = ~ X1,W = W_sim$W_list[[2]], M = sim_M )
#' ZIPG_summary(ZIPG_res)
ZIPG_simulate <- function(M,X,X_star,
                        A=1,d,d_star,
                        parms,N,
                        zi = TRUE,returnU = FALSE){
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

  W_list = list()
  U_list = list()

  for(k in 1:A)
  {
    lambda[,k] = exp(X %*% beta[,k])*M
    theta[,k] = exp(X_star %*% beta_star[,k])
  }

  for(i in 1:N){
    U_simulate_i =matrix(rgamma(n*A,shape = 1/theta,scale = theta),
                         n,A)
    W_simulate_i = matrix(rpois(n*A,lambda*U_simulate_i),n,A)

    if(zi==T){
      zero_sim = matrix(rbinom(n*A,size = 1 ,prob = p_0),n,A,byrow = T)
      W_simulate_i = W_simulate_i*(zero_sim==0)
    }
    W_list[[i]] = W_simulate_i
    U_list[[i]] = U_simulate_i
  }
  if(returnU==T){
    return(list(
      W_list = W_list,
      U_list = U_list,
      lambda = lambda,
      theta = theta
    ))
  } else {
    return(list(
      W_list = W_list,
      U_list = NULL,
      lambda = lambda,
      theta = theta
    ))
  }
}






#' Run EM algorithm to fit ZIPG model.
#'
#' @param W
#' @param M
#' @param X
#' @param X_star
#' @param optim_method
#' @param init
#' @param return_model
#' @param ref_init
#'
#' @return
ZIPG_main_EM <- function(W,M,X,X_star,
                       optim_method = 'BFGS',init=NULL,
                       return_model = T,
                       ref_init = NULL){
  test_index=NULL
  doNBZIMM_init = F
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
                         cat("pscl ERROR :",conditionMessage(e))
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
      print('pscl init optim error')
      print(parms_pscl)
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
      logli1 = ZIPG_logli(W,M,X,X_star,A,d,d_star,res1$par,T)
      t_ls1 = pwald(res1,test_index)
      LR1 = logli1-pscl_init$init_logli
      LR1 = 2*c(LR1,sum(LR1))
      LRT_pval1 = pchisq(LR1,df = d_star,lower.tail = F)
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
      print('pscl_init optim error')
      print(parms)
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
      logli1 = ZIPG_logli(W,M,X,X_star,A,d,d_star,res1$par,T)
    }
    r_df = list(pscl_init = list(init = parms,
                                 res = res1,
                                 wald_test = t_ls1,
                                 logli = logli1))

  }
  return(r_df$pscl_init)

}


#' pscl function to get p-value from pscl result
#'
#' @param object
#'
#' @return
pval_getpval<- function(object){
  object$residuals <- residuals(object, type = "pearson")
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
  pval <- 2 * pnorm(-abs(zstat))
  return(list(SE = se,
              pval = pval))
}


#' Fit zero-inflated poisson-gamma model via EM Algorithm
#'
#' @param data Data.frame for covariates of interest
#' @param W Count data
#' @param M Sequencing depth, ZIPG use log(M) as offset by default
#' @param X Formula of covariates of differential abundance
#' @param X_star Formula of covariates of differential varibility
#' @param return_model whether return full complete imfomation for fitted model
#' @param pbWald_list A list of arguments for parameteric bootstrap Wald test,
#' B for bootstrap sample size, X0 and X_star0 for formula of covariates included in H0
#' @param bWald_list A list of arguments for non-parameteric bootstrap Wald test, B for bootstrap sample size,
#'
#' @return A list of ZIPG fitted model. Use ZIPG_summary() for a quick look at the results.
#' @export
#'
#' @examples
#' data(Dietary)
#' dat = Dietary
#' ZIPG_res <- ZIPG_main(data = dat$COV,
#' X = ~ALC01+nutrPC1+nutrPC2, X_star = ~ ALC01,
#' W = dat$OTU[,100], M = dat$M )
#' ZIPG_summary(ZIPG_res)
ZIPG_main <- function(data,W,M,X,X_star,
                      return_model = T, pbWald_list = NULL,bWald_list = NULL){
  EM = TRUE
  optim_method='BFGS'

  mf <- match.call(expand.dots = FALSE)
  fX = X
  fX_star = X_star

  if (missing(data)) {
    stop('Missing covariates data frame!')
  }
  tX <- terms(fX, data = data)
  mf = model.frame(data = data,formula = fX,drop.unused.levels = TRUE)
  X <- model.matrix(tX, mf)

  tX <- terms(fX_star, data = data)
  mf = model.frame(data = data,formula = fX_star,drop.unused.levels = TRUE)
  X_star <- model.matrix(tX, mf)

  d = dim(X)[2]-1
  d_star = dim(X_star)[2]-1
  n_sample = length(M)
  test_index = NULL
  A = 1
  res= ZIPG_main_EM(W,M,X,X_star,
                    optim_method='BFGS',init=NULL,
                    return_model = return_model,ref_init = NULL)
  res$info = list(
    d =d,d_star = d_star
  )

  if(!is.null(bWald_list)){
    cat('Running non-parametric bootstrap wald test \n')
    B =bWald_list$B
    useref = ifelse(is.null(bWald_list$useref),
                    T,bWald_list$useref)
    par_b = sapply(1:B,function(b){
      resample = sample(1:n_sample,n_sample,replace = T)
      Wb = W[resample]
      Mb = M[resample]
      if(d!=0){Xb = X[resample,]}
      if(d_star!=0){X_starb = X_star[resample,]}
      if(useref==T){
        init = res$res$par
      } else {init = NULL}
      res_b = ZIPG_main_EM(W = Wb,M = Mb,X = Xb,X_star = X_starb,
                           optim_method='BFGS',
                           init = init,
                           return_model=F,ref_init = res$res$par)
      return(c(res_b$res$par,res_b$wald_test$SE))
    })
    cat('Finish \n')
    n_parms = length(res$res$par)
    par_SE =  t(par_b)[,-c(1:n_parms)]
    par_b = t(par_b)[,1:n_parms]
    varb = var(par_b)
    n_parms = length(res$res$par)
    SE = sqrt(diag(varb))
    twald = res$res$par^2/diag(varb)
    pval = pchisq(twald,df=1,lower.tail = F)

    res$bWald = list(
      B = B,
      par_b = par_b,
      par_SE = par_SE,
      vcov = varb,
      SE= SE ,
      twald = twald,
      pval = pval
    )
  }

  if(!is.null(pbWald_list)){
    cat('Running parametric bootstrap wald test \n')

    fX0 = pbWald_list$X0
    fX_star0 = pbWald_list$X_star0
    all.vars(fX0) %in% all.vars(fX)

    tX <- terms(fX0, data = data)
    mf = model.frame(data = data,formula = fX0,drop.unused.levels = TRUE)
    X0 <- model.matrix(tX, mf)

    tX <- terms(fX_star0, data = data)
    mf = model.frame(data = data,formula = fX_star0,drop.unused.levels = TRUE)
    X_star0 <- model.matrix(tX, mf)

    d0 = dim(X0)[2]-1
    d_star0 = dim(X_star0)[2]-1

    if(FALSE %in% (colnames(X0) %in% colnames(X)))
    {
      stop("H0 not nested in H1 !")
    }
    if(FALSE %in% (colnames(X_star0) %in% colnames(X_star)))
    {
      stop("H0 not nested in H1 !")
    }

    fX_index = (colnames(X) %in% colnames(X0))
    fX_star_index = (colnames(X_star) %in% colnames(X_star0))
    H0_index = which(c(fX_index,fX_star_index)==FALSE)

    useref = ifelse(is.null(bWald_list$useref),
                    T,bWald_list$useref)


    res0 = ZIPG_main_EM(W,M,X0,X_star0,
                        optim_method='BFGS',
                        init = NULL,return_model = F)

    LRTdf = d+d_star-(d0+d_star0)
    res$res0 = res0
    LR = 2*(res$logli- res0$logli)
    pval_LRT = pchisq(LR,df =  LRTdf,lower.tail = F)
    res$LRT  = list(LR = LR,pval =pval_LRT)

    if(!is.null(pbWald_list$B)){
      B = pbWald_list$B
      if(length(H0_index)==1){
        TWald0 = res$wald_test$twald[H0_index]
      }else{
        TWald0 = pwald(res$res,H0_index)$twald[1]
      }
      W_sim = ZIPG_simulate(M,X0,X_star0,
                            A=1,d0,d_star0,
                            res0$res$par,N=B,
                            zi = T,returnU = F)$W_list
      par_b = c()
      TLRT_b = c()
      TWald_b = c()
      if(useref==T){
        init = res$res$par
        init0 = res$res$par[-H0_index]
      } else {
        init = NULL
        init0 = NULL
      }
      TLRT_b = NA
      for(b in 1:B){
        Wb = W_sim[[b]]
        res_b = ZIPG_main_EM(Wb,M,X,X_star,
                             optim_method='BFGS',
                             init =init,
                             return_model=F,
                             ref_init = res$res$par)

        if(length(H0_index)==1){
          TWald_b[b] = res_b$wald_test$twald[H0_index]
          par_b[b] = res_b$res$par[H0_index]
        }else{
          TWald_b[b] = pwald(res_b$res,H0_index)$twald[1]
        }

      }
      cat('Finish \n')

      pbWald = (1+sum(TWald_b>=TWald0))/(1+B)
      res$pbWald = list(
        par_b= par_b,
        pval_pbWald = pbWald,
        H0_index = H0_index
      )
    }
  }
  return(res)

}


#' Summary for ZIPG_main() result.
#'
#' @param ZIPG_res Result from ZIPG_main()
#' @param type Type of hypothesis testing method, 'Wald','bWald' or 'pbWald'.
#'
#' @return
#' @export
#'
#' @examples
#' data(Dietary)
#' dat = Dietary
#' ZIPG_res <- ZIPG_main(data = dat$COV,
#' X = ~ALC01+nutrPC1+nutrPC2, X_star = ~ ALC01,
#' W = dat$OTU[,100], M = dat$M )
#' ZIPG_summary(ZIPG_res)
ZIPG_summary<-function(ZIPG_res,type='Wald'){
  d = ZIPG_res$info$d
  d_star = ZIPG_res$info$d_star
  n_parms = length(ZIPG_res$res$par)
  if(type =='Wald'){
  star1=as.character(symnum(ZIPG_res$wald_test$pval,
                            cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                            symbols = c("***", "**", "*", ".", " ")))
  Wald_mat = data.frame(cbind(ZIPG_res$res$par,
                   ZIPG_res$wald_test$SE,
                   ZIPG_res$wald_test$pval),star1)
  cat('          ZIPG Wald  \n')
  colnames(Wald_mat)<- c('Estimation','SE','pval','')
  rownames(Wald_mat) = c(paste('beta', c(0:d),sep = ''),
                         paste('beta', c(0:d_star),'*',sep = ''),
                         'gamma')

  print(Wald_mat,digits = 3)
  return(Wald_mat)
  }else if(type =='bWald'){
    star1=as.character(symnum(ZIPG_res$bWald$pval,
                              cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                              symbols = c("***", "**", "*", ".", " ")))
    Wald_mat = data.frame(cbind(ZIPG_res$res$par,
                                ZIPG_res$bWald$SE,
                                ZIPG_res$bWald$pval),star1)
    cat('          ZIPG bWald   \n')
    colnames(Wald_mat)<- c('Estimation','SE','pval','')
    rownames(Wald_mat) = c(paste('beta', c(0:d),sep = ''),
                           paste('beta', c(0:d_star),'*',sep = ''),
                           'gamma')

    print(Wald_mat,digits = 3)
    return(Wald_mat)
  }else if(type =='pbWald'){
    H0_index = ZIPG_res$pbWald$H0_index
    parnames = c(paste('beta', c(0:d),sep = ''),
                 paste('beta', c(0:d_star),'*',sep = ''),
                 'gamma')
    cat('   ZIPG pbWald \n H0:',parnames[H0_index],'= 0 \n')
    cat(' pvalue = ', signif(ZIPG_res$pbWald$pval_pbWald,3),'\n')
    return(ZIPG_res$pbWald$pval_pbWald)
  } else {
    stop("Type must be one of 'Wald', 'bWald' or 'pbWald' ! ")
  }
}

#' Get confidence interval from ZIPG model
#'
#' @param ZIPG_res Result from ZIPG_main()
#' @param type Type of hypothesis testing method, 'Wald' or 'bWald'.
#' @param CI_type Type of confidence interval, 'Wald','bWald' or 'pbWald'.
#' @param alpha We construct (1- alpha)% confidence interval by alpha/2 and (1-alpha/2).
#'
#' @return
#' @export
#'
#' @examples
#' data(Dietary)
#' dat = Dietary
#' ZIPG_res <- ZIPG_main(data = dat$COV,
#' X = ~ALC01+nutrPC1+nutrPC2, X_star = ~ ALC01,
#' W = dat$OTU[,100], M = dat$M )
#' ZIPG_CI(ZIPG_res)
ZIPG_CI<-function(ZIPG_res,type='Wald',CI_type = 'normal',alpha = 0.05){
  d = ZIPG_res$info$d
  d_star = ZIPG_res$info$d_star
  n_parms = length(ZIPG_res$res$par)
  if(type =='Wald'){
    lb = qnorm(alpha/2)
    ub = qnorm(1-alpha/2)

    CI_mat = data.frame(
      ZIPG_res$res$par,
      lb = ZIPG_res$res$par + lb*ZIPG_res$wald_test$SE,
      ub = ZIPG_res$res$par + ub*ZIPG_res$wald_test$SE
    )

    cat('        ZIPG Wald Confidence interval \n')
    colnames(CI_mat)<- c('Estimation','lb','ub')
    rownames(CI_mat) = c(paste('beta', c(0:d),sep = ''),
                           paste('beta', c(0:d_star),'*',sep = ''),
                           'gamma')
    print(CI_mat,digits = 3)
    return(CI_mat)
  }else if(type =='bWald'){
    if(CI_type =='normal'){
      lb = qnorm(alpha/2)
      ub = qnorm(1-alpha/2)

      CI_mat = data.frame(
        ZIPG_res$res$par,
        lb = ZIPG_res$res$par + lb*ZIPG_res$bWald$SE,
        ub = ZIPG_res$res$par + ub*ZIPG_res$bWald$SE
      )

      cat('        ZIPG Wald Confidence interval \n')
      colnames(CI_mat)<- c('Estimation','lb','ub')
      rownames(CI_mat) = c(paste('beta', c(0:d),sep = ''),
                           paste('beta', c(0:d_star),'*',sep = ''),
                           'gamma')
      print(CI_mat,digits = 3)
      return(CI_mat)
    }else if (CI_type =='percentile'){
      warning('Normality-based CI or BCa method are more recommended !')
        ub = c()
        lb = c()
      for(index in 1:length(ZIPG_res$res$par)){
        sortb = sort(ZIPG_res$bWald$par_b[,index])
        B = length(sortb)
        ub[index] = sortb[ B*(1-alpha/2)]
        lb[index] = sortb[B*alpha/2+1]
      }
        CI_mat = data.frame(
          ZIPG_res$res$par,
          lb = lb,
          ub = ub
        )

        cat('        ZIPG Wald percentlile Confidence interval \n')
        colnames(CI_mat)<- c('Estimation','lb','ub')
        rownames(CI_mat) = c(paste('beta', c(0:d),sep = ''),
                             paste('beta', c(0:d_star),'*',sep = ''),
                             'gamma')
        print(CI_mat,digits = 3)
        return(CI_mat)
    }else{
      stop("Type must be one of 'normal' or 'percentile'! ")
    }
  } else {
    stop("Type must be one of 'Wald' or 'bWald'! ")
  }
}

