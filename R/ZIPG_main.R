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
#' @export ZIPG_main
#'
#' @examples
#' data(Dietary)
#' dat = Dietary
#' ZIPG_res <- ZIPG_main(data = dat$COV,
#' X = ~ALC01+nutrPC1+nutrPC2, X_star = ~ ALC01,
#' W = dat$OTU[,100], M = dat$M )
#' ZIPG_summary(ZIPG_res)
ZIPG_main <- function(data,W,M,X,X_star,
                      return_model = TRUE, pbWald_list = NULL,bWald_list = NULL){
  EM = TRUE
  optim_method='BFGS'

  mf <- match.call(expand.dots = FALSE)
  fX = X
  fX_star = X_star

  if (missing(data)) {
    stop('Missing covariates data frame!')
  }
  tX <- stats::terms(fX, data = data)
  mf = stats::model.frame(data = data,formula = fX,drop.unused.levels = TRUE)
  X <- stats::model.matrix(tX, mf)

  tX <- stats::terms(fX_star, data = data)
  mf = stats::model.frame(data = data,formula = fX_star,drop.unused.levels = TRUE)
  X_star <- stats::model.matrix(tX, mf)

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
    message('Running non-parametric bootstrap wald test \n')
    B =bWald_list$B
    useref = ifelse(is.null(bWald_list$useref),
                    TRUE,bWald_list$useref)
    par_b = sapply(1:B,function(b){
      resample = sample(1:n_sample,n_sample,replace = TRUE)
      Wb = W[resample]
      Mb = M[resample]
      if(d!=0){Xb = X[resample,]}
      if(d_star!=0){X_starb = X_star[resample,]}
      if(useref==TRUE){
        init = res$res$par
      } else {init = NULL}
      res_b = ZIPG_main_EM(W = Wb,M = Mb,X = Xb,X_star = X_starb,
                           optim_method='BFGS',
                           init = init,
                           return_model=FALSE,ref_init = res$res$par)
      return(c(res_b$res$par,res_b$wald_test$SE))
    })
    message('Finish \n')
    n_parms = length(res$res$par)
    par_SE =  t(par_b)[,-c(1:n_parms)]
    par_b = t(par_b)[,1:n_parms]
    varb = stats::var(par_b)
    n_parms = length(res$res$par)
    SE = sqrt(diag(varb))
    twald = res$res$par^2/diag(varb)
    pval = stats::pchisq(twald,df=1,lower.tail = FALSE)

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
    message('Running parametric bootstrap wald test \n')

    fX0 = pbWald_list$X0
    fX_star0 = pbWald_list$X_star0
    all.vars(fX0) %in% all.vars(fX)

    tX <- stats::terms(fX0, data = data)
    mf = stats::model.frame(data = data,formula = fX0,drop.unused.levels = TRUE)
    X0 <- stats::model.matrix(tX, mf)

    tX <- stats::terms(fX_star0, data = data)
    mf = stats::model.frame(data = data,formula = fX_star0,drop.unused.levels = TRUE)
    X_star0 <- stats::model.matrix(tX, mf)

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
                    TRUE,bWald_list$useref)


    res0 = ZIPG_main_EM(W,M,X0,X_star0,
                        optim_method='BFGS',
                        init = NULL,return_model = FALSE)

    LRTdf = d+d_star-(d0+d_star0)
    res$res0 = res0
    LR = 2*(res$logli- res0$logli)
    pval_LRT = stats::pchisq(LR,df =  LRTdf,lower.tail = FALSE)
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
                            zi = TRUE,returnU = FALSE)$W_list
      par_b = c()
      TLRT_b = c()
      TWald_b = c()
      if(useref==TRUE){
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
                             return_model=FALSE,
                             ref_init = res$res$par)

        if(length(H0_index)==1){
          TWald_b[b] = res_b$wald_test$twald[H0_index]
          par_b[b] = res_b$res$par[H0_index]
        }else{
          TWald_b[b] = pwald(res_b$res,H0_index)$twald[1]
        }

      }
      message('Finish \n')

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
#' @return pvalue
#' @export ZIPG_summary
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
  star1=as.character(stats::symnum(ZIPG_res$wald_test$pval,
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
    star1=as.character(stats::symnum(ZIPG_res$bWald$pval,
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
#' @param alpha We construct (1- alpha)\% confidence interval by alpha/2 and (1-alpha/2).
#'
#' @return Table of confidence interval
#' @export ZIPG_CI
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
    lb = stats::qnorm(alpha/2)
    ub = stats::qnorm(1-alpha/2)

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
      lb = stats::qnorm(alpha/2)
      ub = stats::qnorm(1-alpha/2)

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
#' @export ZIPG_simulate
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
    U_simulate_i =matrix(stats::rgamma(n*A,shape = 1/theta,scale = theta),
                         n,A)
    W_simulate_i = matrix(stats::rpois(n*A,lambda*U_simulate_i),n,A)

    if(zi==TRUE){
      zero_sim = matrix(stats::rbinom(n*A,size = 1 ,prob = p_0),n,A,byrow = TRUE)
      W_simulate_i = W_simulate_i*(zero_sim==0)
    }
    W_list[[i]] = W_simulate_i
    U_list[[i]] = U_simulate_i
  }
  if(returnU==TRUE){
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





