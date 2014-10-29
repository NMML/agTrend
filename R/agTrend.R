

###
### Functions
###

#' @title Create list of prior parameters for \code{mcmc.aggregate(...)} trend estimation.
#' 
#' @description This function takes user input to create the list of prior paremters necessary for 
#' Bayesian inference of aggregated trends. The prior parameter list is a necessary argument for 
#' the \code{\link{mcmc.aggregate}} function.
#' 
#' @param trend.limit Bounds on the individual site trends.
#' @param model.data A data set which provides the individual site models used in \code{\link{mcmc.aggregate}}.
#' @param gamma.mean Prior mean of the (multivariate) normal prior distribution for gamma
#' @param gamma.prec Precision (inverse of variance) of the gamma prior distribution
#' 
#' @return A named list with entries for each of the parameters for which prior specification is necessary.
#' 
#' @details Using the site model data set \code{model.data} and a soft bound on the individual site
#' trends (\code{trend.limit}) this function creates a sensible default for the list of prior 
#' parameters necessary for Bayesian inference of aggregated trends. The value of 
#' \code{trend.limit} should be in terms of % growth (e.g., \code{trend.limit=0.2} implies an upper limit 
#' of 20 percent growth at each individual site). The lower bound is calculated as 
#' 1/(1+trend.limit). This is the ower bound that is symmetric about zero on the log-scale.
#' For \code{trend.limit=0.2}, the lower bound is approximately -0.17. 
#' 
#' Specifically, the prior distribution for the linear trend parameters (beta) 
#' in log abundance at each site is set to a bivariate normal with mean \code{m.i=c(0,0)}
#' and precision matrix \code{Q.i=diag(c(0,q))}, where \code{q} is chosen such that
#' Pr[1/(1+trend.limit) < exp(beta[2])-1 < trend.limit] = 0.95.
#' 
#' 
#' @export
#' @import Matrix

defaultPriorList = function(trend.limit=NULL, model.data, gamma.mean, gamma.prec){
  if(is.null(trend.limit) & any(c("RW","lin")%in%model.data$trend)) stop("A value must be provided for 'trend.limit' when using 'lin' and 'RW' models!")
  if(!is.null(trend.limit)) Q.trend = 1/(log(1+trend.limit)/2)^2
  use.zi = !is.null(model.data$avail)
  if(!use.zi) use.alpha=FALSE
  if(use.zi) use.alpha = any(model.data$avail=="RW")
  if(missing(gamma.mean) & missing(gamma.prec)){
    gamma=NULL
  } else {
    if(missing(gamma.mean)) {
      gamma.0 = rep(0,sqrt(length(gamma.prec)))
      gamma=list(gamma.0=gamma.0, Q.gamma=gamma.prec)
    } else if(missing(gamma.prec)) {
      Q.gamma=matrix(0,nrow=length(gamma.0))
      gamma=list(gamma.0=gamma.mean, Q.gamma=Q.gamma)
    } else {gamma=list(gamma.0=gamma.mean, Q.gamma=gamma.prec)}
  }
  beta=list(
    beta.0=rep(0,sum(model.data$trend=="const")+(sum(model.data$trend!="const")*2)),
    Q.beta=Matrix(diag(
      unlist(lapply(model.data$trend, 
                  function(x){
                    if(x!="const") {
                      return(c(0,Q.trend))
                    }else {
                      return(c(0))
                    }
                  }))))
  ) 

  if(any(model.data$trend=="RW")){
    tau=list(a.tau=0.5, b.tau=5.0E-5)
  } else tau=NULL
  zeta=list(a.zeta=0.5, b.zeta=5.0E-5)
  if(use.zi & use.alpha) {
    phi=list(a.phi=0.5, b.phi=5.0E-5)
  } else phi=NULL
  if(use.zi){
    tmp=model.data$avail[model.data$avail!="none"]
    theta.0=rep(0,sum(tmp=="const")+(sum(tmp!="const")*2))
    theta=list(
      theta.0=theta.0,
      Q.theta=Matrix(diag(
        unlist(lapply(tmp, 
                      function(x){
                        if(x!="const") {
                          return(c(0,0))
                        }else {
                          return(c(0))
                        }
                      }))))
      
    )
  } else theta=NULL
  return(list(gamma=gamma, beta=beta, theta=theta, zeta=zeta, tau=tau, phi=phi)) 
}

#' @title Unscaled precision matrix for a random walk of order \code{p}.
#' 
#' @description This function calculates the entries for a precision matrix for a random walk of order \code{p}.
#' It is used in the MCMC augementation with \code{p=2}, but, it is presented here for general use.
#' 
#' @param n The length of the RW(p) process. The length must be greater than or equal to \code{2*p}.
#' @param p The order of the random walk process. 
#' 
#' @return A matrix
#' 
#' @references H. Rue and L. Held (2005) Gaussian Markov Random Fields. Chapman & Hall/CRC. 263 pp.
#' @export
#' @import mgcv
iar.Q <- function(n,p)
{
  if(n<2*p) stop("n must be >= 2*p\n")
  tmp1 <- out <- matrix(0,n,n)
  tmp2 <- ((-1)^c(0:p))*choose(p,c(0:p))
  for(i in (p+1):n)
  {
    tmp1[,i] <- c(rep(0,i-p-1), tmp2, rep(0, n-i))
  }
  for(i in n:(p+1))
  {
    tmp4 <- tmp1[,c(1:n)[tmp1[i,]!=0]]
    tmp5 <- tmp1[i,tmp1[i,]!=0]
    tmp6 <- t(t(tmp4) * tmp5)
    tmp6[i,] <- 0
    out[i,] <- -rowSums(tmp6)
  }
  out[1:p,] <- out[n:(n-p+1),n:1]
  return(diag(apply(out,1,sum)) - out)
}

getSiteREInits <- function(data, abund.name, site.name, time.name, ln.adj, 
                           Q0.omega.s, r0.omega.s, a.tau, b.tau, newdata, model.data){
  #require(mgcv)
  data$y <- log(data[,abund.name]+ln.adj)
  omega <- NULL
  tau <- NULL
  zeta <- NULL
  b <- NULL
  num <- 1
  for(i in levels(data[,site.name])){
    models <- model.data[model.data[,site.name]==i,]
    if(models$trend%in%c("const","lin")){
      if(models$trend=="const") fit <- lm(data[data[,site.name]==i,"y"]~1, )
      else fit <- lm(data[data[,site.name]==i,"y"]~data[data[,site.name]==i,time.name])
      omega <- c(omega, rep(0,nrow(newdata)))
      b <- c(b, as.vector(coef(fit)))
      zeta <- c(zeta, 1/var(residuals(fit)))
    }
    else{
      form <- as.formula(paste("y ~ s(", time.name,")", sep=""))
      fit <- gam(form, data=data, subset=(data[,site.name]==i))
      omega.tmp <- predict(fit, newdata=newdata)
      lm.fit <- lm(omega.tmp~newdata[,1])
      omega.tmp <- matrix(residuals(lm.fit), ncol=1)
      tau <- c(tau, (a.tau+r0.omega.s/2 -1)/(b.tau + as.vector(crossprod(omega.tmp, Q0.omega.s)%*%omega.tmp)/2))
      omega <- c(omega, omega.tmp)
      b <- c(b, as.vector(coef(lm.fit)))
      zeta <- c(zeta, 1/var(residuals(fit)))
    }
    num <- num+1
  }
  zeta <- ifelse(zeta==Inf, max(zeta[zeta<Inf]), zeta)
  return(list(omega=omega, zeta=zeta, tau=tau, b=b))
}

getSiteZIInits <- function(data.orig, abund.name, site.name, time.name,
                           Q0.alpha.s, r0.alpha.s, a.phi, b.phi, newdata, model.data){
  #require(mgcv)
  data.orig$y <- 1.0*(data.orig[,abund.name]>0)
  alpha <- NULL
  phi <- NULL
  theta <- NULL
  num <- 1
  for(i in levels(data.orig[,site.name])){
    models <- model.data[model.data[,site.name]==i,]
    if(models$avail=="none") next
    else if(models$avail%in%c("const","lin")){
      if(models$avail=="const") fit <- glm(data.orig[data.orig[,site.name]==i,"y"]~1, family=binomial(link="probit"))
      else fit <- glm(data.orig[data.orig[,site.name]==i,"y"]~data.orig[data.orig[,site.name]==i,time.name], family=binomial("probit"))
      alpha <- c(alpha, rep(0,nrow(newdata)))
      theta <- c(theta, as.vector(coef(fit)))
    }
    else{
      form <- as.formula(paste("y ~ s(", time.name,")", sep=""))
      fit <- gam(form, data=data.orig[data.orig[,site.name]==i,],family=binomial(link="probit"))
      alpha.tmp <- predict(fit, newdata=newdata)
      lm.fit <- lm(alpha.tmp~newdata[,1])
      alpha.tmp <- matrix(residuals(lm.fit), ncol=1)
      phi <- c(phi, (a.phi+r0.alpha.s/2 -1)/(b.phi + as.vector(crossprod(alpha.tmp, Q0.alpha.s)%*%alpha.tmp)/2))
      alpha <- c(alpha, alpha.tmp)
      theta <- c(theta, as.vector(coef(lm.fit)))
    }
    num <- num+1
  }
  return(list(alpha=alpha, phi=phi, theta=theta))
}

#' @title 
#' Posterior predictive sampling, aggregtion of abundance counts, and linear trend summary
#' 
#' @description 
#' Function for sampling from the posterior predictive distribution of abundance (counts) at 
#' individual sites. Then aggregating the counts over the specified aggregation variable. 
#' 
#' @param start  The starting time for trend estimation
#' @param end  The end time for trend estimation
#' @param data  A \code{data.frame} that contains the abundance survey data.  
#' @param obs.formula  A formula object specifying the model for the observation data
#' @param aggregation  A factor variable. Aggregation is performed over each level of the factor.
#' @param model.data  A data frame giving the augmentation model for each site. See 'Details'
#' @param rw.order A names list, e.g., \code{list(omega=2, alpha=2)}, that gives the order of the RW process for the trend and ZI models.
#' @param abund.name  A character string giving the name of the data to be aggregated
#' @param time.name  A character string giving the name of the time variable
#' @param site.name  A character string giving the name of the site variable. The variable should be a factor
#' @param sig.abund  A numeric vector the same length as \code{nrow(data)} which contains the known observation error standard deviations.
#' @param incl.zeros If \code{incl.zeros=TRUE} (default), and zero inflation models are used, then the zeros become part of the 'true' abundance process and are used for trend estimation and abundance prediction. 
#' @param forecast A logical indicating whether to allow forecasting aggregations past the last observed time.
#' @param ln.adj  The adjustment for taking logs of counts if zeros are present, e.g., log(n + ln.adj).
#' @param upper  A data frame containing the upper bounds for augmentation of each site. See 'Details'
#' @param lower  A data frame containing the lower bounds for augmentation of each site. See 'Details'
#' @param burn  The length of burnin for the MCMC augmentation and aggregation.
#' @param iter  The number of MCMC iterations 
#' @param thin  The amount of thinning of the MCMC sample. e.g., \code{thin=5} implies keeping every 5th MCMC sample for inference
#' @param prior.list  A named list containing the prior distributions for the parameters and random effects
#' @param keep.site.abund  Logical. Should the augmented site abundance be retained.
#' @param keep.site.param  Logical. Should the site augmentation parameters be retianed.
#' @param keep.obs.param  Logical. Should the observation parameters (gamma) be retianed. 
#' 
#' @details 
#' This function is the workhorse for the \code{agTrend} package. It performs 
#' MCMC sampling of the posterior predictive distribution of the abundance at 
#' each site at each time, \eqn{N_{st}}{N_st}. The abundance at each site is 
#' modeled, in its most general form, with a zero-inflated, nonparameteric model,
#' \deqn{z_{st} = \beta_{s0} + \beta_{s1}t + \omega_{st} + \delta_{st} \mbox{ if } N_{st}>0,}{z_st = beta_s0 + beta_s1 * t + omega_st + delta_st if N_st > 0,}
#' where \eqn{\beta_{s0}+\beta_{s1}t}{beta_s0 + beta_s1 * t} is the linear 
#' trend, \eqn{\omega}{omega} is a random walk (of order 1 or 2) (RW), and \eqn{\delta_{st}}{delta_st} is an iid normal error variable. The zero-inflation part is added via 
#' the probit regression model
#' \deqn{\mbox{probit}\{P(N_{st}>0)\} = \theta_{s0} + \theta_{s1}t + \alpha_{st},}{probit{P(N_st > 0)} = theta_s0 + theta_s1 * t + alpha_st,}
#' where \eqn{\theta_{s0}}{theta_s0} and \eqn{\theta_{s1}}{theta_s1} are linear regression coefficients and \eqn{\alpha}{alpha} is a RW model.
#' 
#' In order to account for observation effects or changing methodology through time one can specify an \code{obs.model}. The \code{obs.model} is a R formula object 
#' that specifies variables in \code{data} that can account for differences due to sampling methodology alone. If \code{obs.model} is provided, the observation model is specified as
#' \deqn{y_{st} = x_{st}'\gamma + z_{st} + \epsilon_{st},}{y_st = x_st gamma + z_st + eps_st,}
#' where \eqn{y_{st}}{y_st} is the raw observed log abundance, 
#' \eqn{x_{st}}{x_st} is a vector of covariates used to standardize the observed abundance, \eqn{\gamma}{gamma} is a vector of coefficients, and 
#' \eqn{[\epsilon_{st}]=N(0,\sigma^2_{st})}{[eps_st]=N(0,sigma_st^2)}. Currently, \eqn{\sigma_{st}}{sigma_st} is considered to be known and is specified
#' as a column in \code{data} by the \code{sig.abund} argument. Thus, \eqn{z_{st}}{z_st} represents the standardized (wrt the survey method) abundance. See \code{demo(wdpsNonpups)} for an example.
#' 
#' For each iteration of the MCMC sampler, the complete data is aggregated (i.e., summed) over all sites 
#' within a specified region (defined by the \code{aggregation} argument). Thus, one can obtain a posterior predictive sample from the aggregated abundance for every time between the first
#' time in the data set and the last. By using the posterior predictive distribution, we can account for parameter uncertainty at each site as well 
#' as sampling variability of the survey. Even though we are using the Bayesian inference paradigm, we still capture the essence of frequentist inference by accounting for the 'replication'
#' of the survey effort by using the predictive distribution even for times and places where we have survey data. Using the aggregations, the average linear 
#' trend is calculated for all years from \code{start} to \code{end} for each MCMC iteration.  
#' 
#' The \code{model.data} data.frame can be provided to reduce the most general model given above to submodels when there is not adequate data to fit the full model
#' or, if zero-inflation is not necessary. The \code{model.data} must be a \code{data.frame} with columns 
#' \itemize{
#' \item{\code{site.name}- Column of all the sites in \code{data}. The name must be given by \code{site.name}}
#' \item{\code{trend}- A column of all augmentation models to be used at each site, can be one, any only one, of the following: \code{c("const","lin","RW")}
#' for a constant mean, linear mean, RW model respectively. Each generalization includes all previous models, e.g., an RW model contains a linear slope and 
#' an intercept.}
#' \item{\code{avail}- The same as the \code{trend} column, except this specifies the model for the zero-inflation component and can additionally include \code{"none"}
#' if a zero-inflation model is not desired.}
#' }
#' 
#' Examples of this function's use can be seen in the \code{demo(package="agTrend")} files.
#' 
#'  
#'   
#' @return
#' A named list with the following elements:
#' \item{trend.summary}{Summary of the posterior predictive linear trend}
#' \item{aggregation.summary}{Summary of the site aggregations for every time between the first and last}
#' \item{site.summary}{A summary of the abundance augmentation for every site and every time.}
#' \item{mcmc.sample}{A named list containing all of the MCMC sample after thinning.}
#' \item{original.data}{The original data in \code{data}.}
#' 
#' @export
#' @import Matrix
#' @import coda
#' @import truncnorm
#' 
mcmc.aggregate <- function(start, end, data, obs.formula=NULL, aggregation, model.data, rw.order=NULL,
                           abund.name, time.name, site.name, sig.abund, incl.zeros=TRUE, forecast=FALSE, 
                           ln.adj=0, upper=Inf, lower=-Inf,
                           burn, iter, thin=1, prior.list=NULL, 
                           keep.site.abund=FALSE, keep.site.param=FALSE, keep.obs.param=FALSE
                           ){
  #require(Matrix)
  #require(coda)
  use.trunc <- !is.null(model.data$avail) | (any(upper!=Inf) | any(lower!=-Inf))
  #if(use.trunc) require(truncnorm)
  use.zi <- !is.null(model.data$avail)
  if(use.zi){
    use.zi <- any(model.data$avail!="none")
  }
  if(use.zi) {
    use.alpha <- any(model.data$avail=="RW")
  } else use.alpha <- FALSE
  use.gam <- !is.null(obs.formula)
  use.omega <- any(model.data$trend=="RW")

  if(use.zi & !all(model.data$avail%in%c("none","const","lin","RW"))){
    stop("\n Error: currently the only models for availability are: 'none', 'const', 'lin', or 'RW'\n")
  }
  if(!all(model.data$trend%in%c("const","lin","RW"))){
    stop("\n Error: currently the only models for trend are: 'const', 'lin', or 'RW'\n")
  }

  
  # MODEL ###
  
  #(Note: Normals are parameterized with precision matrices)
  # y = Xy %*% gamma + M %*% z + eps
  # eps ~ N(0,Qy), Qy is known
  # Z = Xz %*% beta + omega + delta
  # delta ~ N(0,Qdelta), Qdelta is diagonal with elements zeta
  # omega ~ N(0, Qomega); Qomega is a block diagonal RW(p) model with precision tau
  
  # DATA MANIPULATION / PREPARATION ###
  if(use.zi) ln.adj <- 0
  d.start <-  min(data[,time.name])
  if(!forecast){
    d.end <- max(data[,time.name])
  } else {d.end <- max(max(data[,time.name]), end)}
  d.yrs <- c(d.start:d.end)
  if(missing(start)) start <- d.start
  if(missing(end)) end <- d.end
  if(start<d.start) start <- d.start
  if(end>d.end & !forecast) end <- d.end
  data <- data[order(data[,site.name],data[ ,time.name]),]
  yrs <- c(start:end)
  siteNms <- levels(data[,site.name])
  n.site <- length(siteNms)
  site.idx <- rep(levels(data[,site.name]), each=length(d.yrs))
  yr.idx <- rep(d.yrs, n.site)
  n.year <- length(d.yrs)
  bigN <- n.site*n.year
  Isite <- Diagonal(length(siteNms))
  if(missing(aggregation)){
    data$Total <- "Total"
    aggregation <- "Total"
  }
# Construct indexes for the entire augemented data set
  #aggregation index
  data[,aggregation] <- factor(data[,aggregation])
  ag.data <- unique(data[,c(site.name,aggregation)])
  ag.data <- ag.data[order(ag.data[,aggregation]),]
  agg.idx <- merge(data.frame(site.idx), ag.data, by.y=site.name, by.x=1)[,aggregation]
  # upper and lower bounds
  if(is.data.frame(upper)) upper <- log(merge(data.frame(site.idx), upper, by.y=site.name, by.x=1)$upper + ln.adj)
  if(is.data.frame(lower)) lower <- log(merge(data.frame(site.idx), lower, by.y=site.name, by.x=1)$lower + ln.adj)
  # omega
  if(use.omega){
    omegarw.data <- merge(data.frame(site.idx), model.data, by.y=site.name, by.x=1)
    omegarw <- c(omegarw.data$trend=="RW")
  }
  # alpha and fixed q
  q.data <- merge(data.frame(site.idx), model.data, by.y=site.name, by.x=1)
  if(use.alpha){
    qrw <- c(q.data$avail=="RW")[q.data[,1]%in%model.data[model.data$avail!="none",site.name]]
    Nqrw <- sum(qrw)
  }
  qzi <- c(q.data$avail!="none")
  Nzi <- sum(qzi)
  # upper and lower bounds for zi draws
  data.orig <- data
  if(use.zi){
    obs.idx <- (as.numeric(data[,site.name])-1)*length(d.yrs) + data[,time.name]-min(data[,time.name])+1
    a.zi <- rep(-Inf, bigN)
    a.zi[obs.idx[data[,abund.name]>0]] <- 0
    b.zi <- rep(Inf, bigN)
    b.zi[obs.idx[data[,abund.name]==0]] <- 0
    data <- data[data[,abund.name]>0,]
  }
  obs.idx <- (as.numeric(data[,site.name])-1)*length(d.yrs) + data[,time.name]-min(data[,time.name])+1
  # create observation matrix 
  M <- Diagonal(bigN)[obs.idx,]
  y <- log(data[,abund.name] + ln.adj)
  # create design matrix for gamma
  if(use.gam){
    Xy <- Matrix(model.matrix(obs.formula, data=data))
  } else{Xy <- NULL}
  # create design matrix for beta
  Xz <- .bdiag(lapply(model.data$trend, 
                   function(x){
                      if(x=="const"){return(matrix(rep(1,n.year), ncol=1))}
                      else return(cbind(rep(1,n.year),d.yrs))                         
                     }))
  # create constraint matrix for omega
  if(!is.null(rw.order$omega)){
    omega.order=rw.order$omega
  } else{omega.order=1}
  if(use.omega){
    if(omega.order==2){
    Xz.rw <- .bdiag(lapply(model.data$trend[model.data$trend=="RW"], 
                            function(x){
                              if(x=="const"){return(matrix(rep(1,n.year), ncol=1))}
                              else return(cbind(rep(1,n.year),d.yrs))                         
                            }))
    } else{
      Xz.rw <- .bdiag(lapply(model.data$trend[model.data$trend=="RW"], 
                             function(x){
                               if(x=="const"){return(matrix(rep(1,n.year), ncol=1))}
                               else return(rep(1,n.year))                         
                             }))
    }
    A <- solve(crossprod(Xz.rw), t(Xz.rw))
  }
  # create design matrix for q fixed-only sites
  if(!is.null(rw.order$alpha)){
    alpha.order=rw.order$alpha
  } else{alpha.order=1}
  if(use.zi){
    Xq <- .bdiag(lapply(model.data$avail[model.data$avail!="none"], 
                            function(x){
                              if(x=="const"){return(matrix(rep(1,n.year), ncol=1))}
                              else return(cbind(rep(1,n.year),d.yrs))                         
                            }))
  }
  # Precision for epsilon
  if(!missing(sig.abund)) {
    Qeps <- Diagonal(x=1/log(1+(data[,sig.abund]/data[,abund.name])^2))
  } else {Qeps <- Diagonal(x=rep(1.0E8, length(y)))}
  # Precision for omega and alpha
  Q0.omega.s <- Matrix(iar.Q(n.year, omega.order))
  Q0.alpha.s <- Matrix(iar.Q(n.year, alpha.order))
  r0.omega.s <- dim(Q0.omega.s)[1] - omega.order
  r0.alpha.s <- dim(Q0.alpha.s)[1] - alpha.order
  if(use.omega) suppressMessages(Qomega.0 <- kronecker(Diagonal(sum(model.data$trend=="RW")), Q0.omega.s))
  if(use.alpha) suppressMessages(Qalpha.0 <- kronecker(Diagonal(sum(model.data$avail=="RW")), Q0.alpha.s))
  
  # PRIORS ###
  #gamma
  if(use.gam){
    if(is.null(prior.list$gamma)){
      Qg <- matrix(0,ncol(Xy))
      g0 <- rep(0,ncol(Xy))
    } else{
      Qg <- prior.list$gamma$Q.gamma
      g0 <- prior.list$gamma$gamma.0
    }
  }
  # beta
  if(is.null(prior.list$beta)){
    Qb <- matrix(0,ncol(Xz), ncol(Xz))
    b0 <- rep(0,ncol(Xz))
  } else{
    Qb <- prior.list$beta$Q.beta
    b0 <- prior.list$beta$beta.0
  }
  #tau
  if(is.null(prior.list$tau)){
    a.tau <- 0.5
    b.tau <- 5.0E-5
  } else{
    a.tau <- prior.list$tau$a.tau
    b.tau <- prior.list$tau$b.tau
  }
  #zeta
  if(is.null(prior.list$zeta)){
    a.zeta <- 0.5
    b.zeta <- 5.0E-5
  } else{
    a.zeta <- prior.list$zeta$a.zeta
    b.zeta <- prior.list$zeta$b.zeta
  }
  #phi
  if(use.alpha){
    if(is.null(prior.list$phi)){
      a.phi <- 0.5
      b.phi <- 5.0E-5
    } else{
      a.phi <- prior.list$phi$a.phi
      b.phi <- prior.list$phi$b.phi
    }
  }
  # theta
  if(use.zi){
    if(is.null(prior.list$theta)){
      Qt <- matrix(0,ncol(Xq), ncol(Xq))
      t0 <- rep(0,ncol(Xq))
    } else{
      Qt <- prior.list$theta$Q.theta
      t0 <- prior.list$theta$theta.0
    }
  }
  
  
  # INITIAL VALUES AND PRECALCULATIONS ###
  # gamma
  if(use.gam){
    g <- g0
    Xyg <- Xy%*%g
    XyQepsM <- crossprod(Xy, Qeps) %*% M
    D.g <- crossprod(Xy, Qeps) %*% Xy + Qg
    D.g.inv <- solve(D.g)
    cholD.g.inv <- t(chol(D.g.inv))
    XyQepsy.Qgg0 <- crossprod(Xy, Qeps) %*% y + Qg %*% g0
  } else Xyg <- 0
  #Call init making functions
  newdata <- data.frame(d.yrs)
  colnames(newdata) <- time.name
  
  if(use.zi){
    inits.zi <- suppressWarnings(getSiteZIInits(data.orig, abund.name, site.name, time.name, 
                          Q0.alpha.s, r0.alpha.s, a.phi, b.phi, newdata, model.data))
    a <- inits.zi$alpha
    theta <- inits.zi$theta
    phi <- inits.zi$phi
    muq <- as.vector(suppressMessages(Xq%*%theta))
    q <- as.vector(Xq%*%theta + a)
    if(use.alpha){
      Phi <- Diagonal(x=rep(phi, each=n.year))
      Qalpha <- Phi%*%Qalpha.0
      Ialpha <- Diagonal(nrow(Qalpha))
    }
    #z.01.tilde <- rep(1, length(qzi))
    z.01.tilde <- rtruncnorm(Nzi, a=a.zi[qzi], b=b.zi[qzi], mean=q, sd=1)
  }
  inits <- getSiteREInits(data, abund.name, site.name, time.name, 
                          ln.adj, Q0.omega.s, r0.omega.s, a.tau, b.tau, newdata, model.data)
  
  #beta
  b <- inits$b
  Qbb0 <- Qb%*%b
  muz <- suppressMessages(Xz %*% b)
  #tau and omega
  if(use.omega){
    tau <- inits$tau
    Tau <- Diagonal(x=rep(tau, each=n.year))
    Qomega <- suppressMessages(Tau %*% Qomega.0)
  }
  omega <- inits$omega
  #Qdelta and zeta
  zeta <- inits$zeta
  Qdelta <-  Diagonal(x=rep(zeta, each=n.year))
  #z 
  MQepsM <- crossprod(M,Qeps) %*% M
  MQeps <- crossprod(M,Qeps)
  D.z <- MQepsM + Qdelta
  d.z<- MQeps%*%(y-Xyg) + Qdelta%*%(muz + omega)
  z <- as.vector(solve(D.z, d.z))
  z.pred <- as.vector(muz + omega + rnorm(bigN, 0, 1/sqrt(Qdelta@x)))
  #aggregation
  ag.df <- expand.grid(yrs, levels(ag.data[,aggregation]))
  if(length(unique(ag.data[,aggregation]))>1) ag.mm <- model.matrix(~(ag.df[,2]+0) + (ag.df[,1]:ag.df[,2]+0))
  else ag.mm <- cbind(rep(1,nrow(ag.df)), yrs)
  ag.nms <- apply(ag.df, 1, paste, collapse="-")
  H <- solve(crossprod(ag.mm))%*%t(ag.mm)
  site.nms <- apply(expand.grid(yrs, unique(data[,site.name])), 1, paste, collapse="-")
    
  # STORAGE MATRICES ###
  #prediction and trends
  ag.pred.abund.stor <- matrix(nrow=iter, ncol=length(ag.nms))
  colnames(ag.pred.abund.stor) <- ag.nms
  ag.abund.stor <- matrix(nrow=iter, ncol=length(ag.nms))
  colnames(ag.abund.stor) <- ag.nms
  if(keep.site.abund){
    pred.site.abund.stor <- matrix(nrow=iter, ncol=length(site.nms))
    real.site.abund.stor <- matrix(nrow=iter, ncol=length(site.nms))
    colnames(pred.site.abund.stor) <- colnames(real.site.abund.stor) <- site.nms
  }
  ag.trend.stor <- matrix(nrow=iter, ncol=2*length(levels(data[,aggregation])))
  colnames(ag.trend.stor) <- c(paste(levels(data[,aggregation]),"(Intercept)",sep=":"), 
                                 paste(levels(data[,aggregation]), "(Trend)", sep=":"))
  ag.pred.trend.stor <- ag.trend.stor  
  #beta
  if(keep.site.param){
    b.stor <- matrix(nrow=iter, ncol=ncol(Xz))
    colnames(b.stor) <- unlist(mapply( 
                                      function(x,s){
                                        if(x!="const") return(paste(c("(Intercept)", "(Trend)"), c(s,s), sep=":"))
                                        else return(paste("(Intercept)",s, sep=":"))
                                      },
                                      x=model.data$trend, s=as.character(model.data[,site.name])
                              ))
  }
  #gamma
  if(keep.obs.param & use.gam){
    g.stor <- matrix(nrow=iter, ncol=ncol(Xy))
    colnames(g.stor) <- colnames(Xy)
  }
  #p and phi
  if(keep.site.param & use.zi){
    p.stor <- matrix(nrow=iter, ncol=n.year*sum(model.data$avail!="none")) 
    colnames(p.stor) <- apply(expand.grid(d.yrs, as.character(model.data[model.data$avail!="none",site.name])),1,paste,collapse=":")
    if(use.alpha){
      phi.stor <- matrix(nrow=iter, ncol=sum(model.data$avail=="RW"))
      colnames(phi.stor) <- as.character(model.data[model.data$avail=="RW", site.name])
    }
  }
  #tau and zeta
  if(keep.site.param){
    if(use.omega){
      tau.stor = matrix(nrow=iter, ncol=sum(model.data$trend=="RW"))
      colnames(tau.stor) <- as.character(model.data[model.data$trend=="RW", site.name])
    }
    zeta.stor <- matrix(nrow=iter, ncol=n.site)
    colnames(zeta.stor) <- as.character(levels(data.orig[,site.name]))
  }
  # MCMC UPDATES ###
  st.mcmc <- Sys.time()

  for(i in 1:(burn + iter*thin)){
    
    #update z
    D.z <- MQepsM + Qdelta
    d.z<- MQeps%*%(y-Xyg) + Qdelta%*%(muz + omega)
    if(!use.trunc){
      z <- as.vector(solve(D.z, d.z) + solve(chol(D.z), rnorm(bigN,0,1)))
    } else {
      z <- rtruncnorm(bigN, a=lower, b=upper, mean=as.vector(solve(D.z,d.z)), sd=1/sqrt(diag(D.z)))
    }
    
    #update alpha, phi, and z.01.tilde
    if(use.zi){
      #alpha, phi, and z.01.tilde
      if(use.alpha){
        D.alpha <- Ialpha + Qalpha
        d.alpha <- (z.01.tilde[qrw]-muq[qrw])
        a[qrw] <- as.vector(solve(D.alpha,d.alpha) + solve(chol(D.alpha), rnorm(Nqrw,0,1)))
        b1.phi <- b.phi + aggregate(a[qrw], list(site.idx[qzi][qrw]), FUN=function(x,Q0.s){crossprod(x, Q0.s)%*%x}, Q0.s=as.matrix(Q0.alpha.s))$x/2
        a1.phi <- a.phi + r0.alpha.s/2
        phi <- rgamma(length(b1.phi), a1.phi, b1.phi)
        Phi <- Diagonal(x=rep(phi, each=n.year))
        Qalpha <- Phi %*% Qalpha.0
      }
      #theta
      D.t <- crossprod(Xq) + Qt
      d.t <- suppressMessages(crossprod(Xq, z.01.tilde-a) + Qt%*%t0)
      theta <- solve(D.t,d.t) + solve(chol(D.t), rnorm(ncol(Xq),0,1))
      muq <- Xq%*%theta
      q <- as.vector(muq + a)
      #z.01.tilde
      z.01.tilde <- rtruncnorm(Nzi, a=a.zi[qzi], b=b.zi[qzi], mean=q, sd=1)
    }
    
    #update gamma
    if(use.gam){
      d.g <- XyQepsy.Qgg0 - XyQepsM%*%z
      g <- D.g.inv %*% d.g + cholD.g.inv %*% rnorm(ncol(Xy),0,1)
    }
  
    #update beta
    D.b <- suppressMessages(crossprod(Xz, Qdelta) %*% Xz + Qb)
    d.b <- crossprod(Xz,Qdelta)%*%(z-omega) + Qbb0
    b <- solve(D.b,d.b) + solve(chol(D.b), rnorm(ncol(Xz),0,1))
    muz <- Xz%*%b
    
    if(use.omega){
      #update omega
      D.omega <- Qdelta[omegarw,omegarw] + Qomega
      d.omega <- Qdelta[omegarw,omegarw]%*%(z[omegarw]-muz[omegarw])
      omega[omegarw] <- as.vector(solve(D.omega,d.omega) + solve(chol(D.omega), rnorm(sum(omegarw),0,1)))
      omega[omegarw] <- as.vector(omega[omegarw] - (D.omega%*%t(A)) %*% solve(A%*%D.omega%*%t(A)) %*% (A%*%omega[omegarw]))
      
      #update tau
      b1.tau <- b.tau + aggregate(as.vector(omega[omegarw]), list(site.idx[omegarw]), FUN=function(x,Q0.s){crossprod(x, Q0.s)%*%x}, Q0.s=as.matrix(Q0.omega.s))$x/2
      a1.tau <- a.tau + r0.omega.s/2
      tau <- rgamma(sum(model.data$trend=="RW"), a1.tau, b1.tau)
      Tau <- Diagonal(x=rep(tau, each=n.year))
      Qomega <- Tau %*% Qomega.0
    }
    
    #update zeta
    res.z2 <- as.vector((z-muz-omega)^2)
    b1.zeta <- b.zeta + aggregate(res.z2, list(site.idx), FUN=sum)$x/2
    a1.zeta <- a.zeta + n.year/2
    zeta <- rgamma(n.site, a1.zeta, b1.zeta)
    Qdelta <-  Diagonal(x=rep(zeta, each=n.year))
        
    #aggregate abundence trends
    N.obs <- exp(as.vector(z))-ln.adj
    if(use.zi){
      if(incl.zeros) N.obs[qzi] <- ifelse(z.01.tilde>0, 1, 0)*N.obs[qzi]
    }
    ag.abund <- aggregate(N.obs, list(yr.idx, agg.idx), FUN=sum)
    ag.abund <- log(ag.abund[ag.abund[,1]>=start & ag.abund[,1]<=end,"x"] + ln.adj)
    ag.trend <- H%*%ag.abund
    
    #posterior predictive trend
    if(!use.trunc) {
      z.pred <- as.vector(muz + omega + rnorm(bigN, 0, 1/sqrt(Qdelta@x)))
    } else z.pred <- rtruncnorm(bigN, a=lower, b=upper, mean=as.vector(muz+omega), sd=1/sqrt(Qdelta@x))
    N.pred <- exp(z.pred)-ln.adj
    if(use.zi){
      if(incl.zeros) N.pred[qzi] <- ifelse(rnorm(sum(qzi), mean=q, sd=1)>0,1,0)*N.pred[qzi]
    }
    ag.abund.pred <- aggregate(N.pred, list(yr.idx, agg.idx), FUN=sum)
    ag.abund.pred <- log(ag.abund.pred[ag.abund.pred[,1]>=start & ag.abund.pred[,1]<=end,"x"] + ln.adj)
    ag.pred.trend <- H%*%ag.abund.pred
    
    #Store MCMC output
    if(i > burn){
      if((i-burn)%%thin==0){
        j <- (i-burn)/thin
        ag.pred.abund.stor[j,] <- pmax(0,exp(as.vector(ag.abund.pred))-ln.adj) 
        ag.abund.stor[j,] <- pmax(0,exp(as.vector(ag.abund))-ln.adj) 
        ag.trend.stor[j,] <- as.vector(ag.trend)
        ag.pred.trend.stor[j,] <- as.vector(ag.pred.trend)
        if(keep.site.abund){
          pred.site.abund.stor[j,] <- N.pred[yr.idx>=start & yr.idx<=end]
          real.site.abund.stor[j,] <- N.obs[yr.idx>=start & yr.idx<=end]
        }
        if(keep.site.param){
          b.stor[j,] <- as.vector(b)
          if(use.omega) tau.stor[j,] <- tau
          zeta.stor[j,] <- zeta
        }
        if(keep.obs.param & use.gam) g.stor[j,] <- as.numeric(g)
        if(use.zi & keep.site.param){
          p.stor[j,] <- as.vector(pnorm(q))
          if(use.alpha) phi.stor[j,] <- phi
        }
      }
    }
    
    #Estimate time to completion
    if(i==15){
      tme <- Sys.time()
      tpi <- difftime(tme,st.mcmc, units="mins")/15
      ttc <- round((iter*thin + burn - 15)*tpi, 1)
      cat("\nApproximate time to completion ", ttc, "minutes... \n\n")  
      pb <- txtProgressBar(min = 0, max = burn + iter*thin, style = 3)
    }
    
    if(i>15) setTxtProgressBar(pb, i)
    
  }
  close(pb)
  if(keep.site.abund){
    pred.site.abund <- mcmc(pred.site.abund.stor)
    real.site.abund <- mcmc(real.site.abund.stor)
  }
  else pred.site.abund <- NULL
  if(keep.site.param){
    site.param=list(
      beta=mcmc(b.stor),
      zeta=mcmc(zeta.stor)
    )
    if(use.omega) site.param = append(site.param, list(tau=mcmc(tau.stor)))
  }

  else site.param=NULL
  if(keep.obs.param & use.gam){
    gamma=mcmc(g.stor)
  }
  else gamma=NULL
  if(use.zi) prob.avail <- mcmc(p.stor)
  else  prob.avail <- NULL
  mcmc.sample <- list(pred.trend=mcmc(ag.pred.trend.stor), trend=mcmc(ag.trend.stor), 
                      aggregated.pred.abund=mcmc(ag.pred.abund.stor), aggregated.real.abund=mcmc(ag.abund.stor),
                      pred.site.abund=pred.site.abund, real.site.abund=real.site.abund,
                      site.param=site.param, gamma=gamma,
                      prob.avail=prob.avail)
  summ.dat1 <- expand.grid(d.yrs, levels(ag.data[,2]))
  summ.dat1 <- summ.dat1[summ.dat1[,1]>=start & summ.dat1[,1]<=end,]
  summ.dat3 <- summ.dat1
  
  summ.dat1$post.median.abund <- apply(mcmc.sample$aggregated.pred.abund, 2, median)
  summ.dat1$low90.hpd <- HPDinterval(mcmc.sample$aggregated.pred.abund, 0.9)[,1]
  summ.dat1$hi90.hpd <- HPDinterval(mcmc.sample$aggregated.pred.abund, 0.9)[,2]
  
  summ.dat3$post.median.abund <- apply(mcmc.sample$aggregated.real.abund, 2, median)
  summ.dat3$low90.hpd <- HPDinterval(mcmc.sample$aggregated.real.abund, 0.9)[,1]
  summ.dat3$hi90.hpd <- HPDinterval(mcmc.sample$aggregated.real.abund, 0.9)[,2]
  
  if(keep.site.abund){
    summ.dat2 <- expand.grid(d.yrs, unique(as.character(data[,site.name])))
    summ.dat2 <- summ.dat2[summ.dat2[,1]>=start & summ.dat2[,1]<=end,]
    summ.dat2$post.median.abund <- apply(mcmc.sample$pred.site.abund, 2, median)
    summ.dat2$low90.hpd <- HPDinterval(mcmc.sample$pred.site.abund, 0.9)[,1]
    summ.dat2$hi90.hpd <- HPDinterval(mcmc.sample$pred.site.abund, 0.9)[,2]
  }
  else summ.dat2 <- NULL

  colnames(summ.dat1)[1:2] <- c(time.name, aggregation)
  colnames(summ.dat3)[1:2] <- c(time.name, aggregation)
  colnames(summ.dat2)[1:2] <- c(time.name, site.name)
  output <- list(
    trend.summary=summary(mcmc.sample$pred.trend), 
    aggregation.pred.summary=summ.dat1, 
    aggregation.real.summary=summ.dat3,
    site.summary=summ.dat2, 
    mcmc.sample=mcmc.sample, 
    original.data=data.orig
    )
  attr(output, "site.name") <- site.name
  attr(output, "time.name") <- time.name
  attr(output, "aggregation") <- aggregation
  attr(output, "ln.adj") <- ln.adj
  
  return(output)
  
}

##############################################################
#' @title Recalculate posterior predictive trend with a different time scale.
#' @param x An mcmc augmentation object produced by a call to \code{\link{mcmc.aggregate}} or 
#' an element of the list produced by a call to \code{\link{newAggregation}}.
#' @param start A new start value for the time span
#' @param end A new end value for the time span
#' @param type The type of trend calculated. Use \code{"pred"} for posterior predictive trends
#' and \code{"real"} to use the estimated, realized abumndance aggregation.
#' @param order The order of trend calculated. Can be one of \code{"lin"}, for linear trends, 
#' or, \code{"const"}, for mean log-abundence.
#' @export
#'  
updateTrend <- function(x, start, end, type="pred", order="lin"){
  #require(coda)
  if(type=="pred" | is.null(x$mcmc.sample$aggregated.real.abund)) smp <- x$mcmc.sample$aggregated.pred.abund
  else if(type=="real") smp <- x$mcmc.sample$aggregated.real.abund
  else stop("Unknown 'type', must be 'pred' or 'real'.\n")
  nms <- strsplit(colnames(smp),"-")
  time <- as.numeric(sapply(nms, function(x)x[[1]]))
  start <- max(start, min(time))
  end <- min(end, max(time))
  time.idx <- time>=start & time<=end
  agg <- sapply(nms, function(x)x[[2]])
  if(any(is.na(suppressWarnings(as.numeric(agg))))) agg <- factor(agg)
  else agg <- factor(as.numeric(agg))
  nms.agg <- as.character(levels(agg))
  
  yrs <- c(start:end)
  ag.df <- expand.grid(yrs, levels(agg))
  if(order=="lin"){
    if(length(nms.agg)>1) {
      ag.mm <- model.matrix(~(ag.df[,2]+0) + (ag.df[,1]:ag.df[,2]+0))
    } else {
      ag.mm <- cbind(rep(1,nrow(ag.df)), yrs)
    }
    H <- solve(crossprod(ag.mm))%*%t(ag.mm)
    y <- log(smp[,time.idx] + attr(x, "ln.adj"))
    tsmp <- t(apply(y, 1, function(y){H%*%y}))
    colnames(tsmp) <- c(paste(nms.agg, "(Intercept)"), paste(nms.agg, "(Trend)"))
  } else if(order=="const"){
    if(length(nms.agg)>1) ag.mm <- model.matrix(~(ag.df[,2]+0))
    else ag.mm <- matrix(rep(1,nrow(ag.df)), ncol=1)
    H <- solve(crossprod(ag.mm))%*%t(ag.mm)
    y <- log(smp[,time.idx] + attr(x, "ln.adj"))
    tsmp <- matrix(t(apply(y, 1, function(y){H%*%y})), ncol=ncol(ag.mm))
    colnames(tsmp) <- c(paste(nms.agg, "(Intercept)"))
  }
  else stop("Unknown 'type' specified! See ?updateTrend.")
  return(mcmc(tsmp))
}

##############################################################
#' @title Recalculate new site aggregations from a previous MCMC aggregation
#' 
#' @description 
#' If the site abundence sample was retained in a call to \code{\link{mcmc.aggregate}} the 
#' sites can be re-aggregated according to different region specifications.
#' 
#' @param fit The output list from a previous call to \code{\link{mcmc.aggregate}}.
#' In order to use this function, \code{keep.site.abund = TRUE} had to be used in the original creation of \code{fit}. Else,
#' there is nothing to be aggregated!
#' @param aggregation.data  A data frame with the sites in one column 
#' (with the same name as \code{site.name} used in the original call to create \code{fit}). The other columns
#' are factor variables defining other site aggregations.
#' @param type Which site abundance augmentation should be used, \code{"pred"} for posterior
#' predictive or \code{"real"} for realized (just the posterior).
#' 
#' @return 
#' A named list with names equal to the variables in \code{aggregation.data}. Each
#' element of the list is another list with elements:
#' \item{aggregated.abund}{The MCMC sample of the new aggregation}
#' \item{aggregation.summary}{A summary of the aggregation MCMC}
#' @export
#'  
newAggregation <- function(fit, aggregation.data, type="pred"){
  #require(coda)
  if(type == "pred") xxx <- fit$mcmc.sample$pred.site.abund
  if(type == "real") xxx <- fit$mcmc.sample$real.site.abund
  if(is.null(xxx)) stop("Site abundance data was not retained in the call to 'mcmc.aggreation()'\n Please re-run with 'keep.site.abund=TRUE'\n")
  site.name <- attr(fit,"site.name")
  time.name <- attr(fit,"time.name")
  site.idx <- data.frame(sapply(strsplit(colnames(xxx),"-"), function(x){paste(x[-1],collapse="-")}),
                    as.numeric(sapply(strsplit(colnames(xxx),"-"), function(x){x[[1]]}))) 
  colnames(site.idx) <- c(site.name, time.name)
  m1 <- merge(site.idx, aggregation.data, all=TRUE)
  ag.names <- colnames(aggregation.data)[colnames(aggregation.data)!=site.name]
  outlist <- vector("list", length(ag.names))
  #Tmat <- cbind(rep(1,length(unique(site.idx[,2]))), unique(site.idx[,2]))
  names(outlist) <- ag.names
    for(i in 1:length(ag.names)){
      cat("Processing", type, "aggregation:", ag.names[i], "...\n")
      a1 <- mcmc(t(apply(as.matrix(xxx), 1, FUN=function(v){aggregate(v, list(m1[,time.name],m1[,ag.names[i]]), FUN=sum)$x})))
      colnames(a1) <- apply(expand.grid(unique(site.idx[,time.name]), levels(factor(aggregation.data[,ag.names[i]]))),1,paste, collapse="-")
      ag.summary <- expand.grid(unique(site.idx[,time.name]), levels(factor(aggregation.data[,ag.names[i]])))
      colnames(ag.summary) <- c(time.name, ag.names[i])
      ag.summary$post.median.abund <- apply(a1,2,median)
      hpd <- HPDinterval(a1, 0.9)
      ag.summary$low90.hpd <- hpd[,1]
      ag.summary$hi90.hpd <- hpd[,2]
      if(type=="pred") {
        outlist[[i]] <- list(mcmc.sample=list(aggregated.pred.abund=a1), aggregation.pred.summary=ag.summary)
      } else  {
        outlist[[i]] <- list(mcmc.sample=list(aggregated.real.abund=a1), aggregation.real.summary=ag.summary)
      }
      attr(outlist[[i]], "ln.adj") <- attr(fit, "ln.adj")
    }  
  return(outlist)
}

