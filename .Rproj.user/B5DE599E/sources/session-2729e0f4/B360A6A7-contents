## Computes log-likelihood for a single study
dc_loglik <- function(E.new,W_s,P,M_s,A_s,fac) {
  term <- M_s*log(W_s*P%*%A_s%*%E.new)
  term[M_s==0&is.na(term)] <- 0
  return(-1*sum(W_s*(P%*%A_s%*%E.new)-term+
                  fac[M_s+1]))
}

## Computes log-likelihood across all studies
dc_loglik.all <- function(P,A,E,W,M,fac,S) {
  llk <- sum(sapply(1:S,function(s) dc_loglik(E[[s]],W[[s]],P,M[[s]],diag(A[s,]),fac)))
  return(llk)
}

## Sample one column of A
dc_sample_A_marg <- function(n,A,M,W,P,E,p,S,fac,gamma) {
  # For each study...
  for (s in 1:S) {

    if (gamma[1]==(-1)) {
      # Create a version of A where study contains the signature, and where it doesn't
      Ar1 <- A
      Ar0 <- A
      Ar1[s,n] <- 1
      Ar0[s,n] <- 0

      # Compute likelihoods under both scenarios
      p1 <- dc_loglik(E[[s]],W[[s]],P,M[[s]],diag(Ar1[s,]),fac)+log(p[s,n])
      p0 <- dc_loglik(E[[s]],W[[s]],P,M[[s]],diag(Ar0[s,]),fac)+log(1-p[s,n])
      prob <- exp(p1-sumLog(c(p1,p0)))

      # Draw according to these probabilities
      A[s,n] <- rbinom(1,1,prob)
    } else {
      # Propose a value from a beta distribution
      prop <- rbeta(1,gamma,3*gamma)
      Ar1 <- A
      Ar0 <- A
      Ar0[s,n] <- prop

      # Compute acceptance ratio
      p0 <- dc_loglik(E[[s]],W[[s]],P,M[[s]],diag(Ar0[s,]),fac)
      p1 <- dc_loglik(E[[s]],W[[s]],P,M[[s]],diag(Ar1[s,]),fac)
      ratio <- (exp((p0-p1)))
      if (runif(1)<ratio) {
        A[s,n] <- prop
      }
    }

  }
  return(A)
}

## Sample latent count
dc_sample_zijs <- function(i,j,s,M_s,P,E_s,A,N) {
  probs <- P[i,]*A[s,]*t(E_s[,j])
  if (sum(is.na(probs))==0) {
    return(rmultinom(1,M_s[i,j],prob=probs))
  } else {
    return(rep(0,length(probs)))
  }
}

## Sample latent count (vectorized)
dc_sample_zjs <- function(j,s,M_s,P,E_s,N) {
  probs <- sweep(P,2,E_s[,j],'*')
  return(sapply(1:nrow(probs),function(x) rmultinom(1,M_s[x,j],probs[x,])))
}

## Sample signature value (vectorized)
dc_sample_pn <- function(n,W,E,A,Z,alpha_p,beta_p,S,G,gamma) {
  shape <- alpha_p[,n]+1+rowSums(sapply(1:S,function(s) gamma*A[s,n]*rowSums(Z[[s]][[n]])))
  rate <- beta_p[,n]+gamma*sum(sapply(1:S,function(s) sum(A[s,n]*E[[s]][n,]*W[[s]][1,1])))
  return(rgamma(nrow(alpha_p),shape,rate))
}

## Sample exposure value
dc_sample_ens <- function(n,s,W_s,P,A,Z_s,alpha_e,beta_e,gamma) {
  shape <- alpha_e[[s]][n,]+gamma*A[s,n]*colSums(Z_s[[n]])+1
  rate <- beta_e[[s]][n,]+gamma*sum(P[,n]*A[s,n]*W_s[1,1])
  return(rgamma(ncol(alpha_e[[s]]),shape,rate))
}

## Sample signature hyperparameter beta
dc_sample_betap <- function(P,alpha_p,a_p,b_p,gamma) {
  return(array(rgamma(nrow(alpha_p)*ncol(alpha_p),a_p+gamma*alpha_p+gamma,b_p+gamma*P),
               dim=dim(alpha_p)))
}

## Posterior density for signature hyperparameter alpha
dc_alphainp_density <- function(i,n,P,y,beta_p,lambda_p,gamma) {
  return(log(lambda_p)-lambda_p*y+gamma*((y+1)*log(beta_p[i,n])+y*log(P[i,n])-lgamma(y+1)))
}

## Gamma proposal density
dc_gamma_density <- function(y,x) {
  return(dgamma(y,x,1,log=TRUE))
}

## Sample signature hyperparameter alpha
dc_sample_alphainp <- function(i,n,P,alpha_p,beta_p,lambda_p,gamma) {
  x <- alpha_p[i,n]
  y <- rgamma(1,x,1)
  alpha <- dc_alphainp_density(i,n,P,y,beta_p,lambda_p,gamma)-
    dc_alphainp_density(i,n,P,x,beta_p,lambda_p,gamma)+
    dc_gamma_density(x,y)-dc_gamma_density(y,x)
  if (runif(1) <= min(1,exp(alpha))) {
    return(y)
  } else {
    return(x)
  }
}

## Posterior density for exposure hyperparameter alpha
dc_alphanjse_density <- function(n,s,E_s,y,beta_e,lambda_e,gamma) {
  return(log(lambda_e)-lambda_e*y+sum(gamma*((y+1)*log(beta_e[[s]][n,]))+
                                        (y*log(E_s[n,]))-lgamma(y+1)))
}

## Sample exposure hyperparameter alpha
dc_sample_alphanse <- function(n,s,E_s,alpha_e,beta_e,lambda_e,gamma) {
  x <- alpha_e[[s]][n,1]
  y <- rgamma(1,x,1)
  alpha <- dc_alphanjse_density(n,s,E_s,y,beta_e,lambda_e,gamma)-
    dc_alphanjse_density(n,s,E_s,x,beta_e,lambda_e,gamma)+
    dc_gamma_density(x,y)-dc_gamma_density(y,x)
  if (runif(1) <= min(1,exp(alpha))) {
    return(y)
  } else {
    return(x)
  }
}

## Sample mu0
dc_sample_mu0 <- function(s,n,G,E,A,alpha_e,beta_e,cov,delta,mu0_alpha,mu0_beta,gamma=1) {
  return(rgamma(1,mu0_alpha+sum(alpha_e[[s]][n,]+1),
                mu0_beta+sum(sapply(1:G[s],function(j)
                  E[[s]][n,j]*exp(-1*t(delta[[s]][,n,drop=F])%*%cov[[s]][,j,drop=F])))))
}

## Sample indicators
dc_sample_indic <- function(delta,tau,c) {
  array(rbinom(length(delta),1,(delta^2/(tau^2*c^2))*dnorm(delta,0,c*tau)/
           ((delta^2/(tau^2*c^2))*dnorm(delta,0,c*tau)+dnorm(delta,0,tau))),dim=dim(delta))
}

## Sample truncation points
dc_sample_lambda <- function(delta,tau,c) {
  array(runif(length(delta),0,delta^2/(tau^2*c^2)),dim=dim(delta))
}

## Posterior for delta
dc_delta_f <- function(x,s,l,n,a,tau,E,cov,mu0,alpha_e,A,delta) {
  delta[[s]][l,n] <- x
  return(sum(-(alpha_e[[s]][n,]+1)*t(delta[[s]][,n,drop=F])%*%cov[[s]]-
           (E[[s]][n,]*mu0[[s]][n,]/(exp(t(delta[[s]][,n,drop=F])%*%cov[[s]]))))-
           (1/(2*a[[s]][l,n]^2*tau^2))*x^2)
}

### Full Sampler
discovery_covariates <- function(M,cov,lambda_p=0.5,a_p=5,b_p=0.05,lambda_e=0.1,a_e=5,b_e=0.01,
                                 Ntot=50,inits=NULL,fixed=F,A.fixed=NULL,first=T) {
  # Identify values
  K <- dim(M[[1]])[1]
  S <- length(M)
  G <- sapply(1:S,function(s) dim(M[[s]])[2])
  fac <- c(0,cumsum(log(seq_len(max(sapply(1:S,function(s) max(M[[s]])))))))
  W <- lapply(M,function(x) matrix(rep(colSums(x),K),nrow=K,ncol=ncol(x),byrow=T))
  for (i in 1:length(W)) {
    W[[i]][W[[i]]>0] <- median(W[[i]])
  }
  lambda_e <- rep(lambda_e,S)
  a_e <- rep(a_e,S)
  b_e <- rep(b_e,S)
  mu0_alpha <- a_e[1]
  mu0_beta <- b_e[1]

  # Initialize signature indicator matrix
  if (!fixed) {
    A <- array(0.5,dim=c(S,Ntot))
  } else {
    A <- A.fixed
  }
  N <- ncol(A)

  # Initialize Bernoulli parameters
  p <- array(1/4,dim=c(S,N))

  # Initialize hyperparameters
  alpha_p <- array(rexp(K*N,lambda_p),dim=c(K,N))
  beta_p <- array(rgamma(K*N,a_p,b_p),dim=c(K,N))
  alpha_e <- lapply(1:S,function(s) array(rexp(N*G[s],lambda_e[s]),dim=c(N,G[s])))
  mu0 <- lapply(1:S,function(s) array(1,dim=c(N,G[s])))
  delta <- lapply(1:S,function(s) array(0,dim=c(nrow(cov[[1]]),N)))
  indic <- lapply(1:S,function(s) array(0,dim=c(nrow(cov[[1]]),N)))
  tau <- 0.3
  c <- 10
  an <- lapply(1:S,function(s) array(1,dim=c(nrow(cov[[1]]),N)))
  lambda <- lapply(1:S,function(s) array(runif(length(delta[[s]]),0,delta[[s]]^2/(tau*c)),
                                         dim=dim(delta[[s]])))
  beta_e <- lapply(1:S, function(s) mu0[[s]]*exp(-1*t(delta[[s]])%*%cov[[s]]))

  # Initialize signatures and exposures
  P <- array(rgamma(K*N,alpha_p+1,beta_p),dim=c(K,N))
  P <- as.matrix(P)
  E <- lapply(1:S,function(s) array(rgamma(N*G[s],1,1),dim=c(N,G[s])))
  if (!is.null(inits)) {
    P <- inits[[1]]
    E <- inits[[2]]
  }

  # Initialize latent counts
  Z <- lapply(1:S,function(s) lapply(1:N,function(n) array(0,dim=c(K,G[s]))))
  for (s in 1:S) {
    for (i in 1:K) {
      for (j in 1:G[s]) {
        z <- dc_sample_zijs(i,j,s,M[[s]],P,E[[s]],A,N)
        for (n in 1:N) {
          Z[[s]][[n]][i,j] <- z[n]
        }
      }
    }
  }

  # Initialize chains
  A_chain <- list()
  P_chain <- list()
  E_chain <- list()
  indic_chain <- list()
  delta_chain <- list()

  if (first) {
    temps <- c(10^(rep(seq(-2,0,length.out=68),each=100)),rep(1,6200))
    params <- c(rep(seq(1,0.20,by=-0.05),each=100),
                rep(seq(0.25,0.01,length.out=50),each=100),
                rep(seq(0.005,0.0005,length.out=10),each=100),
                rep(-1,5300))
    pb1 = txtProgressBar(min = 0, max = 7700, initial = 0, style=3)
    pb2 = txtProgressBar(min = 7701, max = 13000, initial = 0, style=3)
  } else {
    params <- rep(-1,5000)
    pb = txtProgressBar(min = 0, max = 5000, initial = 0, style=3)
  }

  gamma <- 1
  for (l in 1:length(params)) {

    if (l==1 & first) {
      cat("\nTempering phase:\n")
    }
    if (l==7701 & first) {
      cat("\nSampling phase:\n")
    }
    if (l==1 & !first) {
      cat("\nEvaluating solution:\n")
    }

    # Update exposures
    for (s in 1:S) {
      for (n1 in sample(1:N,N,replace=FALSE)) {
        E[[s]][n1,] <- dc_sample_ens(n1,s,W[[s]],P,A,Z[[s]],alpha_e,beta_e,gamma)
      }
    }

    # Update signatures
    if (params[l]==(-1)) {
      for (n in 1:N) {
        P[,n] <- dc_sample_pn(n,W,E,A,Z,alpha_p,beta_p,S,G,gamma)
      }
    } else {
      for (n1 in sample(1:N,N,replace=FALSE)) {
        maxit <- 1
        curr <- Inf
        while (curr > 0.9 & maxit < 10) {
          P[,n1] <- dc_sample_pn(n1,W,E,A,Z,alpha_p,beta_p,S,G,temps[l])
          curr <- max(cosineDist(t(P[,n1,drop=F]),t(P[,-n1])))
          maxit <- maxit+1
        }
      }
    }

    # Update columns of A
    if (!fixed) {
      sigs.order <- sample(1:N,N,replace=FALSE)
      for (n in sigs.order) {
        A <- dc_sample_A_marg(n,A,M,W,P,E,p,S,fac,params[l])
      }
    }

    # Update latent counts
    if (params[l]==(-1)) {
      for (s in 1:S) {
        nonzero.sigs <- which(A[s,]==1)
        for (j in 1:G[s]) {
          nonzero <- which(M[[s]][,j]!=0)
          if (length(nonzero)>0) {
            z <- dc_sample_zjs(j,s,M[[s]][nonzero,,drop=F],P[nonzero,nonzero.sigs,drop=F],
                               E[[s]][nonzero.sigs,,drop=F],N)
            if (length(nonzero.sigs)==1) {
              Z[[s]][[nonzero.sigs]][nonzero,j] <- z
            } else {
              for (n1 in 1:length(nonzero.sigs)) {
                Z[[s]][[nonzero.sigs[n1]]][nonzero,j] <- z[n1,]
              }
            }
          }
        }
      }
    } else {
      for (s in 1:S) {
        for (i in 1:K) {
          for (j in 1:G[s]) {
            z <- dc_sample_zijs(i,j,s,M[[s]],P,E[[s]],A,N)
            for (n in 1:N) {
              Z[[s]][[n]][i,j] <- z[n]
            }
          }
        }
      }
    }

    # Update exposure hyperparameters alpha
    for (s in 1:S) {
      for (n in 1:N) {
        alpha_e[[s]][n,] <- dc_sample_alphanse(n,s,E[[s]],alpha_e,beta_e,lambda_e[s],gamma)
      }
    }

    # Update signature hyperparameters beta
    beta_p <- dc_sample_betap(P,alpha_p,a_p,b_p,gamma)

    # Update signature hyperparameters alpha
    for (i in 1:K) {
      for (n in 1:N) {
        alpha_p[i,n] <- dc_sample_alphainp(i,n,P,alpha_p,beta_p,lambda_p,gamma)
      }
    }

    # Update covariate-related parameters
    for (s in 1:S) {
      for (n in 1:N) {
        mu0[[s]][n,] <- dc_sample_mu0(s,n,G,E,A,alpha_e,beta_e,cov,delta,mu0_alpha,
                                      mu0_beta,gamma=1)
      }
    }

    for (s in 1:S) {
      lambda[[s]] <- dc_sample_lambda(delta[[s]],tau,c)
    }

    for (s in 1:S) {
      for (n in 1:N) {
        for (f in 1:nrow(cov[[1]])) {
          delta_f_inst <- function(x) { dc_delta_f(x,s,f,n,an,tau,E,cov,mu0,alpha_e,A,delta) }
          if (an[[s]][f,n]==1) {
            delta[[s]][f,n] <- armspp::arms(1,delta_f_inst,-10,10,metropolis=T)
          } else {
            thresh <- c(-1*sqrt(lambda[[s]][f,n]*tau^2*c^2),sqrt(lambda[[s]][f,n]*tau^2*c^2))
            delta[[s]][f,n] <- armspp::arms(1,delta_f_inst,-10,10,metropolis=T)
            cnter <- 0
            while ((!(as.numeric(delta[[s]][f,n])<thresh[1]|as.numeric(delta[[s]][f,n])>thresh[2]))) {
              delta[[s]][f,n] <- armspp::arms(1,delta_f_inst,-10,10,metropolis=T)
              cnter <- cnter+1
              if (cnter==10) {
                const <- max(c(delta_f_inst(thresh[1]),delta_f_inst(thresh[2])))
                delta_f_inst2 <- function(x) { sapply(x,function(val)
                  exp(-1*const+dc_delta_f(val,s,f,n,an,tau,E,cov,mu0,alpha_e,A,delta))) }
                comp1 <- integrate(delta_f_inst2,thresh[2],10)$value
                comp2 <- integrate(delta_f_inst2,-10,thresh[1])$value
                pick <- rbinom(1,1,comp1/(comp1+comp2))
                if (pick==1) {
                  delta[[s]][f,n] <- as.numeric(armspp::arms(1,delta_f_inst,thresh[2],10,
                                                             metropolis=T))
                } else {
                  delta[[s]][f,n] <- as.numeric(armspp::arms(1,delta_f_inst,-10,thresh[1],
                                                             metropolis=T))
                }
              }
            }
          }
          delta[[s]] <- array(as.numeric(delta[[s]]),dim=c(nrow(cov[[1]]),N))
        }
      }
      delta[[s]] <- array(as.numeric(delta[[s]]),dim=c(nrow(cov[[1]]),N))
    }

    for (s in 1:S) {
      indic[[s]] <- dc_sample_indic(delta[[s]],tau,c)
      an[[s]] <- array(ifelse(indic[[s]]==0,1,c),dim=dim(indic[[s]]))
    }
	  beta_e <- lapply(1:S, function(s) mu0[[s]]*exp(-1*t(delta[[s]])%*%cov[[s]]))

    # Store in chains
    A_chain[[l]] <- A
    P_chain[[l]] <- P
    E_chain[[l]] <- E
    indic_chain[[l]] <- indic
    delta_chain[[l]] <- delta

    if (first) {
      setTxtProgressBar(pb1,l)
      setTxtProgressBar(pb2,l)
    } else {
      setTxtProgressBar(pb,l)
    }
  }

  if (first) {
    close(pb1)
    close(pb2)
    As <- sapply(7701:length(A_chain),function(x)
      paste0(sort(sapply(1:N,function(y) bitsToInt(A_chain[[x]][,y])),decreasing=T),
             collapse=""))
    As.tab <- sort(table(As),decreasing=TRUE)
    As.vec <- names(As.tab)[1:min(5,sum(As.tab>1))]
    inds <- sapply(As.vec,function(x) which(As==x)[1]+7700)
    As.list <- lapply(inds,function(x) A_chain[[x]][,colSums(A_chain[[x]])>0])
    Ps.list <- lapply(inds,function(x) P_chain[[x]][,colSums(A_chain[[x]])>0])
    Es.list <- lapply(inds,function(x) lapply(1:S,function(s) E_chain[[x]][[s]][colSums(A_chain[[x]])>0,]))
    return(list(As.list,Ps.list,Es.list))
  } else {
    close(pb)
    liks <- sapply(1001:5000,function(x) dc_loglik.all(P_chain[[x]],A_chain[[x]],
                                                    E_chain[[x]],W,M,fac,S))
    bic <- (nrow(P)*ncol(P)+sum(rowSums(A.fixed)*G))*log(nrow(P)*sum(G))-2*max(liks)

    colnames(A.fixed) <- paste0('sig',1:ncol(A.fixed))
    rownames(A.fixed) <- paste0('study',1:nrow(A.fixed))

    Pest <- array(0,dim=c(dim(P),4000))
    for (i in 1001:5000) {
      Pest[,,(i-1000)] <- P_chain[[i]]
    }
    Pest <- apply(Pest,c(1,2),median)
    colnames(Pest) <- colnames(A.fixed)

    Eest <- list()
    for (s in 1:S) {
      Eest.s <- array(0,dim=c(dim(E[[s]][which(A[s,]==1),,drop=F]),4000))
      for (i in 1001:5000) {
        Eest.s[,,(i-1000)] <- E_chain[[i]][[s]][which(A[s,]==1),,drop=F]
      }
      Eest[[s]] <- apply(Eest.s,c(1,2),median)
      rownames(Eest[[s]]) <- colnames(A.fixed)[A[s,]==1]
    }
    names(Eest) <- rownames(A.fixed)

    indics <- lapply(1:S,function(s) {
      out <- t(sapply(1:nrow(cov[[1]]),function(f)
        rowMeans(sapply(1001:5000,function(x) indic_chain[[x]][[s]][f,]))))
      out <- out[,which(A[s,]==1),drop=F]
      colnames(out) <- colnames(A.fixed)[A[s,]==1]
      out
      })
    names(indics) <- rownames(A.fixed)

    deltas <- lapply(1:S,function(s) {
      out <- t(sapply(1:nrow(cov[[1]]),function(f)
        rowMedians(sapply(1001:5000,function(x)
          ifelse(indic_chain[[x]][[s]][f,]==1,delta_chain[[x]][[s]][f,],NA)))))
      out <- out[,which(A[s,]==1),drop=F]
      colnames(out) <- colnames(A.fixed)[A[s,]==1]
      out
    })
    names(deltas) <- rownames(A.fixed)

    return(list(BIC=bic,A=A.fixed,P=Pest,E=Eest,indics=indics,coefs=deltas,W=W))
  }
}

#' Run the Discovery-Only Sampler with Covariates
#'
#' This runs the fully unsupervised discovery-only sampler for multi-study NMF, with covariates.
#'
#' @param M List of matrices corresponding to each study, in features-by-samples format.
#' @param cov List of covariates matrices corresponding to each study,
#' in covariates-by-samples format; the same covariates in the same order should be supplied for
#' each study. Covariates must be numeric.
#' @param lambda_p Lambda_p, set to 0.5 by default.
#' @param a_p a_p, set to 5 by default.
#' @param b_p b_p, set to 0.05 by default.
#' @param lambda_e lambda_e, set to 0.1 by default.
#' @param a_e a_e, set to 5 by default.
#' @param b_e b_e, set to 0.01 by default.
#' @param Ntot Total possible number of signatures to consider, set to 50 by default.
#' @return A named list describing the best solution:
#' \itemize{
#'  \item{BIC: }{the BIC value of the solution.}
#'  \item{A: }{the point estimate of the A matrix, which indicates whether each signature is
#'  present in each study.}
#'  \item{P: }{the posterior median of the signatures matrix.}
#'  \item{E: }{a list of the posterior medians of the exposures matrices for each study s,
#'  subsetted to just the signatures indicated by A to belong to that study.}
#'  \item{indics: }{a list of the posterior inclusion probabilities for each study s,
#'  describing the probability that each covariate affects each signature indicated by A
#'  to belong to that study.}
#'  \item{coefs: }{a list of the coefficients, conditional on inclusion, for each study s,
#'  describing the contribution of each covariate to each signature indicated by A to belong
#'  to that study.}
#'  \item{W: }{a list of the normalizing constants for each study; note that the input counts
#'  matrix M_s for each study s can be modeled as Poisson(W_s * P_s E_s), where P_s and E_s
#'  contain just the signatures that belong to study s.}
#'  }
#' @export
discovery_cov_fit <- function(M,cov,lambda_p=0.5,a_p=5,b_p=0.05,lambda_e=0.1,a_e=5,b_e=0.01,
                              Ntot=50) {

  initial <- discovery_covariates(M,cov,lambda_p,a_p,b_p,lambda_e,a_e,b_e,
                                  Ntot,inits=NULL)
  outs <- lapply(1:length(initial[[1]]),function(x)
    discovery_covariates(M,cov,lambda_p,a_p,b_p,lambda_e,a_e,b_e,
                      inits=list(initial[[2]][[x]],initial[[3]][[x]]),
                      fixed=T,A.fixed=initial[[1]][[x]],first=F))
  pick <- which.min(sapply(outs,function(x) x[[1]]))
  return(outs[[pick]])
}
