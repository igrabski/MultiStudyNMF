## Compute log-likelihood of one study
rdc_loglik <- function(s,M_s,W_s,Pr,Er_s,Ar,P,E_s,A,fac) {
  term <- M_s*log(W_s*P%*%diag(A[s,])%*%E_s+W_s*Pr%*%diag(Ar[s,])%*%Er_s)
  term[M_s==0&is.na(term)] <- 0
  return(-1*sum(W_s*(P%*%diag(A[s,])%*%E_s+Pr%*%diag(Ar[s,])%*%Er_s)-term+
                  fac[M_s+1]))
}

## Sample one column of A
rdc_sample_A <- function(n,A,M,W,P,E,p,S,N,fac,Ar,Pr,Er,gamma) {
  for (s in 1:S) {
    if (gamma==(-1)) {
      Ar1 <- A
      Ar0 <- A
      Ar1[s,n] <- 1
      Ar0[s,n] <- 0
      p1 <- rdc_loglik(s,M[[s]],W[[s]],Pr,Er[[s]],Ar,P,E[[s]],Ar1,fac)+log(p[s,n])
      p0 <- rdc_loglik(s,M[[s]],W[[s]],Pr,Er[[s]],Ar,P,E[[s]],Ar0,fac)+log(1-p[s,n])
      prob <- exp(p1-sumLog(c(p1,p0)))
      A[s,n] <- rbinom(1,1,prob)
    } else {
      scale <- 9
      prop <- rbeta(1,gamma,scale*gamma)
      Ar1 <- A
      Ar0 <- A
      Ar0[s,n] <- prop

      # Compute acceptance ratio
      p0 <- rdc_loglik(s,M[[s]],W[[s]],Pr,Er[[s]],Ar,P,E[[s]],Ar0,fac)
      p1 <- rdc_loglik(s,M[[s]],W[[s]],Pr,Er[[s]],Ar,P,E[[s]],Ar1,fac)
      ratio <- (exp((p0-p1)))
      if (runif(1)<ratio) {
        A[s,n] <- prop
      }
    }
  }
  return(A)
}

## Sample one column of Ar
rdc_sample_Ar <- function(n,Ar,M,W,Pr,Er,P,E,A,pr,S,fac,gamma) {
  for (s in 1:S) {
    if (gamma==(-1)) {
      Ar1 <- Ar
      Ar0 <- Ar
      Ar1[s,n] <- 1
      Ar0[s,n] <- 0
      p1 <- rdc_loglik(s,M[[s]],W[[s]],Pr,Er[[s]],Ar1,P,E[[s]],A,fac)+log(pr[s,n])
      p0 <- rdc_loglik(s,M[[s]],W[[s]],Pr,Er[[s]],Ar0,P,E[[s]],A,fac)+log(1-pr[s,n])
      prob <- exp(p1-sumLog(c(p1,p0)))
      Ar[s,n] <- rbinom(1,1,prob)
    } else {
      scale <- 3
      prop <- rbeta(1,gamma,scale*gamma)
      Ar1 <- Ar
      Ar0 <- Ar
      Ar0[s,n] <- prop

      # Compute acceptance ratio
      p0 <- rdc_loglik(s,M[[s]],W[[s]],Pr,Er[[s]],Ar0,P,E[[s]],A,fac)
      p1 <- rdc_loglik(s,M[[s]],W[[s]],Pr,Er[[s]],Ar1,P,E[[s]],A,fac)
      ratio <- (exp((p0-p1)))
      if (runif(1)<ratio) {
        Ar[s,n] <- prop
      }
    }
  }
  return(Ar)
}

## Sample latent count
rdc_sample_zijs <- function(i,j,s,M_s,Pstar,Estar,Astar) {
  probs <- Pstar[i,]*Astar[s,]*t(Estar[,j])
  if (sum(is.na(probs))==0) {
    return(rmultinom(1,M_s[i,j],prob=probs))
  } else {
    return(rep(0,length(probs)))
  }
  return(rmultinom(1,M_s[i,j],prob=probs))
}

## Sample latent count (vectorized)
rdc_sample_zjs <- function(j,s,M_s,P,E_s) {
  probs <- sweep(P,2,E_s[,j],'*')
  return(sapply(1:nrow(probs),function(x) rmultinom(1,M_s[x,j],probs[x,])))
}

## Sample signature value
rdc_sample_pn <- function(n,W,E,A,Z,alpha_p,beta_p,S,G,gamma) {
  shape <- alpha_p[,n]+1+rowSums(sapply(1:S,function(s) gamma*A[s,n]*rowSums(Z[[s]][[n]])))
  rate <- beta_p[,n]+gamma*sum(sapply(1:S,function(s) sum(A[s,n]*E[[s]][n,]*W[[s]][1,1])))
  return(rgamma(nrow(alpha_p),shape,rate))
}

## Sample exposure value
rdc_sample_ens <- function(n,s,W_s,P,A,Z_s,alpha_e,beta_e,gamma) {
  shape <- alpha_e[[s]][n,]+gamma*A[s,n]*colSums(Z_s[[n]])+1
  rate <- beta_e[[s]][n,]+gamma*sum(P[,n]*A[s,n]*W_s[1,1])
  return(rgamma(ncol(alpha_e[[s]]),shape,rate))
}

## Sample signature hyperparameter beta
rdc_sample_betap <- function(P,alpha_p,a_p,b_p,gamma) {
  return(array(rgamma(nrow(alpha_p)*ncol(alpha_p),a_p+gamma*alpha_p+gamma,b_p+gamma*P),
               dim=dim(alpha_p)))
}

## Posterior density for signature hyperparameter alpha
rdc_alphainp_density <- function(i,n,P,y,beta_p,lambda_p,gamma) {
  return(log(lambda_p)-lambda_p*y+gamma*((y+1)*log(beta_p[i,n])+y*log(P[i,n])-lgamma(y+1)))
}

## Gamma proposal density
rdc_gamma_density <- function(y,x) {
  return(dgamma(y,x,1,log=TRUE))
}

## Sample signature hyperparameter alpha
rdc_sample_alphainp <- function(i,n,P,alpha_p,beta_p,lambda_p,gamma) {
  x <- alpha_p[i,n]
  y <- rgamma(1,x,1)
  alpha <- rdc_alphainp_density(i,n,P,y,beta_p,lambda_p,gamma)-
    rdc_alphainp_density(i,n,P,x,beta_p,lambda_p,gamma)+
    rdc_gamma_density(x,y)-rdc_gamma_density(y,x)
  if (runif(1) <= min(1,exp(alpha))) {
    return(y)
  } else {
    return(x)
  }
}

## Posterior density for exposure hyperparameter alpha
rdc_alphanjse_density <- function(n,s,E_s,y,beta_e,lambda_e,gamma) {
  return(log(lambda_e)-lambda_e*y+sum(gamma*((y+1)*log(beta_e[[s]][n,]))+
                                        (y*log(E_s[n,]))-lgamma(y+1)))
}

## Sample exposure hyperparameter alpha
rdc_sample_alphanse <- function(n,s,E_s,alpha_e,beta_e,lambda_e,gamma) {
  x <- alpha_e[[s]][n,1]
  y <- rgamma(1,x,1)
  alpha <- rdc_alphanjse_density(n,s,E_s,y,beta_e,lambda_e,gamma)-
    rdc_alphanjse_density(n,s,E_s,x,beta_e,lambda_e,gamma)+
    rdc_gamma_density(x,y)-rdc_gamma_density(y,x)
  if (runif(1) <= min(1,exp(alpha))) {
    return(y)
  } else {
    return(x)
  }
}

## Sample mu0
rdc_sample_mu0 <- function(s,n,G,E,A,alpha_e,beta_e,cov,delta,mu0_alpha,mu0_beta,gamma=1) {
  return(rgamma(1,mu0_alpha+sum(alpha_e[[s]][n,]+1),
                mu0_beta+sum(sapply(1:G[s],function(j)
                  E[[s]][n,j]*exp(-1*t(delta[[s]][,n,drop=F])%*%cov[[s]][,j,drop=F])))))
}

## Sample indicators
rdc_sample_indic <- function(delta,tau,c) {
  array(rbinom(length(delta),1,(delta^2/(tau^2*c^2))*dnorm(delta,0,c*tau)/
                 ((delta^2/(tau^2*c^2))*dnorm(delta,0,c*tau)+dnorm(delta,0,tau))),dim=dim(delta))
}

## Sample truncation points
rdc_sample_lambda <- function(delta,tau,c) {
  array(runif(length(delta),0,delta^2/(tau^2*c^2)),dim=dim(delta))
}

## Posterior for delta
rdc_delta_f <- function(x,s,l,n,a,tau,E,cov,mu0,alpha_e,A,delta) {
  delta[[s]][l,n] <- x
  return(sum(-(alpha_e[[s]][n,]+1)*t(delta[[s]][,n,drop=F])%*%cov[[s]]-
               (E[[s]][n,]*mu0[[s]][n,]/(exp(t(delta[[s]][,n,drop=F])%*%cov[[s]]))))-
           (1/(2*a[[s]][l,n]^2*tau^2))*x^2)
}

### Full Sampler
recovery_discovery_covariates <- function(M,cov,lambda_p=0.5,a_p=5,b_p=0.05,lambda_e=0.1,
                                          a_e=5,b_e=0.01,a_er=3,b_er=0.01,
                                          lambda_er=0.5,inits=NULL,fixed=FALSE,A.fixed=NULL,first=T) {
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
  lambda_er <- rep(lambda_er,S)
  a_er <- rep(a_er,S)
  b_er <- rep(b_er,S)
  data('signature_priors')
  sigs.r <- signature_priors

  # Initialize signature indicator matrices
  A <- array(0,dim=c(S,20))
  if (fixed) {
    A <- A.fixed[[1]]
  }
  N <- ncol(A)
  Ar <- array(0.25, dim=c(S,ncol(sigs.r[[1]])))
  colnames(Ar) <- colnames(sigs.r[[1]])
  if (fixed) {
    Ar <- A.fixed[[2]]
    to.keep <- which(colSums(Ar)>0)
    sigs.r[[1]] <- sigs.r[[1]][,to.keep]
    sigs.r[[2]] <- sigs.r[[2]][,to.keep]
    Ar <- Ar[,to.keep]
  }
  Nr <- ncol(Ar)

  # Initialize Bernoulli parameters
  pr <- array(0.25,dim=c(S,Nr))
  p <- array(0.1,dim=c(S,N))

  # Initialize hyperparameters
  alpha_p <- array(rexp(K*N,lambda_p),dim=c(K,N))
  beta_p <- array(rgamma(K*N,a_p,b_p),dim=c(K,N))
  alpha_e <- lapply(1:S,function(s) array(rexp(N*G[s],lambda_e[s]),dim=c(N,G[s])))
  alpha_pr <- as.matrix(sigs.r[[1]])
  beta_pr <- as.matrix(sigs.r[[2]])
  alpha_er <- lapply(1:S,function(s) array(rexp(Nr*G[s],lambda_er[s]),dim=c(Nr,G[s])))

  mu0 <- lapply(1:S,function(s) array(1,dim=c(N,G[s])))
  delta <- lapply(1:S,function(s) array(0,dim=c(nrow(cov[[1]]),N)))
  indic <- lapply(1:S,function(s) array(0,dim=c(nrow(cov[[1]]),N)))
  tau <- 0.3
  c <- 10
  an <- lapply(1:S,function(s) array(1,dim=c(nrow(cov[[1]]),N)))
  lambda <- lapply(1:S,function(s) array(runif(length(delta[[s]]),0,delta[[s]]^2/(tau*c)),
                                         dim=dim(delta[[s]])))
  beta_e <- lapply(1:S, function(s) mu0[[s]]*exp(-1*t(delta[[s]])%*%cov[[s]]))

  mu0r <- lapply(1:S,function(s) array(1,dim=c(Nr,G[s])))
  deltar <- lapply(1:S,function(s) array(0,dim=c(nrow(cov[[1]]),Nr)))
  indicr <- lapply(1:S,function(s) array(0,dim=c(nrow(cov[[1]]),Nr)))
  anr <- lapply(1:S,function(s) array(1,dim=c(nrow(cov[[1]]),Nr)))
  lambdar <- lapply(1:S,function(s) array(runif(length(deltar[[s]]),0,deltar[[s]]^2/(tau*c)),
                                          dim=dim(deltar[[s]])))
  beta_er <- lapply(1:S, function(s) mu0r[[s]]*exp(-1*t(deltar[[s]])%*%cov[[s]]))

  # Initialize signatures and exposures
  P <- array(rgamma(K*N,alpha_p+1,beta_p),dim=c(K,N))
  P <- P+min(P[P!=0])
  Pr <- array(rgamma(K*Nr,alpha_pr,beta_pr),dim=c(K,Nr))
  colnames(Pr) <- colnames(sigs.r[[1]])
  E <- lapply(1:S,function(s) array(rgamma(N*G[s],1,1),dim=c(N,G[s])))
  Er <- lapply(1:S,function(s) array(rgamma(Nr*G[s],1,1),dim=c(Nr,G[s])))
  if (!is.null(inits)) {
    P <- inits[[1]]
    E <- inits[[2]]
    Pr <- inits[[3]]
    Er <- inits[[4]]
    if (fixed) {
      Pr <- Pr[,to.keep]
      Er <- lapply(1:S,function(s) Er[[s]][to.keep,])
    }
  }

  # Initialize latent counts
  Z <- lapply(1:S,function(s) lapply(1:N,function(n) array(0,dim=c(K,G[s]))))
  Zr <- lapply(1:S,function(s) lapply(1:Nr,function(n) array(0,dim=c(K,G[s]))))
  Pstar <- cbind(P,Pr)
  Astar <- cbind(A,Ar)
  for (s in 1:S) {
    Estar <- rbind(E[[s]],Er[[s]])
    for (i in 1:K) {
      for (j in 1:G[s]) {
        z <- rdc_sample_zijs(i,j,s,M[[s]],Pstar,Estar,Astar)
        for (n in 1:N) {
          Z[[s]][[n]][i,j] <- z[n]
        }
        for (n in (N+1):length(z)) {
          Zr[[s]][[(n-N)]][i,j] <- z[n]
        }
      }
    }
  }

  # Initialize chains
  A_chain <- list()
  P_chain <- list()
  Ar_chain <- list()
  Pr_chain <- list()
  E_chain <- list()
  Er_chain <- list()
  indic_chain <- list()
  delta_chain <- list()
  indicr_chain <- list()
  deltar_chain <- list()

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

  st <- Sys.time()
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

    # Update recovery signatures
    if (params[l]==(-1)) {
      sigs.order <- sample(1:Nr,Nr,replace=FALSE)
      for (n1 in sigs.order) {
        Pr[,n1] <- rdc_sample_pn(n1,W,Er,Ar,Zr,alpha_pr-1,beta_pr,S,G,gamma)
      }
    } else {
      for (n1 in sample(1:Nr,Nr,replace=FALSE)) {
        maxit <- 1
        curr <- Inf
        while (curr > 0.9 & maxit < 10) {
          Pr[,n1] <- rdc_sample_pn(n1,W,Er,Ar,Zr,alpha_pr-1,beta_pr,S,G,temps[l])
          curr <- max(cosineDist(t(Pr[,n1,drop=F]),t(cbind(P,Pr[,-n1]))))
          maxit <- maxit+1
        }
      }
    }

    # Update recovery exposures
    sigs.order <- sample(1:Nr,Nr,replace=FALSE)
    for (s in 1:S) {
      for (n1 in sigs.order) {
        Er[[s]][n1,] <- rdc_sample_ens(n1,s,W[[s]],Pr,Ar,Zr[[s]],alpha_er,beta_er,gamma)
      }
    }

    # Update discovery signatures
    if (params[l]==(-1)) {
      sigs.order <- sample(1:N,N,replace=FALSE)
      for (n1 in sigs.order) {
        P[,n1] <- rdc_sample_pn(n1,W,E,A,Z,alpha_p,beta_p,S,G,gamma)
      }
    } else {
      for (n1 in sample(1:N,N,replace=FALSE)) {
        maxit <- 1
        curr <- Inf
        while (curr > 0.9 & maxit < 10) {
          P[,n1] <- rdc_sample_pn(n1,W,E,A,Z,alpha_p,beta_p,S,G,temps[l])
          curr <- max(cosineDist(t(P[,n1,drop=F]),t(cbind(Pr,P[,-n1]))))
          maxit <- maxit+1
        }
      }
    }

    # Update discovery exposures
    for (s in 1:S) {
      for (n1 in sample(1:N,N,replace=FALSE)) {
        E[[s]][n1,] <- rdc_sample_ens(n1,s,W[[s]],P,A,Z[[s]],alpha_e,beta_e,gamma)
      }
    }

    # Update Ar and A
    if (!fixed) {
      sigs.order <- sample(1:Nr,Nr,replace=FALSE)
      for (n in sigs.order) {
        Ar <- rdc_sample_Ar(n,Ar,M,W,Pr,Er,P,E,A,pr,S,fac,params[l])
      }

      sigs.order <- sample(1:N,N,replace=FALSE)
      for (n in sigs.order) {
        A <- rdc_sample_A(n,A,M,W,P,E,p,S,N,fac,Ar,Pr,Er,params[l])
      }
    }

    # Update latent counts
    if (params[l]==(-1)) {
      A.comb <- cbind(A,Ar)
      P.comb <- cbind(P,Pr)
      E.comb <- lapply(1:S,function(s) rbind(E[[s]],Er[[s]]))
      Z.comb <- lapply(1:S,function(s) lapply(1:(N+Nr),function(n)
        array(0,dim=c(nrow(M[[s]]),ncol(M[[s]])))))
      for (s in 1:S) {
        nonzero.sigs <- which(A.comb[s,]==1)
        for (j in 1:G[s]) {
          nonzero <- which(M[[s]][,j]!=0)
          if (length(nonzero)>0) {
            z <- rdc_sample_zjs(j,s,M[[s]][nonzero,,drop=F],P.comb[nonzero,nonzero.sigs,drop=F],
                            E.comb[[s]][nonzero.sigs,,drop=F])
            if (length(nonzero.sigs)==1) {
              Z.comb[[s]][[nonzero.sigs]][nonzero,j] <- z
            } else {
              for (n1 in 1:length(nonzero.sigs)) {
                Z.comb[[s]][[nonzero.sigs[n1]]][nonzero,j] <- z[n1,]
              }
            }
          }
        }
      }
      Zr <- lapply(1:S,function(s) lapply((N+1):(N+Nr),function(n) Z.comb[[s]][[n]]))
      Z <- lapply(1:S,function(s) lapply(1:N,function(n) Z.comb[[s]][[n]]))
    } else {
      Pstar <- cbind(P,Pr)
      Astar <- cbind(A,Ar)
      for (s in 1:S) {
        Estar <- rbind(E[[s]],Er[[s]])
        for (i in 1:K) {
          for (j in 1:G[s]) {
            z <- rdc_sample_zijs(i,j,s,M[[s]],Pstar,Estar,Astar)
            for (n in 1:N) {
              Z[[s]][[n]][i,j] <- z[n]
            }
            for (n in (N+1):length(z)) {
              Zr[[s]][[(n-N)]][i,j] <- z[n]
            }
          }
        }
      }
    }

    # Sample $\beta_{in}^p$
    beta_p <- rdc_sample_betap(P,alpha_p,a_p,b_p,gamma)

    # Sample $\alpha_{in}^p$
    for (i in 1:K) {
      for (n in 1:N) {
        alpha_p[i,n] <- rdc_sample_alphainp(i,n,P,alpha_p,beta_p,lambda_p,gamma)
      }
    }

    for (s in 1:S) {
      for (n in 1:N) {
        alpha_e[[s]][n,] <- rdc_sample_alphanse(n,s,E[[s]],alpha_e,beta_e,lambda_e[s],gamma)
      }
    }

    # Sample alpha^e_r
    for (s in 1:S) {
      for (n in 1:Nr) {
        alpha_er[[s]][n,] <- rdc_sample_alphanse(n,s,Er[[s]],alpha_er,beta_er,lambda_er[s],gamma)
      }
    }

    # Update covariate-related parameters in discovery component
    for (s in 1:S) {
      for (n in 1:N) {
        mu0[[s]][n,] <- rdc_sample_mu0(s,n,G,E,A,alpha_e,beta_e,cov,delta,mu0_alpha,
                                       mu0_beta,gamma=1)
      }
    }

    for (s in 1:S) {
      lambda[[s]] <- rdc_sample_lambda(delta[[s]],tau,c)
    }

    for (s in 1:S) {
      for (n in 1:N) {
        for (f in 1:nrow(cov[[1]])) {
          delta_f_inst <- function(x) { rdc_delta_f(x,s,f,n,an,tau,E,cov,mu0,alpha_e,A,delta) }
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
                  exp(-1*const+rdc_delta_f(val,s,f,n,an,tau,E,cov,mu0,alpha_e,A,delta))) }
                comp1 <- integrate(delta_f_inst2,thresh[2],10)$value
                comp2 <- integrate(delta_f_inst2,-10,thresh[1])$value
                pick <- rbinom(1,1,comp1/(comp1+comp2))
                if (pick==1) {
                  delta[[s]][f,n] <- as.numeric(armspp::arms(1,delta_f_inst,thresh[2],
                                                             10,metropolis=T))
                } else {
                  delta[[s]][f,n] <- as.numeric(armspp::arms(1,delta_f_inst,-10,
                                                             thresh[1],metropolis=T))
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
      indic[[s]] <- rdc_sample_indic(delta[[s]],tau,c)
      an[[s]] <- array(ifelse(indic[[s]]==0,1,c),dim=dim(indic[[s]]))
    }
    beta_e <- lapply(1:S, function(s) mu0[[s]]*exp(-1*t(delta[[s]])%*%cov[[s]]))

    # Update covariate-related parameters for the recovery component
    for (s in 1:S) {
      for (n in 1:Nr) {
        mu0r[[s]][n,] <- rdc_sample_mu0(s,n,G,Er,Ar,alpha_er,beta_er,cov,deltar,
                                    mu0_alpha,mu0_beta,gamma=1)
      }
    }

    for (s in 1:S) {
      lambdar[[s]] <- rdc_sample_lambda(deltar[[s]],tau,c)
    }

    for (s in 1:S) {
      for (n in 1:Nr) {
        for (f in 1:nrow(cov[[1]])) {
          delta_f_inst <- function(x) { rdc_delta_f(x,s,f,n,anr,tau,Er,cov,mu0r,alpha_er,Ar,deltar) }
          if (anr[[s]][f,n]==1) {
            deltar[[s]][f,n] <- armspp::arms(1,delta_f_inst,-10,10,metropolis=T)
          } else {
            thresh <- c(-1*sqrt(lambdar[[s]][f,n]*tau^2*c^2),sqrt(lambdar[[s]][f,n]*tau^2*c^2))
            deltar[[s]][f,n] <- armspp::arms(1,delta_f_inst,-10,10,metropolis=T)
            cnter <- 0
            while ((!(as.numeric(deltar[[s]][f,n])<thresh[1]|as.numeric(deltar[[s]][f,n])>thresh[2]))) {
              deltar[[s]][f,n] <- armspp::arms(1,delta_f_inst,-10,10,metropolis=T)
              cnter <- cnter+1
              if (cnter==10) {
                const <- max(c(delta_f_inst(thresh[1]),delta_f_inst(thresh[2])))
                delta_f_inst2 <- function(x) { sapply(x,function(val)
                  exp(-1*const+rdc_delta_f(val,s,f,n,anr,tau,Er,cov,mu0r,alpha_er,Ar,deltar))) }
                comp1 <- integrate(delta_f_inst2,thresh[2],10)$value
                comp2 <- integrate(delta_f_inst2,-10,thresh[1])$value
                pick <- rbinom(1,1,comp1/(comp1+comp2))
                if (pick==1) {
                  deltar[[s]][f,n] <- as.numeric(armspp::arms(1,delta_f_inst,thresh[2],
                                                              10,metropolis=T))
                } else {
                  deltar[[s]][f,n] <- as.numeric(armspp::arms(1,delta_f_inst,-10,thresh[1],
                                                              metropolis=T))
                }
              }
            }
          }
          deltar[[s]] <- array(as.numeric(deltar[[s]]),dim=c(nrow(cov[[1]]),Nr))
        }
      }
      deltar[[s]] <- array(as.numeric(deltar[[s]]),dim=c(nrow(cov[[1]]),Nr))
    }

    for (s in 1:S) {
      indicr[[s]] <- rdc_sample_indic(deltar[[s]],tau,c)
      anr[[s]] <- array(ifelse(indicr[[s]]==0,1,c),dim=dim(indicr[[s]]))
    }
    beta_er <- lapply(1:S, function(s) mu0r[[s]]*exp(-1*t(deltar[[s]])%*%cov[[s]]))

    A_chain[[l]] <- A
    P_chain[[l]] <- P
    Ar_chain[[l]] <- Ar
    Pr_chain[[l]] <- Pr
    E_chain[[l]] <- E
    Er_chain[[l]] <- Er
    indic_chain[[l]] <- indic
    indicr_chain[[l]] <- indicr
    delta_chain[[l]] <- delta
    deltar_chain[[l]] <- deltar

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
      paste0(paste0(sort(sapply(1:N,function(y) bitsToInt(A_chain[[x]][,y])),decreasing=T),
                    collapse=""),'_',
             paste0((sapply(1:Nr,function(y) bitsToInt(Ar_chain[[x]][,y]))),
                    collapse="")))
    As.tab <- sort(table(As),decreasing=TRUE)
    As.vec <- names(As.tab)[1:min(5,sum(As.tab>1))]
    inds <- sapply(As.vec,function(x) which(As==x)[1]+7700)
    Asr.list <- lapply(inds,function(x) Ar_chain[[x]])
    Psr.list <- lapply(inds,function(x) Pr_chain[[x]])
    Esr.list <- lapply(inds,function(x) lapply(1:S,function(s) Er_chain[[x]][[s]]))
    As.list <- lapply(inds,function(x) A_chain[[x]][,unique(c(1,2,which(colSums(A_chain[[x]])>0))),drop=F])
    Ps.list <- lapply(inds,function(x) P_chain[[x]][,unique(c(1,2,which(colSums(A_chain[[x]])>0))),drop=F])
    Es.list <- lapply(inds,function(x) lapply(1:S,function(s) E_chain[[x]][[s]][unique(c(1,2,which(colSums(A_chain[[x]])>0))),]))
    return(list(Asr.list,Psr.list,Esr.list,As.list,Ps.list,Es.list))
  } else {
    close(pb)
    liks <- sapply(1001:5000,function(x) sum(sapply(1:S,function(s)
      rdc_loglik(s,M[[s]],W[[s]],Pr_chain[[x]],Er_chain[[x]][[s]],
             Ar_chain[[x]],P_chain[[x]],E_chain[[x]][[s]],A_chain[[x]],fac))))
    bic <- (nrow(P)*sum(colSums(A)>0)+sum(rowSums(A)*G)+
              nrow(Pr)*ncol(Pr)+sum(rowSums(Ar)*G))*log(nrow(P)*sum(G))-2*max(liks)

    Pest <- array(0,dim=c(dim(P),4000))
    for (i in 1001:5000) {
      Pest[,,(i-1000)] <- P_chain[[i]]
    }
    Pest <- apply(Pest,c(1,2),median)

    Prest <- array(0,dim=c(dim(Pr),4000))
    for (i in 1001:5000) {
      Prest[,,(i-1000)] <- Pr_chain[[i]]
    }
    Prest <- apply(Prest,c(1,2),median)

    Eest <- list()
    for (s in 1:S) {
      if (sum(A[s,]==1)==0) {
        Eest[[s]] <- 0
      } else {
        Eest.s <- array(0,dim=c(dim(E[[s]][which(A[s,]==1),,drop=F]),4000))
        for (i in 1001:5000) {
          Eest.s[,,(i-1000)] <- E_chain[[i]][[s]][which(A[s,]==1),,drop=F]
        }
        Eest[[s]] <- apply(Eest.s,c(1,2),median)
      }
    }
    names(Eest) <- paste0('study',1:nrow(A))
    Erest <- list()
    for (s in 1:S) {
      Erest.s <- array(0,dim=c(dim(Er[[s]][which(Ar[s,]==1),,drop=F]),4000))
      for (i in 1001:5000) {
        Erest.s[,,(i-1000)] <- Er_chain[[i]][[s]][which(Ar[s,]==1),,drop=F]
      }
      Erest[[s]] <- apply(Erest.s,c(1,2),median)
    }
    names(Erest) <- paste0('study',1:nrow(A))

    indics <- lapply(1:S,function(s) {
      out <- t(sapply(1:nrow(cov[[1]]),function(f)
        rowMeans(sapply(1001:5000,function(x) indic_chain[[x]][[s]][f,]))))
      out <- out[,which(A[s,]==1),drop=F]
    })
    names(indics) <- paste0('study',1:nrow(A))
    indicsr <- lapply(1:S,function(s) {
      out <- t(sapply(1:nrow(cov[[1]]),function(f)
        rowMeans(sapply(1001:5000,function(x) indicr_chain[[x]][[s]][f,]))))
      out <- out[,which(Ar[s,]==1),drop=F]
    })
    names(indicsr) <- paste0('study',1:nrow(A))

    deltas <- lapply(1:S,function(s) {
      out <- t(sapply(1:nrow(cov[[1]]),function(f)
        rowMedians(sapply(1001:5000,function(x)
          ifelse(indic_chain[[x]][[s]][f,]==1,delta_chain[[x]][[s]][f,],NA)))))
      out <- out[,which(A[s,]==1),drop=F]
    })
    names(deltas) <- paste0('study',1:nrow(A))
    deltasr <- lapply(1:S,function(s) {
      out <- t(sapply(1:nrow(cov[[1]]),function(f)
        rowMedians(sapply(1001:5000,function(x)
          ifelse(indicr_chain[[x]][[s]][f,]==1,deltar_chain[[x]][[s]][f,],NA)))))
      out <- out[,which(Ar[s,]==1),drop=F]
    })
    names(deltasr) <- paste0('study',1:nrow(A))

    Pest <- Pest[,colSums(A)>0,drop=F]
    A <- A[,colSums(A)>0,drop=F]

    if (sum(A)>0) {
      colnames(A) <- paste0('d_sig',1:ncol(A))
      rownames(A) <- paste0('study',1:nrow(A))
      colnames(Pest) <- colnames(A)
    }

    colnames(Ar) <- paste0('r_sig',1:ncol(Ar))
    rownames(Ar) <- paste0('study',1:nrow(Ar))

    colnames(Prest) <- colnames(Ar)
    for (s in 1:S) {
      if (sum(A[s,]!=0)) {
        rownames(Eest[[s]]) <- colnames(A)[A[s,]==1]
        colnames(indics[[s]]) <- colnames(A)[A[s,]==1]
        colnames(deltas[[s]]) <- colnames(A)[A[s,]==1]
      }
      rownames(Erest[[s]]) <- colnames(Ar)[Ar[s,]==1]
      colnames(indicsr[[s]]) <- colnames(Ar)[Ar[s,]==1]
      colnames(deltasr[[s]]) <- colnames(Ar)[Ar[s,]==1]
    }

    return(list(BIC=bic,A=list(discovery=A,recovery=Ar),
                P=list(discovery=Pest,recovery=Prest),
                E=list(discovery=Eest,recovery=Erest),
                indics=list(discovery=indics,recovery=indicsr),
                coefs=list(discovery=deltas,recovery=deltasr),W=W))
  }
}

#' Run the Recovery-Discovery Sampler with Covariates
#'
#' This runs the recovery-discovery sampler for multi-study NMF using COSMIC v3.2 priors,
#' with covariates. Importantly, the order of mutational categories of the input data has
#' to match that of the COSMIC reference, see data('cosmic').
#'
#' @param M List of matrices corresponding to each study, in features-by-samples format. The order
#' of the mutational categories must match that of data('cosmic').
#' @param cov List of covariates matrices corresponding to each study,
#' in covariates-by-samples format; the same covariates in the same order should be supplied for
#' each study. Covariates must be numeric.
#' @param lambda_p Lambda_p, set to 0.5 by default.
#' @param a_p a_p, set to 5 by default.
#' @param b_p b_p, set to 0.05 by default.
#' @param lambda_e lambda_e, set to 0.1 by default.
#' @param a_e a_e, set to 5 by default.
#' @param b_e b_e, set to 0.01 by default.
#' @param a_er a_er, set to 3 by default.
#' @param b_er b_er, set to 0.01 by default.
#' @param lambda_er lambda_er, set to 0.5 by default.
#' @return A named list describing the best solution:
#' \itemize{
#'  \item{BIC: }{the BIC value of the solution.}
#'  \item{A: }{a list of the point estimates of the A matrices, indicating whether each signature is
#'  present in each study, corresponding to the discovery and recovery components respectively.}
#'  \item{P: }{a list of the posterior medians of the signatures matrices, corresponding to the
#'  discovery and recovery components respectively.}
#'  \item{E: }{a list of lists of the posterior medians of the exposures matrices for each study s,
#'  subsetted to just the signatures indicated by the A matrices to belong to that study,
#'  corresponding to the discovery and recovery components respectively.}
#'  \item{indics: }{a list of lists of the posterior inclusion probabilities for each study s,
#'  describing the probability that each covariate affects each signature indicated by A
#'  to belong to that study, corresponding to the discovery and recovery components respectively.}
#'  \item{coefs: }{a list of lists of the coefficients, conditional on inclusion, for each study s,
#'  describing the contribution of each covariate to each signature indicated by A to belong
#'  to that study, corresponding to the discovery and recovery components respectively.}
#'  \item{W: }{a list of the normalizing constants for each study; note that the input counts
#'  matrix M_s for each study s can be modeled as Poisson(W_s * P_s E_s), where P_s and E_s
#'  contain just the signatures that belong to study s.}
#'  }
#' @export
recovery_discovery_cov_fit <- function(M,cov,lambda_p=0.5,a_p=5,b_p=0.05,lambda_e=0.1,a_e=5,
                                       b_e=0.01,a_er=3,b_er=0.01,lambda_er=0.5) {

  message("Make sure the order of mutational categories in the input data matches
          the COSMIC reference, see data('cosmic')")

  initial <- recovery_discovery_covariates(M,cov,lambda_p,a_p,b_p,lambda_e,a_e,b_e,
                                  a_er,b_er,lambda_er,inits=NULL)
  outs <- lapply(1:length(initial[[1]]),function(x)
    recovery_discovery_covariates(M,cov,lambda_p,a_p,b_p,lambda_e,a_e,b_e,
                                  a_er,b_er,lambda_er,
                                  inits=list(initial[[5]][[x]],initial[[6]][[x]],
                                             initial[[2]][[x]],initial[[3]][[x]]),
                                  fixed=T,A.fixed=list(initial[[4]][[x]],initial[[1]][[x]]),first=F))
  pick <- which.min(sapply(outs,function(x) x[[1]]))
  return(outs[[pick]])
}
