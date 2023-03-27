library(SADISA)
library(armspp)

### Parameters
## M = S-list of K x G_s count matrices 
## W = S-list of K x G_s frequency matrices 
## P = K x N signature loadings
## E = N x G_s genome loadings 
## A = S x N signature indicator matrix
## Z = S-list of N-list of K x G_s latent variable matrices 
## alpha_p, beta_p = K x N hyperparameter matrices
## alpha_e, beta_e = S-lists of N x G_s hyperparameter matrices 
## alpha_IBP, beta_IBP = scalar hyperparameters for A
## lambda_p, a_p, b_p, lambda_e, a_e, b_e = scalar hyperprior parameters 

## Helper function to write binary matrix in condensed form
bitsToInt<-function(x) {
  packBits(rev(c(rep(FALSE, 32-length(x)%%32), as.logical(x))), "integer")
}

## Helper function to compute log of sum
sumLog <- function(x) {
  ord <- sort(x,decreasing=T)
  s <- ord[1]
  for (i in 2:length(x)) {
    s <- s + log(1+exp(ord[i]-s))
  }
  return(s)
}

## Helper function to find medians of each row
rowMedians <- function(x) {
  sapply(1:nrow(x),function(i) median(x[i,],na.rm=T))
}

## Compute log-likelihood of one study
loglik <- function(s,M_s,W_s,Pr,Er_s,Ar,P,E_s,A,fac) {
  term <- M_s*log(W_s*P%*%diag(A[s,])%*%E_s+W_s*Pr%*%diag(Ar[s,])%*%Er_s)
  term[M_s==0&is.na(term)] <- 0
  return(-1*sum(W_s*(P%*%diag(A[s,])%*%E_s+Pr%*%diag(Ar[s,])%*%Er_s)-term+
                  fac[M_s+1]))
}

## Sample one column of A
sample_A <- function(n,A,M,W,P,E,p,S,N,fac,Ar,Pr,Er,gamma) {
  for (s in 1:S) {
    Ar1 <- A
    Ar0 <- A
    Ar1[s,n] <- 1
    Ar0[s,n] <- 0
    p1 <- sum(sapply(1:S,function(s) gamma*loglik(s,M[[s]],W[[s]],Pr,Er[[s]],Ar,P,E[[s]],Ar1,fac)))+log(p[s,n])
    p0 <- sum(sapply(1:S,function(s) gamma*loglik(s,M[[s]],W[[s]],Pr,Er[[s]],Ar,P,E[[s]],Ar0,fac)))+log(1-p[s,n])
    prob <- exp(p1-tryCatch({ sumLog(c(p1,p0)) },error=function(cond) {
      return(max(p1,p0))}))
    A[s,n] <- rbinom(1,1,prob)
  }
  return(A)
}

## Sample p
sample_p <- function(p,A,a,b,S,N,gamma) {
  for (n in 1:N) {
    p[,n] <- rbeta(1,a+sum(A[,n]),b+sum(1-A[,n]))
  }
  return(p)
}

## Sample one column of Ar
sample_Ar <- function(n,Ar,M,W,Pr,Er,P,E,A,pr,S,fac,gamma) {
  for (s in 1:S) {
    Ar1 <- Ar
    Ar0 <- Ar
    Ar1[s,n] <- 1
    Ar0[s,n] <- 0
    p1 <- gamma*loglik(s,M[[s]],W[[s]],Pr,Er[[s]],Ar1,P,E[[s]],A,fac)+log(pr[s,n]) 
    p0 <- gamma*loglik(s,M[[s]],W[[s]],Pr,Er[[s]],Ar0,P,E[[s]],A,fac)+log(1-pr[s,n]) 
    prob <- exp(p1-tryCatch({ sumLog(c(p1,p0)) },error=function(cond) {
      return(max(p1,p0))}))
    Ar[s,n] <- rbinom(1,1,prob)
  }
  return(Ar)
}

## Sample pr
sample_pr <- function(pr,Ar,a,b,S,Nr,gamma) {
  for (n in 1:Nr) {
    pr[,n] <- rbeta(1,a+sum(Ar[,n]),b+sum(1-Ar[,n]))
  }
  return(pr)
}
 
## Sample latent count
sample_zijs <- function(i,j,s,M_s,P,Pr,E_s,Er_s,A,Ar) {
  probs <- cbind(P,Pr)[i,]*cbind(A,Ar)[s,]*t(rbind(E_s,Er_s)[,j])/(sum(P[i,]*A[s,]*t(E_s[,j]))+
                                                                     sum(Pr[i,]*Ar[s,]*t(Er_s[,j])))
  if (sum(is.na(probs))==0) {
    return(rmultinom(1,M_s[i,j],prob=probs))
  } else {
    return(rep(0,length(probs)))
  }
  return(rmultinom(1,M_s[i,j],prob=probs))
}

## Sample latent count (vectorized)
sample_zjs <- function(j,s,M_s,P,E_s) {
  probs <- sweep(P,2,E_s[,j],'*')
  return(sapply(1:nrow(probs),function(x) rmultinom(1,M_s[x,j],probs[x,])))
}

## Sample signature value
sample_pin <- function(i,n,W,E,A,Z,alpha_p,beta_p,S,G,gamma) {
  if (sum(A[,n])>0) {
    shape <- alpha_p[i,n]+gamma*sum(sapply(which(A[,n]==1),function(s) sum(Z[[s]][[n]][i,])))+1
  } else {
    shape <- alpha_p[i,n]+1
  }
  rate <- beta_p[i,n]+gamma*sum(sapply(1:S,function(s) sum(A[s,n]*E[[s]][n,]*W[[s]][i,])))
  return(rgamma(1,shape,rate))
}

## Sample exposure value
sample_enjs <- function(n,j,s,W_s,P,A,Z_s,alpha_e,beta_e,gamma) {
  shape <- alpha_e[[s]][n,j]+gamma*(A[s,n]==1)*sum(Z_s[[n]][,j])+1
  rate <- beta_e[[s]][n,j]+gamma*sum(P[,n]*A[s,n]*W_s[,j])
  return(rgamma(1,shape,rate))
}

## Sample signature hyperparameter beta
sample_betainp <- function(i,n,P,alpha_p,a_p,b_p,gamma) {
  return(rgamma(1,a_p+gamma*alpha_p[i,n]+gamma,b_p+gamma*P[i,n]))
}

## Sample exposure hyperparameter beta
sample_betanjse <- function(n,j,s,E_s,alpha_e,a_e,b_e,gamma) {
  return(rgamma(1,a_e+gamma*alpha_e[[s]][n,j]+gamma,b_e+gamma*E_s[n,j]))
}

## Posterior density for signature hyperparameter alpha
alphainp_density <- function(i,n,P,y,beta_p,lambda_p,gamma) {
  return(log(lambda_p)-lambda_p*y+gamma*((y+1)*log(beta_p[i,n])+y*log(P[i,n])-lgamma(y+1)))
}

## Gamma proposal density 
gamma_density <- function(y,x) {
  return(dgamma(y,x,1,log=TRUE))
}

## Sample signature hyperparameter alpha
sample_alphainp <- function(i,n,P,alpha_p,beta_p,lambda_p,gamma) {
  x <- alpha_p[i,n]
  y <- rgamma(1,x,1)
  alpha <- alphainp_density(i,n,P,y,beta_p,lambda_p,gamma)-alphainp_density(i,n,P,x,beta_p,lambda_p,gamma)+
    gamma_density(x,y)-gamma_density(y,x)
  if (runif(1) <= min(1,exp(alpha))) {
    return(y)
  } else {
    return(x)
  }
}

## Posterior density for exposure hyperparameter alpha
alphanjse_density <- function(n,s,E_s,y,beta_e,lambda_e,gamma) {
  return(log(lambda_e)-lambda_e*y+sum(gamma*((y+1)*log(beta_e[[s]][n,]))+(y*log(E_s[n,]))-lgamma(y+1)))
}

## Sample exposure hyperparameter alpha
sample_alphanse <- function(n,s,E_s,alpha_e,beta_e,lambda_e,gamma) {
  x <- alpha_e[[s]][n,1]
  y <- rgamma(1,x,1)
  alpha <- alphanjse_density(n,s,E_s,y,beta_e,lambda_e,gamma)-
    alphanjse_density(n,s,E_s,x,beta_e,lambda_e,gamma)+
    gamma_density(x,y)-gamma_density(y,x)
  if (runif(1) <= min(1,exp(alpha))) {
    return(y)
  } else {
    return(x)
  }
}

## Sample mu0
sample_mu0 <- function(s,G,E,A,alpha_e,beta_e,cov,delta,mu0_alpha,mu0_beta,gamma=1) {
  Ns <- sum(A[s,])
  if (Ns==0) {
    return(rgamma(1,mu0_alpha,mu0_beta))
  }
  return(rgamma(1,mu0_alpha+Ns*ncol(E[[s]])+sum(sapply(which(A[s,]==1),function(n) 
    sapply(1:G[s],function(j) alpha_e[[s]][n,j]))),
    mu0_beta+sum(sapply(which(A[s,]==1),function(n) 
      sapply(1:G[s],function(j) E[[s]][n,j]*(1/exp(t(cov[[s]][,j,drop=F])%*%delta[[s]][,n,drop=F])))))))
}

## Sample indicators
sample_indic <- function(delta,tau,c) {
  array(rbinom(length(delta),1,(delta^2/(tau^2*c^2))*dnorm(delta,0,c*tau)/
                 ((delta^2/(tau^2*c^2))*dnorm(delta,0,c*tau)+dnorm(delta,0,tau))),dim=dim(delta)) # check
}

## Sample truncation points
sample_lambda <- function(delta,tau,c) {
  array(runif(length(delta),0,delta^2/(tau^2*c^2)),dim=dim(delta)) # check
}

## Posterior for delta
delta_f <- function(x,s,l,n,a,tau,E,cov,mu0,alpha_e,A,delta) {
  delta[[s]][l,n] <- x
  return(A[s,n]*sum(-(alpha_e[[s]][n,]+1)*t(delta[[s]][,n,drop=F])%*%cov[[s]]-
                      (E[[s]][n,]*mu0[[s]][n,]/(exp(t(delta[[s]][,n,drop=F])%*%cov[[s]]))))-
           (1/(2*a[[s]][l,n]^2*tau^2))*x^2)
}

integrand <- function(x,k,alphap,w,betae,alphae,betap) {
  int <- -1*betap*x+(k+alphap-1)*log(x)-(alphae+k)*log(w*x+betae)+lchoose(k+alphae-1,k)+
    alphae*log(betae)+k*log(w)+alphap*log(betap)-lgamma(alphap)
  return(exp(int))
}

integrand2 <- function(x,k,alphap,w,betae,alphae,betap) {
  int <- -1*betap*x+(k+alphap-1)*log(x)-(alphae+k)*log(w*x+betae)+lchoose(k+alphae-1,k)+
    alphae*log(betae)+k*log(w)+alphap*log(betap)-lgamma(alphap)
  return(int)
}

### Full Sampler
recovery_discovery_covariates <- function(M,cov,lambda_p=0.5,a_p=5,b_p=0.05,lambda_e=1,a_e=5,b_e=0.01,a_er=5,b_er=0.01,
                                          lambda_er=1,a=1,b=0.05,inits=NULL,fixed=FALSE,A.fixed=NULL,first=T) {
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
  load('signature_priors.rda')
  sigs.r <- all
 
  # Initialize signature indicator matrices
  A <- array(0,dim=c(S,20))
  if (fixed) {
    A <- A.fixed[[1]]
  } else {
    A[,1:10] <- rbinom(S*10,1,0.5)
    A[rowSums(A)==0,1] <- 1
  }
  N <- ncol(A)
  Ar <- array(1, dim=c(S,ncol(sigs.r[[1]])))
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
  pr <- array(rbeta(Nr*S,1,1),dim=c(S,Nr))
  p <- array(rbeta(S*N,1,1),dim=c(S,N))
  
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
  lambda <- lapply(1:S,function(s) array(runif(length(delta[[s]]),0,delta[[s]]^2/(tau*c)),dim=dim(delta[[s]])))
  beta_e <- lapply(1:S, function(s) mu0[[s]]*exp(-1*t(delta[[s]])%*%cov[[s]]))
  
  mu0r <- lapply(1:S,function(s) array(1,dim=c(Nr,G[s])))
  deltar <- lapply(1:S,function(s) array(0,dim=c(nrow(cov[[1]]),Nr)))
  indicr <- lapply(1:S,function(s) array(0,dim=c(nrow(cov[[1]]),Nr)))
  anr <- lapply(1:S,function(s) array(1,dim=c(nrow(cov[[1]]),Nr)))
  lambdar <- lapply(1:S,function(s) array(runif(length(deltar[[s]]),0,deltar[[s]]^2/(tau*c)),dim=dim(deltar[[s]])))
  beta_er <- lapply(1:S, function(s) mu0r[[s]]*exp(-1*t(deltar[[s]])%*%cov[[s]]))
  
  # Initialize signatures and exposures
  P <- array(rgamma(K*N,alpha_p,beta_p),dim=c(K,N))
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
  for (s in 1:S) {
    for (i in 1:K) {
      for (j in 1:G[s]) {
        z <- sample_zijs(i,j,s,M[[s]],P,Pr,E[[s]],Er[[s]],A,Ar)
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
  llk_chain <- list()
  gamma_sched <- c(rep(0,10),sapply((-9):(-5),function(x) rep(10^x,100)),rep(10^(-4),500),
                   sapply(seq(0.1,8.9,by=0.1),function(x) rep((1+x)*10^(-4),100)),
                   sapply(seq(0,8.9,by=0.1),function(x) rep((1+x)*10^(-3),100)),
                   sapply(seq(0,8.9,by=0.1),function(x) rep((1+x)*10^(-2),100)),
                   sapply(seq(0,8.9,by=0.1),function(x) rep((1+x)*10^(-1),100)),rep(1,13090))
 
  if (first) {
    totaliters <- 50000
  } else {
    totaliters <- 1000
  }
 
  # Iterate
  for (l in 1:totaliters) {
    gamma <- 1
    
    # Update recovery signatures
    sigs.order <- sample(1:Nr,Nr,replace=FALSE)
    for (i in 1:K) {
      for (n1 in sigs.order) {
        Pr[i,n1] <- sample_pin(i,n1,W,Er,Ar,Zr,alpha_pr-1,beta_pr,S,G,gamma)
      }
    }
    
    # Update recovery exposures
    for (s in 1:S) {
      for (j in 1:G[s]) {
        for (n1 in sigs.order) {
          Er[[s]][n1,j] <- sample_enjs(n1,j,s,W[[s]],Pr,Ar,Zr[[s]],alpha_er,beta_er,gamma)
        }
      }
    }
    
    # Update discovery signatures
    sigs.order <- sample(1:N,N,replace=FALSE)
    for (i in 1:K) {
      for (n1 in sigs.order) {
        P[i,n1] <- sample_pin(i,n1,W,E,A,Z,alpha_p,beta_p,S,G,gamma)
      }
    }
    
    # Update discovery exposures
    for (s in 1:S) {
      for (j in 1:G[s]) {
        for (n1 in sigs.order) {
          E[[s]][n1,j] <- sample_enjs(n1,j,s,W[[s]],P,A,Z[[s]],alpha_e,beta_e,gamma)
        }
      }
    }
    
    # Update Ar and A
    if (!fixed) {
      sigs.order <- sample(1:Nr,Nr,replace=FALSE)
      for (n in sigs.order) {
        Ar <- sample_Ar(n,Ar,M,W,Pr,Er,P,E,A,pr,S,fac,gamma=gamma_sched[l])
      }
    
      sigs.order <- sample(1:N,N,replace=FALSE)
      for (n in sigs.order) {
        A <- sample_A(n,A,M,W,P,E,p,S,N,fac,Ar,Pr,Er,gamma=gamma_sched[l])
      }
    }

    # Update latent counts
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
          z <- sample_zjs(j,s,M[[s]][nonzero,,drop=F],P.comb[nonzero,nonzero.sigs,drop=F],
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
    
    # Sample pr and p
    pr <- sample_pr(pr,Ar,a=a,b=b,S,Nr,gamma)
    p <- sample_p(p,A,a=1,b=20,S,N,gamma)
    
    # Sample $\beta_{in}^p$
    for (i in 1:K) {
      for (n in 1:N) {
        beta_p[i,n] <- sample_betainp(i,n,P,alpha_p,a_p,b_p,gamma)
      }
    }
    
    # Sample $\alpha_{in}^p$
    for (i in 1:K) {
      for (n in 1:N) {
        alpha_p[i,n] <- sample_alphainp(i,n,P,alpha_p,beta_p,lambda_p,gamma)
      }
    }
    
    for (s in 1:S) {
      for (n in 1:N) {
        alpha_e[[s]][n,] <- sample_alphanse(n,s,E[[s]],alpha_e,beta_e,lambda_e[s],gamma)
      }
    }   
 
    # Sample alpha^e_r
    for (s in 1:S) {
      for (n in 1:Nr) {
        alpha_er[[s]][n,] <- sample_alphanse(n,s,Er[[s]],alpha_er,beta_er,lambda_er[s],gamma)
      }
    }
    
    # Update covariate-related parameters in discovery component
    for (s in 1:S) {
      mu0[[s]] <- array(sample_mu0(s,G,E,A,alpha_e,beta_e,cov,delta,mu0_alpha,mu0_beta,gamma=1),
                        dim=dim(mu0[[s]]))
    }
    
    for (s in 1:S) {
      lambda[[s]] <- sample_lambda(delta[[s]],tau,c)
    }
    
    for (s in 1:S) {
      for (n in 1:N) {
        for (f in 1:nrow(cov[[1]])) {
          delta_f_inst <- function(x) { delta_f(x,s,f,n,an,tau,E,cov,mu0,alpha_e,A,delta) }
          if (an[[s]][f,n]==1) {
            delta[[s]][f,n] <- arms(1,delta_f_inst,-5,5,metropolis=F)
          } else {
            thresh <- c(-1*sqrt(lambda[[s]][f,n]*tau^2*c^2),sqrt(lambda[[s]][f,n]*tau^2*c^2))
            delta[[s]][f,n] <- arms(1,delta_f_inst,-5,5,metropolis=F)
            cnter <- 0
            while ((!(as.numeric(delta[[s]][f,n])<thresh[1]|as.numeric(delta[[s]][f,n])>thresh[2]))) {
              delta[[s]][f,n] <- arms(1,delta_f_inst,-5,5,metropolis=F)
              cnter <- cnter+1
              if (cnter==10) { 
                delta_f_inst2 <- function(x) { sapply(x,function(val)  
                  exp(delta_f(val,s,f,n,an,tau,E,cov,mu0,alpha_e,A,delta))) }
                comp1 <- tryCatch({ integrate(delta_f_inst2,thresh[2],5)$value }, error=function(cond) {
                  NA
                })
                comp2 <- tryCatch({ integrate(delta_f_inst2,-5,thresh[1])$value },error=function(cond) {
                  NA
                })
                pick <- tryCatch({ rbinom(1,1,comp1/(comp1+comp2)) },error=function(cond) {
                  return(which.max(c(comp1,comp2))) })
                if (is.na(pick)) {
                  delta_f_inst3 <- function(x) { sapply(x,function(val)  
                    delta_f(val,s,f,n,an,tau,E,cov,mu0,alpha_e,A,delta)) }
                  comp1 <- integrate(delta_f_inst3,thresh[2],5)$value
                  comp2 <- integrate(delta_f_inst3,-5,thresh[1])$value
                  pick <- which.max(c(comp1,comp2))
                }
                if (pick==1) {
                  delta[[s]][f,n] <- as.numeric(arms(1,delta_f_inst,thresh[2],5,metropolis=F))
                } else {
                  delta[[s]][f,n] <- as.numeric(arms(1,delta_f_inst,-5,thresh[1],metropolis=F))
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
      indic[[s]] <- sample_indic(delta[[s]],tau,c)
      an[[s]] <- array(ifelse(indic[[s]]==0,1,c),dim=dim(indic[[s]]))
    }
    beta_e <- lapply(1:S, function(s) mu0[[s]]*exp(-1*t(delta[[s]])%*%cov[[s]]))
    
    # Update covariate-related parameters for the recovery component
    for (s in 1:S) {
      mu0r[[s]] <- array(sample_mu0(s,G,Er,Ar,alpha_er,beta_er,cov,deltar,mu0_alpha,mu0_beta,gamma=1),
                        dim=dim(mu0r[[s]])) 
    }
    
    for (s in 1:S) {
      lambdar[[s]] <- sample_lambda(deltar[[s]],tau,c)
    }
    
    for (s in 1:S) {
      for (n in 1:Nr) {
        for (f in 1:nrow(cov[[1]])) {
          delta_f_inst <- function(x) { delta_f(x,s,f,n,anr,tau,Er,cov,mu0r,alpha_er,Ar,deltar) }
          if (anr[[s]][f,n]==1) {
            deltar[[s]][f,n] <- arms(1,delta_f_inst,-5,5,metropolis=F)
          } else {
            thresh <- c(-1*sqrt(lambdar[[s]][f,n]*tau^2*c^2),sqrt(lambdar[[s]][f,n]*tau^2*c^2))
            deltar[[s]][f,n] <- arms(1,delta_f_inst,-5,5,metropolis=F)
            cnter <- 0
            while ((!(as.numeric(deltar[[s]][f,n])<thresh[1]|as.numeric(deltar[[s]][f,n])>thresh[2]))) {
              deltar[[s]][f,n] <- arms(1,delta_f_inst,-5,5,metropolis=F)
              cnter <- cnter+1
              if (cnter==10) { 
                delta_f_inst2 <- function(x) { sapply(x,function(val)  
                  exp(delta_f(val,s,f,n,anr,tau,Er,cov,mu0r,alpha_er,Ar,deltar))) }
                comp1 <- tryCatch({ integrate(delta_f_inst2,thresh[2],5)$value }, error=function(cond) {
                  NA
                })
                comp2 <- tryCatch({ integrate(delta_f_inst2,-5,thresh[1])$value },error=function(cond) {
                  NA
                })
                pick <- tryCatch({ rbinom(1,1,comp1/(comp1+comp2)) },error=function(cond) {
                  return(which.max(c(comp1,comp2))) })
                if (is.na(pick)) {
                  delta_f_inst3 <- function(x) { sapply(x,function(val)  
                    delta_f(val,s,f,n,anr,tau,Er,cov,mu0r,alpha_er,Ar,deltar)) }
                  comp1 <- integrate(delta_f_inst3,thresh[2],5)$value
                  comp2 <- integrate(delta_f_inst3,-5,thresh[1])$value
                  pick <- which.max(c(comp1,comp2))
                }
                if (pick==1) {
                  deltar[[s]][f,n] <- as.numeric(arms(1,delta_f_inst,thresh[2],5,metropolis=F))
                } else {
                  deltar[[s]][f,n] <- as.numeric(arms(1,delta_f_inst,-5,thresh[1],metropolis=F))
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
      indicr[[s]] <- sample_indic(deltar[[s]],tau,c)
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
    
    if (!first) {
      if (l%%10==0) {
        llk <- NULL
        for (s in 1:S) {
          for (n in which(A[s,]>0)) {
            for (i in 1:K) {
              for (j in 1:G[s]) {
                integrand3 <- function(x) { integrand2(x,k=Z[[s]][[n]][i,j],
                                                       alphap=alpha_p[i,n]+1,w=W[[s]][i,j],betae=beta_e[[s]][n,j],
                                                       alphae=alpha_e[[s]][n,j]+1,betap=beta_p[i,n]) }
                llk <- c(llk,integral_peak(integrand3))
              }
            }
          }
        }
        llk <- sum(as.numeric(llk))
        
        llk_chain[[l/10]] <- llk
        
        llk <- NULL
        for (s in 1:S) {
          for (n in which(Ar[s,]>0)) {
            for (i in 1:K) {
              for (j in 1:G[s]) {
                integrand3 <- function(x) { integrand2(x,k=Zr[[s]][[n]][i,j],
                                                       alphap=alpha_pr[i,n]+1,w=W[[s]][i,j],betae=beta_er[[s]][n,j],
                                                       alphae=alpha_er[[s]][n,j]+1,betap=beta_pr[i,n]) }
                llk <- c(llk,integral_peak(integrand3))
              }
            }
          }
        }
        llk <- sum(as.numeric(llk))
        
        llk_chain[[l/10]] <- llk_chain[[l/10]]+llk
      }
    }
  }
  
  if (first) {
    As <- sapply(1:length(Ar_chain),function(x) paste0(c(paste0((bitsToInt(Ar_chain[[x]][,colSums(Ar_chain[[x]])>0])),collapse=""),
                                                  paste0((bitsToInt(A_chain[[x]][,colSums(A_chain[[x]])>0])),collapse="")),collapse=""))
    As.tab <- sort(table(As),decreasing=TRUE)
    As.vec <- names(As.tab)[1:5]
    inds <- sapply(As.vec,function(x) which(As==x)[1])
    Asr.list <- lapply(inds,function(x) Ar_chain[[x]])
    Psr.list <- lapply(inds,function(x) Pr_chain[[x]])
    Esr.list <- lapply(inds,function(x) lapply(1:S,function(s) Er_chain[[x]][[s]]))
    As.list <- lapply(inds,function(x) A_chain[[x]][,unique(c(1,2,which(colSums(A_chain[[x]])>0))),drop=F])
    Ps.list <- lapply(inds,function(x) P_chain[[x]][,unique(c(1,2,which(colSums(A_chain[[x]])>0))),drop=F])
    Es.list <- lapply(inds,function(x) lapply(1:S,function(s) E_chain[[x]][[s]][unique(c(1,2,which(colSums(A_chain[[x]])>0))),]))
    return(list(As.list,Ps.list,Es.list,Asr.list,Psr.list,Esr.list))
  } else {
    marg.lik <- -1*sumLog(sapply(50:100,function(x) -1*llk_chain[[x]]))
    Pest <- array(0,dim=c(dim(P),200))
    for (i in 801:1000) {
      Pest[,,(i-800)] <- P_chain[[i]]
    }
    Pest <- apply(Pest,c(1,2),median)

    Prest <- array(0,dim=c(dim(Pr),200))
    for (i in 801:1000) {
      Prest[,,(i-800)] <- Pr_chain[[i]]
    }
    Prest <- apply(Prest,c(1,2),median)
    
    indics <- lapply(1:S,function(s) {
      out <- t(sapply(1:nrow(cov[[1]]),function(f)
        rowMeans(sapply(801:1000,function(x) indic_chain[[x]][[s]][f,]))))
      out[,which(A[s,]==0)] <- NA
      out
    })
    indicsr <- lapply(1:S,function(s) {
      out <- t(sapply(1:nrow(cov[[1]]),function(f)
        rowMeans(sapply(801:1000,function(x) indicr_chain[[x]][[s]][f,]))))
      out[,which(Ar[s,]==0)] <- NA
      out
    })
    
    deltas <- lapply(1:S,function(s) {
      out <- t(sapply(1:nrow(cov[[1]]),function(f)
        rowMedians(sapply(801:1000,function(x) 
          ifelse(indic_chain[[x]][[s]][f,]==1,delta_chain[[x]][[s]][f,],NA)))))
      out[,which(A[s,]==0)] <- NA
      out
    })
    deltasr <- lapply(1:S,function(s) {
      out <- t(sapply(1:nrow(cov[[1]]),function(f)
        rowMedians(sapply(801:1000,function(x) 
          ifelse(indicr_chain[[x]][[s]][f,]==1,deltar_chain[[x]][[s]][f,],NA)))))
      out[,which(Ar[s,]==0)] <- NA
      out
    })
    
    return(list(marg.lik,list(A,Ar),list(Pest,Prest),list(indics,indicsr),list(deltas,deltasr)))
  }  
}

cosineDist <- function(x,y) {
  (x%*%t(y))/(sqrt(rowSums(x^2) %*% t(rowSums(y^2))))
}

label_A <- function(A,P) {
  true.sigs <- read.table('COSMIC_v3.2_SBS_GRCh38.txt',header=T)
  P <- P[,which(colSums(A)>0),drop=F]
  if (ncol(P)==0) {
    return(0)
  }
  cosD <- cosineDist(t(P),t(true.sigs[,2:ncol(true.sigs)]))
  unl <- 1
  sigs.new <- array(P,dim=c(96,ncol(P)))
  colnames(sigs.new) <- 1:ncol(sigs.new)
  for (i in 1:ncol(sigs.new)) {
    m <- which.max(cosD[i,])
    if (cosD[i,m]>=0.95) {
      colnames(sigs.new)[i] <- colnames(true.sigs[,2:ncol(true.sigs)])[m]
    } else if (cosD[i,m]>=0.7) {
      colnames(sigs.new)[i] <- paste('Pseudo',colnames(true.sigs[,2:ncol(true.sigs)][m]))
    } else {
      colnames(sigs.new)[i] <- paste('Unlabeled',unl)
      unl <- unl + 1
    }
  }
  colnames(sigs.new) <- make.unique(colnames(sigs.new))
  colnames(A) <- 1:ncol(A)
  colnames(A)[which(colSums(A)>0)] <- colnames(sigs.new)
  return(A[,which(colSums(A)>0),drop=F])
}
