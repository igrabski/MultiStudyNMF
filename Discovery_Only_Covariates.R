library(armspp)
library(SADISA)

### Parameters
## M = S-list of K x G_s count matrices 
## W = S-list of K x G_s frequency matrices 
## P = K x N signature loadings
## E = N x G_s genome loadings 
## A = S x N signature indicator matrix
## Z = S-list of N-list of K x G_s latent variable matrices 
## alpha_p, beta_p = K x N hyperparameter matrices
## alpha_e, beta_e = S-lists of N x G_s hyperparameter matrices 
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

## Computes log-likelihood for a single study
loglik <- function(E.new,W_s,P,M_s,fac) {
  term <- M_s*log(W_s*P%*%E.new)
  term[M_s==0&is.na(term)] <- 0
  return(-1*sum(W_s*(P%*%E.new)-term+
                  fac[M_s+1]))
}

## Computes log-likelihood across all studies
loglik.all <- function(P,A,E,W,M,fac,S) {
  llk <- sum(sapply(1:S,function(s) loglik(E[[s]][which(A[s,]==1),,drop=F],W[[s]],P[,which(A[s,]==1),drop=F],M[[s]],fac)))
  return(llk)
}

## Sample one column of A
sample_A_marg <- function(n,A,M,W,P,E,p,S,fac,gamma) {
  # For each study... 
  for (s in 1:S) {
    
    # Create a version of A where study contains the signature, and where it doesn't
    Ar1 <- A
    Ar0 <- A
    Ar1[s,n] <- 1
    Ar0[s,n] <- 0
    
    # Compute likelihoods under both scenarios 
    p1 <- gamma*loglik.all(P,Ar1,E,W,M,fac,S)+log(p[s,n])
    p0 <- gamma*loglik.all(P,Ar0,E,W,M,fac,S)+log(1-p[s,n])
    prob <- exp(p1-tryCatch({ sumLog(c(p1,p0)) },error=function(cond) {
      return(max(p1,p0))}))
    
    # Draw according to these probabilities
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

## Sample latent count
sample_zijs <- function(i,j,s,M_s,P,E_s,A,N) {
  probs <- P[i,]*A[s,]*t(E_s[,j])/sum(P[i,]*A[s,]*t(E_s[,j]))
  if (sum(is.na(probs))==0) {
    return(rmultinom(1,M_s[i,j],prob=probs))
  } else {
    return(rep(0,length(probs)))
  }
}

## Sample latent count (vectorized)
sample_zjs <- function(j,s,M_s,P,E_s,N) {
  probs <- sweep(P,2,E_s[,j],'*')
  return(sapply(1:nrow(probs),function(x) rmultinom(1,M_s[x,j],probs[x,])))
}

## Sample signature value (vectorized)
sample_pi <- function(i,W,E,A,Z,alpha_p,beta_p,S,G,gamma) {
  shape <- alpha_p[i,]+1
  nonzero <- which(colSums(A)>0)
  shape[nonzero] <- shape[nonzero]+sapply(nonzero,function(n) sum(sapply(which(A[,n]==1),function(s) 
    gamma*sum(Z[[s]][[n]][i,]))))
  rate <- beta_p[i,]+sapply(1:ncol(A),function(n) sum(sapply(1:S,function(s) 
    gamma*sum(A[s,n]*E[[s]][n,]*W[[s]][i,]))))
  return(rgamma(ncol(A),shape,rate))
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

### Full Sampler
discovery_covariates <- function(M,cov,lambda_p=0.5,a_p=5,b_p=0.05,lambda_e=1,a_e=5,b_e=0.01,
                                 inits=NULL,fixed=F,A.fixed=NULL,first=T) {
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
  A <- array(0,dim=c(S,50))
  if (!fixed) {
    A[,1:50] <- rbinom(S*50,1,0.5)
    A[rowSums(A)==0,1] <- 1
  } else {
    A <- A.fixed
  }
  N <- ncol(A)
  
  # Initialize Bernoulli parameters
  p <- array(rbeta(S*N,1,1),dim=c(S,N))
  
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
  lambda <- lapply(1:S,function(s) array(runif(length(delta[[s]]),0,delta[[s]]^2/(tau*c)),dim=dim(delta[[s]])))
  beta_e <- lapply(1:S, function(s) mu0[[s]]*exp(-1*t(delta[[s]])%*%cov[[s]]))
  
  # Initialize signatures and exposures
  P <- array(rgamma(K*N,alpha_p,beta_p),dim=c(K,N))
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
        z <- sample_zijs(i,j,s,M[[s]],P,E[[s]],A,N)
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
  llk_chain <- list()
  gamma_sched <- c(rep(0,10),sapply((-9):(-5),function(x) rep(10^x,100)),rep(10^(-4),500),
                   sapply(seq(0.1,8.9,by=0.1),function(x) rep((1+x)*10^(-4),100)),
                   sapply(seq(0,8.9,by=0.1),function(x) rep((1+x)*10^(-3),100)),
                   sapply(seq(0,8.9,by=0.1),function(x) rep((1+x)*10^(-2),100)),
                   sapply(seq(0,8.9,by=0.1),function(x) rep((1+x)*10^(-1),100)),rep(1,13090))

  # Iterate
  if (first) {
    totaliters <- 50000
  } else {
    totaliters <- 1000
  }

  for (l in 1:totaliters) {
    gamma <- 1
    
    # Update exposures
    for (s in 1:S) {
      for (j in 1:G[s]) {
        for (n1 in sample(1:N,N,replace=FALSE)) {
          E[[s]][n1,j] <- sample_enjs(n1,j,s,W[[s]],P,A,Z[[s]],alpha_e,beta_e,gamma)
        }
      }
    }
    
    # Update signatures
    for (i in 1:K) {
      P[i,] <- sample_pi(i,W,E,A,Z,alpha_p,beta_p,S,G,gamma)
    }
    
    # Update columns of A
    if (!fixed) { 
      sigs.order <- sample(1:N,N,replace=FALSE)
      for (n in sigs.order) {
        A <- sample_A_marg(n,A,M,W,P,E,p,S,fac,gamma=gamma_sched[l])
      }
    }

    # Update latent counts 
    for (s in 1:S) {
      nonzero.sigs <- which(A[s,]==1)
      for (j in 1:G[s]) {
        nonzero <- which(M[[s]][,j]!=0)
        if (length(nonzero)>0) {
          z <- sample_zjs(j,s,M[[s]][nonzero,,drop=F],P[nonzero,nonzero.sigs,drop=F],E[[s]][nonzero.sigs,,drop=F],N)
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
    
    # Update exposure hyperparameters alpha
    for (s in 1:S) {
      for (n in 1:N) {
        alpha_e[[s]][n,] <- sample_alphanse(n,s,E[[s]],alpha_e,beta_e,lambda_e[s],gamma)
      }
    }
	  
    # Update Bernoulli parameters
    p <- sample_p(p,A,a=0.8,b=0.8,S,N,gamma)
    
    # Update signature hyperparameters beta
    for (i in 1:K) {
      for (n in 1:N) {
        beta_p[i,n] <- sample_betainp(i,n,P,alpha_p,a_p,b_p,gamma)
      }
    }
    
    # Update signature hyperparameters alpha
    for (i in 1:K) {
      for (n in 1:N) {
        alpha_p[i,n] <- sample_alphainp(i,n,P,alpha_p,beta_p,lambda_p,gamma)
      }
    }

    # Update covariate-related parameters
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

    # Store in chains
    A_chain[[l]] <- A
    P_chain[[l]] <- P
    E_chain[[l]] <- E
    indic_chain[[l]] <- indic
    delta_chain[[l]] <- delta

    if (!first) {
      if (l%%10==0) {
        llk <- NULL
        for (s in 1:S) {
          for (n in which(A[s,]>0)) {
            for (i in 1:K) {
              for (j in 1:G[s]) {
                t <- tryCatch(expr={log(integrate(integrand,lower=0,upper=Inf,k=Z[[s]][[n]][i,j],
                                              alphap=alpha_p[i,n]+1,w=W[[s]][i,j],betae=beta_e[[s]][n,j],
                                              alphae=alpha_e[[s]][n,j]+1,betap=beta_p[i,n])$value)},
                          error=function(e){log(integrate(integrand,lower=0,upper=1000,k=Z[[s]][[n]][i,j],
                                                          alphap=alpha_p[i,n]+1,w=W[[s]][i,j],betae=beta_e[[s]][n,j],
                                                          alphae=alpha_e[[s]][n,j]+1,betap=beta_p[i,n])$value)})
                llk<-c(llk,t)
              }
            }
          }
        }
        llk <- sum(as.numeric(llk))

        llk_chain[[l/10]] <- llk
     }
   }
  }

  if (first) {
    As <- sapply(1:length(A_chain),function(x) paste0((bitsToInt(A_chain[[x]][,colSums(A_chain[[x]])>0])),collapse=""))
    As.tab <- sort(table(As),decreasing=TRUE)
    As.vec <- names(As.tab)[1:5]
    inds <- sapply(As.vec,function(x) which(As==x)[1])
    As.list <- lapply(inds,function(x) A_chain[[x]][,colSums(A_chain[[x]])>0])
    Ps.list <- lapply(inds,function(x) P_chain[[x]][,colSums(A_chain[[x]])>0])
    Es.list <- lapply(inds,function(x) lapply(1:S,function(s) E_chain[[x]][[s]][colSums(A_chain[[x]])>0,]))
    return(list(As.list,Ps.list,Es.list))
  } else {
    marg.lik <- -1*sumLog(sapply(50:100,function(x) -1*llk_chain[[x]]))
    Pest <- array(0,dim=c(dim(P),200))
    for (i in 801:1000) {
      Pest[,,(i-800)] <- P_chain[[i]]
    }
    Pest <- apply(Pest,c(1,2),median)
    
    indics <- lapply(1:S,function(s) {
      out <- t(sapply(1:nrow(cov[[1]]),function(f)
        rowMeans(sapply(801:1000,function(x) indic_chain[[x]][[s]][f,]))))
      out[,which(A[s,]==0)] <- NA
      out
      })
    
    deltas <- lapply(1:S,function(s) {
      out <- t(sapply(1:nrow(cov[[1]]),function(f)
        rowMedians(sapply(801:1000,function(x) 
          ifelse(indic_chain[[x]][[s]][f,]==1,delta_chain[[x]][[s]][f,],NA)))))
      out[,which(A[s,]==0)] <- NA
      out
    })
    
    return(list(marg.lik,A.fixed,Pest,indics,deltas))
  }  
}

cosineDist <- function(x,y) {
  (x%*%t(y))/(sqrt(rowSums(x^2) %*% t(rowSums(y^2))))
}