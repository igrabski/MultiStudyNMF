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

## Helper function to compute cosine similarity
cosineDist <- function(x,y) {
  (x%*%t(y))/(sqrt(rowSums(x^2) %*% t(rowSums(y^2))))
}

#' Label the A matrix by COSMIC Signatures
#'
#' This labels each column of an estimated A matrix by the closest COSMIC v3.2 signature,
#' as long as the cosine similarity is above 0.7. Note that if the cosine similarity is
#' below 0.95, "Pseudo" will be appended to the signature name. Duplicates will be indicated
#' by suffixes.
#'
#' @param A Estimate of the A matrix.
#' @param P Corresponding estimate of the signatures matrix.
#' @return The A matrix with columns labeled as indicated above.
#' @export
label_A <- function(A,P) {
  data('cosmic')
  P <- P[,which(colSums(A)>0),drop=F]
  if (ncol(P)==0) {
    return(0)
  }
  cosD <- cosineDist(t(P),t(cosmic[,2:ncol(cosmic)]))
  unl <- 1
  sigs.new <- array(P,dim=c(96,ncol(P)))
  colnames(sigs.new) <- 1:ncol(sigs.new)
  for (i in 1:ncol(sigs.new)) {
    m <- which.max(cosD[i,])
    if (cosD[i,m]>=0.95) {
      colnames(sigs.new)[i] <- colnames(cosmic[,2:ncol(cosmic)])[m]
    } else if (cosD[i,m]>=0.7) {
      colnames(sigs.new)[i] <- paste('Pseudo',colnames(cosmic[,2:ncol(cosmic)][m]))
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

## Helper function to find medians of each row
rowMedians <- function(x) {
  sapply(1:nrow(x),function(i) median(x[i,],na.rm=T))
}

#' Annotate signatures by proportion with high exposures.
#'
#' This populates each entry of the A matrix by the proportion of samples in that study
#' with relative exposures above a specified threshold.
#'
#' @param out Output of sampler.
#' @param thresh Threshold for relative exposures; default 0.05.
#' @param label Whether the A matrix should be annotated by its closest match to
#' COSMIC v3.2; default TRUE.
#' @return The A matrix, where each entry indicates the proportion of samples in that study
#' with relative exposures above `thresh`.
#' @export
proportion_A <- function(out,thresh=0.05,label=T) {
  if (label) {
    if (is.list(out$A)) {
      if (length(out$A[[1]])>1) {
        A <- label_A(cbind(out$A[[1]],out$A[[2]]),cbind(out$P[[1]],out$P[[2]]))
        P <- cbind(out$P[[1]],out$P[[2]])
      } else {
        A <- label_A(out$A[[2]],out$P[[2]])
        P <- out$P[[2]]
      }
    } else {
      A <- label_A(out$A,out$P)
    }
  } else {
    if (is.list(out$A)) {
      if (length(out$A[[1]])>1) {
        A <- cbind(out$A[[1]],out$A[[2]])
        P <- cbind(out$P[[1]],out$P[[2]])
      } else {
        A <- out$A[[2]]
        P <- out$P[[2]]
      }
    } else {
      A <- out$A
    }
  }

  for (s in 1:nrow(A)) {
    if (is.list(out$A)) {
      if (length(out$E[[1]][[s]])>1) {
        E <- rbind(out$E[[1]][[s]],out$E[[2]][[s]])
        E <- sweep(E,1,colSums(P[,A[s,]==1]),'*')
        E <- sweep(E,2,colSums(E),'/')
        A[s,which(A[s,]==1)] <- rowSums(E>thresh)/ncol(E)
      } else {
        E <- sweep(out$E[[2]][[s]],1,colSums(P[,A[s,]==1]),'*')
        E <- sweep(E,2,colSums(E),'/')
        A[s,which(A[s,]==1)] <- rowSums(E>thresh)/ncol(E)
      }
    } else {
      E <- sweep(out$E[[s]],1,colSums(out$P[,A[s,]==1]),'*')
      E <- sweep(E,2,colSums(E),'/')
      A[s,which(A[s,]==1)] <- rowSums(E>thresh)/ncol(E)
    }
  }
  return(A)
}
