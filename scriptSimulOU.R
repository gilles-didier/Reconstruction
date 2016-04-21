#
#    A comparison of ancestral state reconstruction methods for quantitative characters
#    see http://biorxiv.org/content/early/2016/01/25/037812
#    Copyright (C) 2016  Manuela ROYER-CARENZI and Gilles DIDIER
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

###############################################
##############################################
##  OU Simulations
###############################################
###############################################



##################################################
# Simulation Parameters
###################################################


nomFichiertxt <- paste("TreeWeb.txt")

# Simulation number
N <- 1


# root value 
# "constant" then z0=m0
# or "gaussian" then Z0 with law N(m0, v0)
m0 <- 100
v0 <- 30
root <- "constant"


# sigma 
sigmaM <- c(3,10)
nSM <- length(sigmaM)


# alpha
aMax <- 0.5
naM <- 20
aM <- seq(0, aMax, length=naM)
#


# theta
bM <- c(m0,200)
nbM <- length(bM)
#





##############################################
##  Phylogenetic tree format
###############################################


library(ade4)
library(ape)
library(nlme)
library(ouch)

#
#
#ecrit le noeud d'index ind de l'arbre sur le flux out
ecritnoeud <- function(out, arbre, ind) {
    cour <- arbre$edge[ind,2]
    l <- which(arbre$edge[, 1] == cour)
    if(length(l)>0){
      cat("(", file = out)
      ecritnoeud(out, arbre, l[1])
      if(length(l)>1){
	for(i in 2:length(l)) {
	  cat(", ", file = out)
	  ecritnoeud(out, arbre, l[i])
	}
      }
      cat(")",cour,":",arbre$length[ind], file = out)
    } else {
      cat(cour,":",arbre$edge.length[ind], file = out)
    }
}

#ecrit l'arbre sur le flux out
ecritarbre <- function(out, arbre) {
    root<-racine(arbre)
    l <- which(arbre$edge[, 1] == root)
    cat("(", file = out)
    if(length(l)>0){
      ecritnoeud(out, arbre, l[1])
      if(length(l)>1){
	for(i in 2:length(l)) {
	  cat(", ", file = out)
	  ecritnoeud(out, arbre, l[i])
	}
      }
    }
    cat(")", root, ";", file = out)
}

#renvoie la racine d'arbre
racine <- function(arbre) {
    m <- max(arbre[["edge"]])
    n<-dim(arbre$edge)[1]
    flag <-logical(m)
    for(i in 1:n) {
      flag[arbre$edge[i,2]] <- TRUE
    }
    for(i in 1:m) {
      if(!flag[i]) {
	root<-i
      }
    }
    root
}

nombreInternes <- function(arbre) {
    m <- max(arbre[["edge"]])
    n<-dim(arbre$edge)[1]
    flag <-logical(m)
    for(i in 1:n) {
      flag[arbre$edge[i,1]] <- TRUE
    }
    interne <- 0
    for(i in 1:m) {
      if(flag[i]) {
	interne<-interne+1
      }
    }
    interne
}

# renvoie une liste dont la première composante est l'arbre renumerote
# de telle sorte que l'index d'un enfant est superieur à  celui de son père,
# la seconde compopsante est la fonction de l'index initial vers le second,
# et la troisième son inverse 
# (attention problème pour l'image de 0 mise à  l'image du max)
#
renumeroteArbre <- function(arbre) {
  m <- max(arbre[["edge"]])
  v<-numeric(m)
  t<-numeric(m)
  stack<-numeric(m)
  istack<-1
  stack[istack]<-racine(arbre)
  codeI<-1
  codeL<-nombreInternes(arbre)+1
  while(istack>0){
    cour<-stack[istack]
    istack<-istack-1
    l <- which(arbre$edge[, 1] == cour)
    if(length(l)>0){
      v[cour] <- codeI
      t[codeI] <- cour
      codeI <- codeI+1
      for(i in 1:length(l)) {
	istack<-istack+1
	stack[istack] <- arbre$edge[l[i], 2]
      }
    } else {
      v[cour] <- codeL
      t[codeL] <- cour
      codeL <- codeL+1
    }
  }
  arbrebis<-arbre
  for(i in 1:dim(arbre$edge)[1]) {
    arbrebis$edge[i,1] <- v[arbre$edge[i,1]]
    arbrebis$edge[i,2] <- v[arbre$edge[i,2]]
  }
  l <- list(arbre = arbrebis, cod = v, dec = t)
  l
}

#renvoie un vecteur contenant les temps de la racine vers les noeuds
calculeV_times <- function(arbre) {
  m <- max(arbre[["edge"]])
  V <- numeric(m)
  V[1] <- 0
  for(i in 2:(m)) {
    l <- which(arbre$edge[, 2] == i)
    if(length(l)>0){
      V[i] <- V[arbre$edge[l[1], 1]]+arbre$edge.length[l[1]]
    }
  }
  V
}

#renvoie un vecteur contenant des 1
calculeV_ones <- function(arbre, a) {
  m <- max(arbre[["edge"]])
  V <- numeric(m)
  V[1] <- 1
  for(i in 1:(m)) {
    V[i]<- 1
  }
  V
}

#renvoie un vecteur contenant les exponentielles des temps de la racine vers les noeuds
calculeV_exptimes <- function(arbre, a) {
  m <- max(arbre[["edge"]])
  V <- numeric(m)
  V[1] <- 1
  for(i in 2:(m)) {
    l <- which(arbre$edge[, 2] == i)
    if(length(l)>0){
      V[i] <- V[arbre$edge[l[1], 1]]*exp(-a*arbre$edge.length[l[1]])
    }
  }
  V
}

#renvoie des temps entre chaque noeud
calculeM_times <- function(arbre) {
  m <- max(arbre[["edge"]])
  M <- matrix(Inf,nrow=m,ncol=m)
  for(i in 1:m) {
    M[i,i] = 0
    l <- which(arbre$edge[, 2] == i)
    if(length(l)>0){
      for(j in 1:m) {
	if(j == arbre$edge[l[1], 1] || !is.infinite(M[j,arbre$edge[l[1], 1]])) {
	  M[j,i] <- M[j, arbre$edge[l[1], 1]]+arbre$edge.length[l[1]]
	}
      }
    }
  }
  M
}

#calcule la matrice C selon le modèle ABM
#
calculeC_ABM <- function(arbre) {
  m <- max(arbre[["edge"]])
  C <- matrix(0,nrow=m,ncol=m)
  for(i in 1:(m)) {
    l <- which(arbre$edge[, 2] == i)
    if(length(l)>0){
      for(j in 1:(m)) {
	C[j,i] <- C[j, arbre$edge[l[1], 1]]
      }
    }
    C[i,i]<-1;
  }
  t(C)
}

#calcule la matrice C selon le modèle OU
#
calculeC_OU <- function(arbre, a) {
  m <- max(arbre[["edge"]])
  C <- matrix(0,nrow=m,ncol=m)
  for(i in 1:(m)) {
    l <- which(arbre$edge[, 2] == i)
    if(length(l)>0){
      for(j in 1:(m)) {
	  C[j,i] <- C[j, arbre$edge[l[1], 1]]*exp(-a*arbre$edge.length[l[1]])
      }
    }
    C[i,i]<-1;
  }
  t(C)
}

#calcule la matrice C selon le modèle type qui vaut ABM ou OU
calculeC <- function(type, arbre, a) {
  switch(type, ABM = calculeC_ABM(arbre), OU = calculeC_OU(arbre, a))
}





#########################
#calcul Variance
#########################

getSumSquare <- function(value, arbre) {
 sum <- 0.
 for(eu in 1:dim(arbre$edge)[1])
   sum <- sum + (value[arbre$edge[eu,2]]-value[arbre$edge[eu,1]])^2/arbre$edge.length[eu]
 sum
}



####################################################
# NOUVELLES RECONSTRUCTION ML ET REML
####################################################

fillCoeffRec <- function(u, x, arbre, coeff) {
 child <- which(arbre$edge[, 1] == u)
 if(length(child)>0) {
  beH <- numeric(length(child))
  muH <- numeric(length(child))
  phH <- numeric(length(child))
  for(i in 1:length(child)) {
   v <- arbre$edge[child[i],2]
   if(length(which(arbre$edge[, 1] == v))>0) {
    coeff <- fillCoeffRec(v, x, arbre, coeff)
    beH[i] <- coeff$alpha[v]/(1+coeff$alpha[v]*arbre$edge.length[child[i]])
    muH[i] <- coeff$mu[v]
    phH[i] <- coeff$phi[v]
   } else {
    beH[i] <- 1/(arbre$edge.length[child[i]])
    muH[i] <- x[v]
    phH[i] <- 0
   }
  }
  coeff$alpha[u] <- 0.
  coeff$phi[u] <- 0.
  sum1 <- 0.
  sum2 <- 0.
  sum3 <- 0.
  for(i in 1:length(child)) {
   coeff$alpha[u] <- coeff$alpha[u] + beH[i]
   coeff$phi[u] <- coeff$phi[u] + phH[i]
   sum1 <- sum1 + beH[i]*muH[i]
   sum2 <- sum2 + beH[i]*muH[i]*muH[i]
   sum3 <- sum3 + beH[i]
  }
  coeff$mu[u] <- sum1/sum3
  coeff$phi[u] <- coeff$phi[u] + coeff$alpha[u]*((sum2/sum3)-coeff$mu[u]*coeff$mu[u])
}
coeff
}

fillCoeff <- function(x, arbre) {
 NN <- max(arbre[["edge"]])
 alpha <- numeric(NN)
 mu <- numeric(NN)
 phi <- numeric(NN)
 coeff <- list(alpha = alpha, mu = mu, phi = phi)
 root<-racine(arbre)
 child <- which(arbre$edge[, 1] == root)
 if(length(child)>0)
  for(i in 1:length(child))
   coeff <- fillCoeffRec(root, x, arbre, coeff)
 coeff
}

fillValuesRec <- function(eu, val, arbre, value, x, coeff) {
 u <- arbre$edge[eu,2]
 child <- which(arbre$edge[, 1] == u)
 if(length(child)>0) {
  value[u] <- (val+arbre$edge.length[eu]*coeff$alpha[u]*coeff$mu[u])/(1+arbre$edge.length[eu]*coeff$alpha[u])
  for(i in 1:length(child)) {
   v <- arbre$edge[child[i],2]
    value <- fillValuesRec(child[i], value[u], arbre, value, x, coeff)
  }
 } else
  value[u] = x[u]
 value
}


getMLHessian <- function(value, arbre) {
   sumSqu <- getSumSquare(value, arbre)
   nI <- nombreInternes(arbre)
   nT <- length(arbre$tip.label)
   nE <- nI+nT-1
   sizeH<-nI+1
   hessian <- matrix(0., nrow=sizeH, ncol=sizeH)
   var <- sumSqu/nE
   sd <- sqrt(var)
   hessian[1,1] <- -nE/(2*var^2)+sumSqu/var^3
   for(i in 1:nI) {
     	child <- which(arbre$edge[, 1] == nT+i)
  	if(length(child)>0) {
      		for(j in 1:length(child)) {
     		hessian[1,i+1] <- hessian[1,i+1]-(value[arbre$edge[child[j],2]]-value[nT+i])/arbre$edge.length[child[j]]
      		hessian[i+1,i+1] <- hessian[i+1,i+1]+1./arbre$edge.length[child[j]]
       		if(arbre$edge[child[j],2]>nT)
		  hessian[i+1,arbre$edge[child[j],2]-nT+1] <- -1./(var*arbre$edge.length[child[j]])
      		}
     	}
      	anc <- which(arbre$edge[, 2] == nT+i)
     	if(length(anc)>0) {
     	 	for(j in 1:length(anc)) {
       		hessian[1,i+1] <- hessian[1,i+1]+(value[nT+i]-value[arbre$edge[anc[j],1]])/arbre$edge.length[anc[j]]
     		hessian[i+1,i+1] <- hessian[i+1,i+1]+1./arbre$edge.length[anc[j]]
       		hessian[i+1,arbre$edge[anc[j],1]-nT+1] <- -1./(var*arbre$edge.length[anc[j]])
      		}
     	}
   hessian[1,i+1] <- -hessian[1,i+1]/sd^2
   hessian[i+1,1] <- hessian[1,i+1]
   hessian[i+1,i+1] <- hessian[i+1,i+1]/var
   }
   hessian 
}


newML <- function(x, arbre) {
 if(!inherits(arbre, "phylo"))
  stop("object \"phy\" is not of class \"phylo\"")
 if(is.null(arbre$edge.length))
  stop("tree has no branch lengths")
 nb.tip <- length(arbre$tip.label)
 nb.node <- arbre$Nnode
 if(nb.node != nb.tip - 1)
  stop("\"arbre\" is not rooted AND fully dichotomous.")
 if(length(x) != nb.tip)
  stop("length of phenotypic and of arbrelogenetic data do not match.")
 if(!is.null(names(x))) {
  if(all(names(x) %in% arbre$tip.label))
   x <- x[arbre$tip.label]
  else warning("the names of 'x' and the tip labels of the tree do not match: the former were ignored in the analysis.")
 }
 NN <- max(arbre[["edge"]])
 value <- numeric(NN)
 coeff <- fillCoeff(x, arbre)
 root<-racine(arbre)
 Qmin <- coeff$phi[root]
 value[root] <- coeff$mu[root]
 child <- which(arbre$edge[, 1] == root)
 if(length(child)>0)
  for(i in 1:length(child))
   value <- fillValuesRec(child[i], value[root], arbre, value, x, coeff)
 hessian <- getMLHessian(value, arbre)
 list(ace = value, se = sqrt(diag(solve(hessian))), sigma2 = getSumSquare(value, arbre)/(dim(arbre$edge)[1]), hessian = hessian)
}


getREMLHessian <- function(value, arbre, sigma2) {
   nI <- nombreInternes(arbre)
   nT <- length(arbre$tip.label)
   sizeH<-nI
   hessian <- matrix(0., nrow=sizeH, ncol=sizeH)
   for(i in 1:nI) {
     	child <- which(arbre$edge[, 1] == nT+i)
  	if(length(child)>0) {
      		for(j in 1:length(child)) {
      		hessian[i,i] <- hessian[i,i]+1./arbre$edge.length[child[j]]
       		if(arbre$edge[child[j],2]>nT)
		  hessian[i,arbre$edge[child[j],2]-nT] <- -1./(sigma2*arbre$edge.length[child[j]])
      		}
     	}
      	anc <- which(arbre$edge[, 2] == nT+i)
     	if(length(anc)>0) {
     	 	for(j in 1:length(anc)) {
      		hessian[i,i] <- hessian[i,i]+1./arbre$edge.length[anc[j]]
       		hessian[i,arbre$edge[anc[j],1]-nT] <- -1./(sigma2*arbre$edge.length[anc[j]])
      		}
     	}
      hessian[i,i] <- hessian[i,i]/sigma2
   }
   hessian 
}


newREML <- function(x, arbre) {
 if(!inherits(arbre, "phylo"))
  stop("object \"phy\" is not of class \"phylo\"")
 if(is.null(arbre$edge.length))
  stop("tree has no branch lengths")
 nb.tip <- length(arbre$tip.label)
 nb.node <- arbre$Nnode
 if(nb.node != nb.tip - 1)
  stop("\"arbre\" is not rooted AND fully dichotomous.")
 if(length(x) != nb.tip)
  stop("length of phenotypic and of arbrelogenetic data do not match.")
 if(!is.null(names(x))) {
  if(all(names(x) %in% arbre$tip.label))
   x <- x[arbre$tip.label]
  else warning("the names of 'x' and the tip labels of the tree do not match: the former were ignored in the analysis.")
 }
 NN <- max(arbre[["edge"]])
 value <- numeric(NN)
 coeff <- fillCoeff(x, arbre)
 root<-racine(arbre)
 Qmin <- coeff$phi[root]
 value[root] <- coeff$mu[root]
 child <- which(arbre$edge[, 1] == root)
 if(length(child)>0)
  for(i in 1:length(child))
   value <- fillValuesRec(child[i], value[root], arbre, value, x, coeff)
  minusLogLik <- function(sig2) {
    if (sig2 < 0) return(1e100)
    V <- sig2 * vcv(arbre)
    ## next three lines borrowed from dmvnorm() in 'mvtnorm'
    distval <- stats::mahalanobis(x, center = mu, cov = V)
    logdet <- sum(log(eigen(V, symmetric = TRUE, only.values = TRUE)$values))
    (nb.tip * log(2 * pi) + logdet + distval)/2
  }
  mu <- rep(ace(x, arbre, method="pic")$ace[1], nb.tip)
  out <- nlm(minusLogLik, 1, hessian = TRUE)
  sigma2 <- out$estimate
  se_sgi2 <- sqrt(1/out$hessian)
 hessian <- getREMLHessian(value, arbre, sigma2)
 list(ace = value, se = sqrt(diag(solve(hessian))), sigma2 = sigma2, se_sgi2 = se_sgi2, hessian = hessian)
}





####################################################
# GLS Reconstruction methods
####################################################


#GLS_BM
aceBM <- function (phy, x, CI=TRUE) 
{	
	obj <- list()
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	nbTotN <- nb.tip+nb.node
	sigmaMF <- 1	
	TsTemps <- dist.nodes(phy)
	TempsRacine <- TsTemps[(nb.tip+1),]
	IndicesMRCA <- mrca(arbre, full=T)
	M <- matrix(NA, ncol=nbTotN, nrow=nbTotN)
	for (i in 1:nbTotN)
	{
	for (j in 1:nbTotN)
	{
		M[i,j] <- sigmaMF^2 * TempsRacine[IndicesMRCA[i,j]]
	}
	}
# M = SigmaZ
	varAL <- M[-(1:nb.tip), 1:nb.tip]
	varAA <- M[-(1:nb.tip), -(1:nb.tip)]
        varLL <- M[(1:nb.tip), 1:nb.tip]
        invVarLL <- solve(varLL)
	UL <- rep(1, length=nb.tip)
	UA <- rep(1, length=nb.node)
	TL <- TempsRacine[1:nb.tip] 
	TA <- TempsRacine[(nb.tip+1):(nb.tip+nb.node)]
#
	IVL_Z <- invVarLL %*% x
	IVL_T <- invVarLL %*% TL
	IVL_U <- invVarLL %*% UL
	U_IVL_U <- t(UL) %*% IVL_U
	U_IVL_Z <- t(UL) %*% IVL_Z
	DeltaU <- UA - varAL %*% IVL_U
#
	Racine_chap <- solve(U_IVL_U) %*% U_IVL_Z 
	Racine_chap <- as.numeric(Racine_chap)
	Anc_chap <- Racine_chap * DeltaU + varAL %*% IVL_Z
	Anc_chap <- as.vector(Anc_chap)
	obj$ace <- Anc_chap
	names(obj$ace) <- (nb.tip + 1):(nb.tip + nb.node)
#
 	if (CI) {
		Vec <- x - Racine_chap 
		Num <- t(Vec) %*% invVarLL %*% Vec
		Num <- as.numeric(Num)
		sigma2_chap <- Num / (nb.tip-1)
		obj$sigma2 <- sigma2_chap 
                se <- sqrt((varAA - varAL %*% invVarLL %*% t(varAL))[cbind(1:nb.node, 
                  1:nb.node)])
		se <- se * sqrt(sigma2_chap)
                tmp <- se * qnorm(0.025)
                obj$CI95 <- cbind(obj$ace + tmp, obj$ace - tmp)
            }
	obj
}

#GLS_ABM
aceABM <- function (phy, x, CI=TRUE) 
{	
	obj <- list()
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	nbTotN <- nb.tip+nb.node
	sigmaMF <- 1	
	TsTemps <- dist.nodes(phy)
	TempsRacine <- TsTemps[(nb.tip+1),]
	IndicesMRCA <- mrca(arbre, full=T)
	M <- matrix(NA, ncol=nbTotN, nrow=nbTotN)
	for (i in 1:nbTotN)
	{
	for (j in 1:nbTotN)
	{
		M[i,j] <- sigmaMF^2 * TempsRacine[IndicesMRCA[i,j]]
	}
	}
# M = SigmaZ
	varAL <- M[-(1:nb.tip), 1:nb.tip]
	varAA <- M[-(1:nb.tip), -(1:nb.tip)]
        varLL <- M[(1:nb.tip), 1:nb.tip]
        invVarLL <- solve(varLL)
	UL <- rep(1, length=nb.tip)
	UA <- rep(1, length=nb.node)
	TL <- TempsRacine[1:nb.tip] 
	TA <- TempsRacine[(nb.tip+1):(nb.tip+nb.node)]
#
	IVL_Z <- invVarLL %*% x
	IVL_T <- invVarLL %*% TL
	IVL_U <- invVarLL %*% UL
	U_IVL_U <- t(UL) %*% IVL_U
	T_IVL_T <- t(TL) %*% IVL_T
	U_IVL_T <- t(UL) %*% IVL_T
	U_IVL_Z <- t(UL) %*% IVL_Z
	T_IVL_Z <- t(TL) %*% IVL_Z
	DeltaT <- TA - varAL %*% IVL_T
	DeltaU <- UA - varAL %*% IVL_U
#
	Den <- U_IVL_U * T_IVL_T - U_IVL_T^2
	Den <- as.numeric(Den)
	Mu_chap <- (U_IVL_U * T_IVL_Z - U_IVL_T * U_IVL_Z) / Den
	Mu_chap <- as.numeric(Mu_chap)
	Racine_chap <- (T_IVL_T * U_IVL_Z - U_IVL_T * T_IVL_Z) / Den
	Racine_chap <- as.numeric(Racine_chap)
	Anc_chap <- Mu_chap * DeltaT + Racine_chap * DeltaU + varAL %*% IVL_Z
	Anc_chap <- as.vector(Anc_chap)
	obj$ace <- Anc_chap
	names(obj$ace) <- (nb.tip + 1):(nb.tip + nb.node)
	obj$mu <- Mu_chap
#
 	if (CI) {
		Vec <- x - Racine_chap - Mu_chap * TL
		Num <- t(Vec) %*% invVarLL %*% Vec
		Num <- as.numeric(Num)
		sigma2_chap <- Num / (nb.tip-1)
		obj$sigma2 <- sigma2_chap 
                se <- sqrt((varAA - varAL %*% invVarLL %*% t(varAL))[cbind(1:nb.node, 
                  1:nb.node)])
		se <- se * sqrt(sigma2_chap)
                tmp <- se * qnorm(0.025)
                obj$CI95 <- cbind(obj$ace + tmp, obj$ace - tmp)
            }
	obj
}

#GLS_OU* (theta = z0)
aceOUv1 <- function (phy, x, alpha, CI=TRUE) 
{	
	obj <- list()
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	nbTotN <- nb.tip+nb.node
	sigmaMF <- 1
	alphaM <- alpha
	nbTotN <- nb.tip+nb.node
	TsTemps <- dist.nodes(phy)
	TempsRacine <- TsTemps[(nb.tip+1),]
	IndicesMRCA <- mrca(arbre, full=T)
	M <- matrix(NA, ncol=nbTotN, nrow=nbTotN)
	for (i in 1:nbTotN)
		{
	for (j in 1:nbTotN)
		{
		Tempsm <- TempsRacine[IndicesMRCA[i,j]]
		Tempsi <- TempsRacine[i]
		Tempsj <- TempsRacine[j]
		M[i,j] <- sigmaMF^2 * exp(-alphaM * (Tempsi+Tempsj-2*Tempsm)) * (1-exp(-2*alphaM * Tempsm)) / (2 * alphaM)
		}
		}
# M = SigmaZ
	varAL <- M[-(1:nb.tip), 1:nb.tip]
	varAA <- M[-(1:nb.tip), -(1:nb.tip)]
        varLL <- M[(1:nb.tip), 1:nb.tip]
        invVarLL <- solve(varLL)
	UL <- rep(1, length=nb.tip)
	UA <- rep(1, length=nb.node)
	TL <- TempsRacine[1:nb.tip] 
	TA <- TempsRacine[(nb.tip+1):(nb.tip+nb.node)]
#
	IVL_Z <- invVarLL %*% x
	IVL_T <- invVarLL %*% TL
	IVL_U <- invVarLL %*% UL
	U_IVL_U <- t(UL) %*% IVL_U
	U_IVL_Z <- t(UL) %*% IVL_Z
	DeltaU <- UA - varAL %*% IVL_U
#
	Racine_chap <- solve(U_IVL_U) %*% U_IVL_Z 
	Racine_chap <- as.numeric(Racine_chap)
	Anc_chap <- Racine_chap * DeltaU + varAL %*% IVL_Z
	Anc_chap <- as.vector(Anc_chap)
	obj$ace <- Anc_chap
	names(obj$ace) <- (nb.tip + 1):(nb.tip + nb.node)
#
# vraisemblance
#
	mL <- Racine_chap
	Num <- t(x-mL) %*% invVarLL %*% (x-mL)
	Num <- as.numeric(Num)
	sigma2_chap <- Num / (nb.tip-1)
	obj$sigma <- sqrt(sigma2_chap)
	VL <- sigma2_chap * varLL
	invVL <- invVarLL / sigma2_chap
	LVrais <- - t(x-mL) %*% invVL %*% (x-mL) /2 - nb.tip * log(2*pi)/2 - log(det(VL))/2
	obj$LLik <- as.numeric(LVrais)  
#
 	if (CI) {
                se <- sqrt((varAA - varAL %*% invVarLL %*% t(varAL))[cbind(1:nb.node, 
                  1:nb.node)])
		se <- se * sqrt(sigma2_chap)
                tmp <- se * qnorm(0.025)
                obj$CI95 <- cbind(obj$ace + tmp, obj$ace - tmp)
            }
	obj
}


#GLS_OU
aceOUv2 <- function (phy, x, alpha, CI=TRUE) 
{	
	obj <- list()
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode
	nbTotN <- nb.tip+nb.node
	sigmaMF <- 1
	nbTotN <- nb.tip+nb.node
	TsTemps <- dist.nodes(phy)
	TempsRacine <- TsTemps[(nb.tip+1),]
	IndicesMRCA <- mrca(arbre, full=T)
	M <- matrix(NA, ncol=nbTotN, nrow=nbTotN)
	for (i in 1:nbTotN)
		{
	for (j in 1:nbTotN)
		{
		Tempsm <- TempsRacine[IndicesMRCA[i,j]]
		Tempsi <- TempsRacine[i]
		Tempsj <- TempsRacine[j]
		M[i,j] <- sigmaMF^2 * exp(-alpha * (Tempsi+Tempsj-2*Tempsm)) * (1-exp(-2*alpha * Tempsm)) / (2 * alpha)
		}
		}
# M = SigmaZ
	varAL <- M[-(1:nb.tip), 1:nb.tip]
	varAA <- M[-(1:nb.tip), -(1:nb.tip)]
        varLL <- M[(1:nb.tip), 1:nb.tip]
        invVarLL <- solve(varLL)
	vecW <- exp(-alpha * TempsRacine)
	UL <- vecW[1:nb.tip] 
	UA <- vecW[(nb.tip+1):(nb.tip+nb.node)]
	TL <- 1-UL
	TA <- 1-UA
#
#
	IVL_Z <- invVarLL %*% x
	IVL_T <- invVarLL %*% TL
	IVL_U <- invVarLL %*% UL
	U_IVL_U <- t(UL) %*% IVL_U
	T_IVL_T <- t(TL) %*% IVL_T
	U_IVL_T <- t(UL) %*% IVL_T
	U_IVL_Z <- t(UL) %*% IVL_Z
	T_IVL_Z <- t(TL) %*% IVL_Z
	DeltaT <- TA - varAL %*% IVL_T
	DeltaU <- UA - varAL %*% IVL_U
#
	Den <- U_IVL_U * T_IVL_T - U_IVL_T^2
	Den <- as.numeric(Den)
	Theta_chap <- (U_IVL_U * T_IVL_Z - U_IVL_T * U_IVL_Z) / Den
	Theta_chap <- as.numeric(Theta_chap)
	Racine_chap <- (T_IVL_T * U_IVL_Z - U_IVL_T * T_IVL_Z) / Den
	Racine_chap <- as.numeric(Racine_chap)
	Anc_chap <- Theta_chap * DeltaT + Racine_chap * DeltaU + varAL %*% IVL_Z
	Anc_chap <- as.vector(Anc_chap)
	obj$ace <- Anc_chap
	names(obj$ace) <- (nb.tip + 1):(nb.tip + nb.node)
	obj$theta <- Theta_chap
#
# vraisemblance
#
	mL <- (Racine_chap * UL + Theta_chap * TL)
	Num <- t(x-mL) %*% invVarLL %*% (x-mL)
	Num <- as.numeric(Num)
	sigma2_chap <- Num / (nb.tip-1)
	obj$sigma <- sqrt(sigma2_chap)
	VL <- sigma2_chap * varLL
	invVL <- invVarLL / sigma2_chap
	LVrais <- - t(x-mL) %*% invVL %*% (x-mL) /2 - nb.tip * log(2*pi)/2 - log(det(VL))/2
	obj$LLik <- as.numeric(LVrais)  
#
 	if (CI) {
                se <- sqrt((varAA - varAL %*% invVarLL %*% t(varAL))[cbind(1:nb.node, 
                  1:nb.node)])
		se <- se * sqrt(sigma2_chap)
                tmp <- se * qnorm(0.025)
                obj$CI95 <- cbind(obj$ace + tmp, obj$ace - tmp)
            }
	obj
}




#######################################################
#######################################################
# Integration numerique des fonctions de repartition (NRJ)
#######################################################
#######################################################


EnergyGaussianVsGaussian1 <- function(mu1, si1, mu2, si2) {
  (2./sqrt(pi)) * (sqrt(2.)*sqrt(si1^2+si2^2)*exp(-(mu2-mu1)^2/(2.*(si1^2+si2^2))) -si1-si2)+2.*(mu2-mu1)*(1-2.*pnorm((mu1-mu2)/sqrt(si1^2+si2^2)))
}


# Distance non normalisee
###########################

#
EnergyGaussianVsGaussian <- function(mu1, si1, mu2, si2) {
  GG <-function(x) {
	((pnorm(x, m=mu1, s=si1)-pnorm(x, m=mu2, s=si2))^2)*2.
	}
	borne <- max( abs(mu1-4*si1), abs(mu1+4*si1), abs(mu2-4*si2), abs(mu2+4*si2) )
	if ( sum(is.na(mu1-mu2))==0 )
		integrate(GG, lower=-borne, upper=borne, stop.on.error=F)$value
	else
		NA
}
#


# Formule explicite
#
EspGauss <- function(mu1, si1){
si1*sqrt(2./pi)*exp(-mu1^2/(2.*si1^2))+mu1*(1-2*pnorm(-mu1/si1))
}
#
EnergyGaussianVsGaussian2 <- function(mu1, si1, mu2, si2) {
( 2 * EspGauss(mu1=mu1-mu2, si1=sqrt(si1^2+si2^2)) - 2*si1/sqrt(pi) - 2*si2/sqrt(pi) ) 
}


EnergyGaussianVsStudent <- function(mu1, si1, mu2, si2, k) {
  GS <-function(x) {
	((pnorm(x, m=mu1, s=si1)-pt((x-mu2)/si2,df=k))^2)*2.
	}
	qualite <- pnorm(4) - pnorm(-4)
	quantile <- qt( (1+qualite)/2, df=k)
	borne <- max( abs(mu1-4*si1), abs(mu1+4*si1), abs(mu2-quantile*si2), abs(mu2+quantile*si2) )
	if ( sum(is.na(mu1-mu2))==0 )
		integrate(GS, lower=-borne, upper=borne, stop.on.error=F)$value
	else
		NA
	}


EnergyConstantVsStudent <- function(mu2, si2, cste, k) {
  CS <-function(x) {
	((pnorm(x, m=cste, s=0)-pt((x-mu2)/si2,df=k))^2)*2.
}
	qualite <- pnorm(4) - pnorm(-4)
	quantile <- qt( (1+qualite)/2, df=k)
	borne <- max( abs(mu2-quantile*si2), abs(mu2+quantile*si2), abs(cste) )
	if ( sum(is.na(mu1-mu2))==0 )
		integrate(CS, lower=-borne, upper=borne, stop.on.error=F)$value
	else
		NA
}



# Distance NON normalisee
############################



EnergyGaussianVsGaussiannN <- function(mu1, si1, mu2, si2) {
EnergyGaussianVsGaussian1(mu1=mu1, si1=si1, mu2=mu2, si2=si2)
}

EnergyConstantVsConstantnN <- function(mu1, mu2) {
if ( !is.na(mu1-mu2) )
{
	if (mu1 == mu2)
		val <- 0
	if (mu1 != mu2)
		val <- 2* abs(mu1-mu2)
}
if ( is.na(mu1-mu2) )
		val <- NA	
return(val)
}


EspStud <- function(mu2, si2, k) {
	(2*si2 / sqrt(pi)) * (sqrt(k) / (k-1)) * (gamma((k+1)/2.)/gamma(k/2.)) * ( 1+ mu2^2 / (k*si2^2) )^((1-k)/2) + mu2*(1-2.*pt(-mu2/si2, df=k))
}


EnergyGaussianVsStudentnN <- function(mu1, si1, mu2, si2, k) {
EnergyGaussianVsStudent(mu1=mu1, si1=si1, mu2=mu2, si2=si2, k=k)  
}



EnergyConstantVsStudentnN <- function(mu2, si2, cste, k) {
EnergyConstantVsStudent(mu2=mu2, si2=si2, cste=cste, k=k) 
}





#######################################################
# Calcul de l'erreur (NRJ non normalisee)
#######################################################

zero <- 1e-9
NRJ <- function(muR, muT, sigma2R, sigma2T, loiR, ddl)
{
indNAs <- is.na(sigma2R)
indNAm <- is.na(muR)
n <- length(muR)
valT <- rep(NA, length=n)
#
# NRJ pour 2 lois normales
#
for (i in 1:n) 
{
if ((!indNAs[i]) & (!indNAm[i]))
{
	if ( (sigma2R[i] > zero) & (sigma2T[i] > zero) ) 
	{
		if (loiR[i]=="Norm") 
		{
		valNRJ <- EnergyGaussianVsGaussiannN(mu1=muT[i], si1=sqrt(sigma2T[i]), mu2=muR[i], si2=sqrt(sigma2R[i]))
		}
		#
		#
		if (loiR[i]=="Stud") 
		{
		valNRJ <- EnergyGaussianVsStudentnN(mu1=muT[i], si1=sqrt(sigma2T[i]), mu2=muR[i], si2=sqrt(sigma2R[i]), k=ddl)
		}
		valT[i] <- valNRJ
	}
#
# lois theorique et reconstruite = constantes
#
	if ( (sigma2T[i] <= zero) & (sigma2R[i] <= zero) ) 
	{
	valT[i] <- EnergyConstantVsConstantnN(mu1=muT[i], mu2=muR[i])
	}
	if ( (sigma2R[i] <= zero) & (sigma2T[i] > zero) ) 
	{
	valT[i] <- EnergyGaussianVsGaussiannN(mu1=muT[i], si1=sqrt(sigma2T[i]), mu2=muR[i], si2=0)	
	}
	if ( (sigma2T[i] <= zero) & (sigma2R[i] > zero) ) 
	{
	if (loiR[i]=="Stud") 
		{
		valT[i] <- EnergyConstantVsStudentnN(mu2=muR[i], si2=sqrt(sigma2R[i]), cste=muT[i], k=ddl)
		}
	if (loiR[i]=="Norm") 
		{
		valT[i] <- EnergyGaussianVsGaussiannN(mu1=muT[i], si1=0, mu2=muR[i], si2=sqrt(sigma2R[i]))
		}
	}
}
}
return(valT)
}



#######################################################
#######################################################
# (KS)
#######################################################
#######################################################



KSStudentVsGaussian <- function(mu1, si1, mu2, si2, k, plotc) {
	zero <- 1e-3
	inc <- 100
	maxIter <- 100
	po <- function(x) {
		-x^3+(mu1+2*mu2)*x^2+x*((k+1)*si1^2-2*mu2*mu1-k*si2^2-mu2^2)+mu1*(k*si2^2+mu2^2)-(k+1)*mu2*si1^2
	}
	logfh <- function(x) {
		dnorm(x, mean = mu1, sd = si1, log = TRUE) - dt((x-mu2)/si2, k, log = TRUE) + log(si2)
	}
	dd <- function(x) {
		abs(pnorm(x, mean = mu1, sd = si1)-pt((x-mu2)/si2, k))
	}
# corrige dt en divisant pas si2
	pseudo <- function(x) {
		-(x-mu1)^2/(2*si1^2)+((k+1)/2)*log(1+(x-mu2)^2/(k*si2^2))
	}

	v <- sort(polyroot(c(mu1*(k*si2^2+mu2^2)-(k+1)*mu2*si1^2, (k+1)*si1^2-2*mu2*mu1-k*si2^2-mu2^2, mu1+2*mu2, -1)))

#keep only real roots
	j <- 0
	for(i in 1:length(v)) {
		if(abs(Im(v[i]))<=zero) {
			j <- j+1
			v[j] <- Re(v[i])
		}
	}
	length(v) <- j
	v <- Re(v)
 #keep only roots which change sign
	sgp = 1;
	j <- 0
	if(length(v)>1)
		for(i in 1:(length(v)-1)) {
			if(sgp*po((v[i]+v[i+1])/2) < 0) {
				j <- j+1
				v[j] <- v[i]
				sgp <- -sgp
			}
#			sgp <- po((v[i]+v[i+1])/2)/abs(po((v[i]+v[i+1])/2))

		}
	if(sgp > 0) {
		j <- j+1
		v[j] <- v[length(v)]
	}
	length(v) <- j
	max <- 0.
	xmax = 0.

	binf <- v[1]-inc
	iter <-1
	while( (is.finite(logfh(binf))) & !is.nan(logfh(binf)) & (logfh(binf)>=0) & iter<maxIter) {
		binf <- binf-inc
		iter <- iter+1
	}
	if(iter<maxIter & is.finite(logfh(binf)) & !is.nan(logfh(binf))) {
		for(i in 1:length(v)) {
			bsup <- v[i]
			if(logfh(binf)*logfh(bsup)<0) {
				x <- uniroot(logfh, lower = binf, upper = bsup, tol = 1e-9)
				if(dd(x$root)>max) {
					max = dd(x$root)
					xmax = x$root
				}
			} else {
				if(abs(logfh(bsup))<zero) {
					if(dd(bsup)>max) {
						max = dd(bsup)
						xmax = bsup
					}
				}
			}
			binf = bsup
		}
		iter <-1
		bsup <- bsup+inc
		while( (is.finite(logfh(bsup))) & !is.nan(logfh(bsup)) & (logfh(bsup)>=0)  & iter<maxIter) {
			bsup <- bsup+inc
			iter <- iter+1
		}
		if(iter<maxIter & is.finite(logfh(bsup)) & !is.nan(logfh(bsup))) {   
			if(logfh(binf)*logfh(bsup)<0) {
				x <- uniroot(logfh, lower = binf, upper = bsup, tol = 1e-9)
				if(dd(x$root)>max) {
					max = dd(x$root)
					xmax = x$root
				}
			}
			if (plotc==TRUE) {
				par(mfrow=c(2,2))
				minc <- -50
				maxc <- 100
				curve(po, from = minc, to = maxc)
				abline(h=0)
				for(i in 1:length(v))
					abline(v=v[i], lty=2)
				curve(pseudo, from = minc, to = maxc)
				abline(h=0)
				for(i in 1:length(v))
					abline(v=v[i], lty=2)
				curve(logfh, from = minc, to = maxc)
				abline(h=0)
				for(i in 1:length(v))
					abline(v=v[i], lty=2)
				curve(dd, from = minc, to = maxc)
				abline(v=xmax)
				par(mfrow=c(1,1))
			}
		} else {
			print("NaN")
			print(logfh(bsup))
			max = NaN
		}
	} else {
		print("NaN")
		print(logfh(binf))
		max = NaN
	}
	max
}


KS <- function(muR, muT, sigma2R, sigma2T, loiR, ddl)
{
zero <- 1e-9
indNAs <- is.na(sigma2R)
indNAm <- is.na(muR)
n <- length(muR)
valT <- rep(NA, length=n)
#
# KS pour 2 lois normales
#
for (i in 1:n) 
{
if ((!indNAs[i]) & (!indNAm[i])) 
{
	if ( (sigma2R[i] > zero) & (sigma2T[i] > zero) ) 
	{
		if (loiR[i]=="Norm") 
		{
		cc <- sigma2T[i]*sigma2R[i] * log(sigma2T[i]/sigma2R[i]) + (muT[i]^2*sigma2R[i] - muR[i]^2*sigma2T[i])
		cb <- 2 * (muR[i]*sigma2T[i] - muT[i]*sigma2R[i])
		ca <- sigma2R[i] - sigma2T[i]
		delta <- cb^2 - 4*ca*cc
			if ( (delta>=0) & (ca!=0) )
			{
			r1 <- (-cb + sqrt(delta)) / (2*ca)
			r2 <- (-cb - sqrt(delta)) / (2*ca)
			KSr1 <- abs( pnorm(r1, m=muR[i], s=sqrt(sigma2R[i])) - pnorm(r1,  m=muT[i], s=sqrt(sigma2T[i])) )
			KSr2 <- abs( pnorm(r2, m=muR[i], s=sqrt(sigma2R[i])) - pnorm(r2,  m=muT[i], s=sqrt(sigma2T[i])) )
			valKS <- max(KSr1, KSr2)
			}
			if ( (ca==0) & (cb!=0) )
			{
			r12 <- -cc/cb
			valKS <- abs( pnorm(r12, m=muR[i], s=sqrt(sigma2R[i])) - pnorm(r12,  m=muT[i], s=sqrt(sigma2T[i])) )
			}
			if ( (ca==0) & (cb==0) )
			{
			valKS <- 0
			}
			if (delta<0)
			{
			valKS <- NA
			}
		}
		#
		#
		if (loiR[i]=="Stud") 
		{
		valKS <- KSStudentVsGaussian(mu1=muT[i], si1=sqrt(sigma2T[i]), mu2=muR[i], si2=sqrt(sigma2R[i]), k=ddl, plotc=F )
		}
		if (is.nan(valKS)) 
		{
		valKS <- NA
		}
		valT[i] <- valKS
	}
#
# loi reconstruite = cste (gls à la racine)
#
	if ( (sigma2R[i] <= zero ) & (sigma2T[i] > zero) ) 
	{
	pc <- pnorm( (muR[i]-muT[i])/ sqrt(sigma2T[i]) )
	valT[i] <- max( pc, 1-pc )
	}
#
# loi theorique à la racine = cste 
#
	if ( (sigma2T[i] <= zero) & (sigma2R[i] > zero) ) 
	{
	if (loiR[i]=="Norm") 
	{
		pc <- pnorm( (muT[i]-muR[i])/ sqrt(sigma2R[i]) )
		valT[i] <- max( pc, 1-pc )
	}
	if (loiR[i]=="Stud") 
	{
		pc <- pt( (muT[i]-muR[i])/ sqrt(sigma2R[i]), df=ddl )
		valT[i] <- max( pc, 1-pc )
	}
	}

#
# lois theorique et reconstruite = constantes
#
	if ( (sigma2T[i] <= zero) & (sigma2R[i] <= zero) ) 
	{
		if (abs(muT[i]-muR[i]) <= zero)
		{
		valT[i] <- 0
		}
		if (abs(muT[i]-muR[i]) > zero)
		{
		valT[i] <- 1
		}
	}	
}
}
return(valT)
}



#######################################################
#######################################################
# ESPERANCE QUADRATIQUE
#######################################################
#######################################################



ESPq <- function(muR, muT, sigma2R, sigma2T)
{
#sigma2R[is.na(sigma2R)] <- Inf
#sigma2T[is.na(sigma2T)] <- Inf
valT <- sigma2R + sigma2T + (muR-muT)^2
return(valT)
}



###########################
# arbre renumerote
#########################

new <- read.table(nomFichiertxt)
new <- as.character(new[1,1])
conv <- newick2phylog(new)
arbre <- as.phylo(conv)
nT <- length(arbre$tip.label)
nN <- arbre$Nnode

ecritarbre(file("ancien.newick", open = "w"), arbre)
transf <- renumeroteArbre(arbre)
arbis <- transf$arbre
ecritarbre(file("nouveau.newick", open = "w"), arbis)
vec <- transf$cod[(nT+1):(nT+nN)]
#
nomNoeud <- as.character(1:nN)


# Longueur de branche dans bon ordre des accroissements
#
LongBranche <- arbis$edge.length
Fils <- arbis$edge[,2]
arbis$edge <- arbis$edge[order(Fils),]
arbis$edge.length <- arbis$edge.length[order(Fils)]
Tau <- arbis$edge.length
LongMin <- min(Tau)
LongMax <- max(Tau)








###########################################################"
# Simulations
############################################################


Model <- c("ABM", "OU")
nM <- length(Model)
Im <- 2
#

#
#GLSOU, GLSABM
Methodes <- c("REMLexact", "pic", "GLSLin", "MLexact", "Parsexact", "GLSOUop", "GLSOU", "GLSABM")
nMe <- length(Methodes)
#



ValMuRec <- matrix(0, nrow=N*nSM*naM*nbM, ncol=nMe*nN+3)
colnames(ValMuRec) <- c(paste(rep(Methodes, each=nN), nomNoeud, sep="_"), "aM", "bM", "sigmaM")
#
ValSigmaRec <- matrix(0, nrow=N*nSM*naM*nbM, ncol=nMe*nN+3)
colnames(ValSigmaRec) <- c(paste(rep(Methodes, each=nN), nomNoeud, sep="_"), "aM", "bM", "sigmaM")
#
ErreursOUnrj <- matrix(NA, nrow=N*nSM*naM*nbM, ncol=nMe*nN+3)
colnames(ErreursOUnrj) <- c(paste(rep(Methodes, each=nN), nomNoeud, sep="_"), "aM", "bM", "sigmaM")
#
ErreursOUnrjM <- matrix(NA, nrow=N*nSM*naM*nbM, ncol=nMe*nN+3)
colnames(ErreursOUnrjM) <- c(paste(rep(Methodes, each=nN), nomNoeud, sep="_"), "aM", "bM", "sigmaM")
#
ErreursOUks <- matrix(NA, nrow=N*nSM*naM*nbM, ncol=nMe*nN+3)
colnames(ErreursOUks) <- c(paste(rep(Methodes, each=nN), nomNoeud, sep="_"), "aM", "bM", "sigmaM")
#
ErreursOUmse <- matrix(NA, nrow=N*nSM*naM*nbM, ncol=nMe*nN+3)
colnames(ErreursOUmse) <- c(paste(rep(Methodes, each=nN), nomNoeud, sep="_"), "aM", "bM", "sigmaM")
#
ErreursOUbiaisG <- matrix(NA, nrow=N*nSM*naM*nbM, ncol=nMe*nN+3)
colnames(ErreursOUbiaisG) <- c(paste(rep(Methodes, each=nN), nomNoeud, sep="_"), "aM", "bM", "sigmaM")
#
ErreursOUbiaisAbsG <- matrix(NA, nrow=N*nSM*naM*nbM, ncol=nMe*nN+3)
colnames(ErreursOUbiaisAbsG) <- c(paste(rep(Methodes, each=nN), nomNoeud, sep="_"), "aM", "bM", "sigmaM")

 
#
#
#
for (ia in 1:naM) 
{
print(paste("ia=", ia))

for (ib in 1:nbM) 
{  
print(paste("ib=", ib))
for (j in 1:nSM) 
{
#
# Calcul de muTheo et sigmaTheo au noeud en fonction de ceux du modèle !!
# on suppose racine constante
# 
	if (root == "constant")
	{
	z0 <- m0
	if (aM[ia]==0)
	{
	muAcc <- rep(0,length(Tau))
	sigmaAcc <- sigmaM[j]^2 * Tau
	V2Acc <- diag(sigmaAcc)
  	w <- rep(1, length=(nN+nT-1))
	CABM <- calculeC("ABM", arbis)
	C <- CABM[2:(nN+nT), 2:(nN+nT)]
	}
	if (aM[ia]>0)
	{
	muAcc <- bM[ib] * (1-exp(- aM[ia] *Tau))
	sigmaAcc <- sigmaM[j]^2 * (1-exp(- 2 * aM[ia] *Tau)) / (2 * aM[ia])
	V2Acc <- diag(sigmaAcc)
  	w <- calculeV_exptimes(arbis, aM[ia])
	w <- w[-1]
  	COU <- calculeC("OU", arbis, aM[ia])
	C <- COU[2:(nN+nT), 2:(nN+nT)]
	}
	muNoeuds <- C %*% muAcc + z0*w
	sigmaNoeuds <- C %*% V2Acc %*% t(C)
	mu1 <- muNoeuds[1:(nN-1)]
	mu2 <- muNoeuds[nN:(nN+nT-1)]
	sigma11 <- sigmaNoeuds[(1:(nN-1)), (1:(nN-1))]
	sigma22 <- sigmaNoeuds[(nN:(nN+nT-1)) , (nN:(nN+nT-1)) ]
	sigma12 <- sigmaNoeuds[(1:(nN-1)), (nN:(nN+nT-1)) ]
	sigma21 <- sigmaNoeuds[(nN:(nN+nT-1)) , (1:(nN-1))]
	}
#
#	
# Calcul de muTheo et sigmaTheo au noeud en fonction de ceux du modèle !!
# on suppose racine suit une loi normale
# 
	if (root == "gaussian")
	{
	z0 <- rnorm(1, m=m0, sd=v0)
	if (aM[ia]==0)
	{
	muAcc <- rep(0,length(Tau))
	sigmaAcc <- sigmaM[j]^2 * Tau
	V2Acc <- diag(sigmaAcc)
  	w <- rep(1, length=(nN+nT))
	CABM <- calculeC("ABM", arbis)
	C <- CABM
	}
	if (aM[ia]>0)
	{
	muAcc <- c(m0, bM[ib] * (1-exp(- aM[ia] *Tau)))
	sigmaAcc <- c(v0, sigmaM[j]^2 * (1-exp(- 2 * aM[ia] *Tau)) / (2 * alpha[ia]))
	V2Acc <- diag(sigmaAcc)
	COU <- calculeC("OU", arbis)
	# verifier colonne 1 que des 1
	C <- COU
	}
	muNoeuds <- C %*% muAcc 
	sigmaNoeuds <- C %*% V2Acc %*% t(C)
 	mu1 <- muNoeuds[1:nN]
	mu2 <- muNoeuds[(nN+1):(nN+nT)]
	sigma11 <- sigmaNoeuds[1:nN, 1:nN]
	sigma22 <- sigmaNoeuds[((nN+1):(nN+nT)) , ((nN+1):(nN+nT)) ]
	sigma12 <- sigmaNoeuds[1:nN, ((nN+1):(nN+nT)) ]
	sigma21 <- sigmaNoeuds[((nN+1):(nN+nT)) , 1:nN]
	}
#
#
	nbRec <- 1
	k <- 1
	while (k <= N)
	{
	if (aM[ia]==0)
	{
	simul <- rTraitCont(arbre, model="BM", sigma=sigmaM[j], ancestor=T,  root.value=z0)
	}
	if (aM[ia]>0)
	{
	simul <- rTraitCont2(arbre, model="OU", theta=bM[ib], alpha=aM[ia], sigma=sigmaM[j], ancestor=T,  root.value=z0)
	}

	EtatTheo <- simul[(nT+1):(nN+nT)]
	EtatFeuilles <- simul[1:nT]
	EtatFeuillesN <- EtatFeuilles
#
	for (iF in 1:nT)
	{
	indT <- transf$cod[iF] - nN 
	EtatFeuillesN[ indT ] <- EtatFeuilles[iF]
	names(EtatFeuillesN)[indT]  <- names(EtatFeuilles)[iF]
	}
	ParadisToNous <- function(x)
	{
		NumP <- as.numeric( names(x) )
		res <- x
		for (iPN in 1:nN)
		{
			res[ transf$cod[NumP[iPN]] ] <- x[iPN]
			names(res) <- as.character(1:nN)
		}
		res
		
	}
#
	EtatTheoN <- ParadisToNous(EtatTheo)
#	
	ParadisToNousBis <- function(x, cod) {
		res <- numeric(length(x))
		if(length(x)>0)
		  for (iPN in 1:length(x))
			res[cod[iPN]] <- x[iPN]
		res
	}
#
	muCond <- mu1 + sigma12 %*% solve(sigma22) %*% (EtatFeuillesN-mu2)
	sigmaCond <- sigma11 - sigma12 %*% solve(sigma22) %*% sigma21
#
	indL <- (ia-1)*nbM*nSM*N + (ib-1)*nSM*N + (j-1)*N +k
	ValMuRec[indL, ((nMe*nN+1):(nMe*nN+3)) ] <- c(aM[ia], bM[ib], sigmaM[j])
	ValSigmaRec[indL, ((nMe*nN+1):(nMe*nN+3)) ] <- c(aM[ia], bM[ib], sigmaM[j])
	#
	#
	recREML <- newREML(EtatFeuilles, arbre)
	ValMuRec[indL, (1:nN) ] <- ParadisToNousBis(recREML$ace, transf$cod)[1:nN]
	ectR <- recREML$se
  	ectRn <- ectR
  	names(ectRn) <- as.character( (nT+1):(nT+nN) )
	ValSigmaRec[indL, (1:nN) ] <- ParadisToNous(ectRn)
	ectRcomp <- ectR
	#
	#
	recpic <- ace(EtatFeuilles, arbre, type="continuous", method="pic")
	ValMuRec[indL, ((nN+1):(2*nN)) ] <- ParadisToNous(recpic$ace)
	ectR <- (recpic$CI95[,2] - recpic$ace) / qnorm(0.975)
	ValSigmaRec[indL, ((nN+1):(2*nN)) ] <- ParadisToNous(ectR)      
	ectRcomp <- c(ectRcomp, ectR)
	#
	recGLSLin <- aceBM(phy=arbre, x=EtatFeuilles, CI=T)
	ValMuRec[indL, ((2*nN+1):(3*nN)) ] <- ParadisToNousBis(recREML$ace, transf$cod)[1:nN]
	ectR <- (recGLSLin$CI95[,2] - recGLSLin$ace) / qnorm(0.975)
	ValSigmaRec[indL, ((2*nN+1):(3*nN)) ] <- ParadisToNous(ectR)
	ectRcomp <- c(ectRcomp, ectR)
	# 
	recOpt <- newML(EtatFeuilles, arbre)
	ValMuRec[indL, ((3*nN+1):(4*nN)) ] <- ParadisToNousBis(recOpt$ace, transf$cod)[1:nN]
	value <- recOpt$ace
  	ectR <- recOpt$se[-1]
  	ectRn <- ectR
  	names(ectRn) <- as.character( (nT+1):(nT+nN) )
  	ValSigmaRec[indL, ((3*nN+1):(4*nN)) ] <- ParadisToNous(ectRn)
	ectRcomp <- c(ectRcomp, ectR)
	#
	arbreTo100 <- compute.brlen(arbre, 100)
	recOptSqu <- newML(EtatFeuilles, arbreTo100)
	ValMuRec[indL, ((4*nN+1):(5*nN)) ] <- ParadisToNousBis(recOptSqu$ace, transf$cod)[1:nN]
	ectR <- rep(0, nN)
	ValSigmaRec[indL, ((4*nN+1):(5*nN)) ] <- ectR
	#
# GLSOU avec theta=z0
#
	funOpt1 <- function(alpha)
	{
		-aceOUv1(phy=arbre, x=EtatFeuilles, alpha, CI=TRUE)$LLik 
	}
	calOp <- optim(par=0.25, fn=funOpt1, method="Brent", lower=0.001, upper=1)
	if (calOp$convergence == 0)
	{
		alphaE <- calOp$par
		recGLSOUv1 <- aceOUv1(phy=arbre, x=EtatFeuilles, alpha=alphaE, CI=T)
		ValMuRec[indL, ((5*nN+1):(6*nN)) ] <- ParadisToNous(recGLSOUv1$ace)
		ectR <- (recGLSOUv1$CI95[,2] - recGLSOUv1$ace) / qnorm(0.975)
		ValSigmaRec[indL, ((5*nN+1):(6*nN)) ] <- ParadisToNous(ectR)
	}
	if (calOp$convergence != 0) 
	{
		ValMuRec[indL, ((5*nN+1):(6*nN)) ] <- rep(NA, nN)
		ectR <- rep(NA, nN)
		ValSigmaRec[indL, ((5*nN+1):(6*nN)) ] <- rep(NA, nN)
	}
	ectRcomp <- c(ectRcomp, ectR)
	#
#GLSOU avec theta non égal à z0
#
	funOpt2 <- function(alpha)
	{
		-aceOUv2(phy=arbre, x=EtatFeuilles, alpha, CI=TRUE)$LLik 
	}
	calOp <- optim(par=0.25, fn=funOpt2, method="Brent", lower=0.001, upper=1)
	if (calOp$convergence == 0)
	{
		alphaE <- calOp$par
		recGLSOUv2 <- aceOUv2(phy=arbre, x=EtatFeuilles, alpha=alphaE, CI=T)
		ValMuRec[indL, ((6*nN+1):(7*nN)) ] <- ParadisToNous(recGLSOUv2$ace)
		ectR <- (recGLSOUv2$CI95[,2] - recGLSOUv2$ace) / qnorm(0.975)
		ValSigmaRec[indL, ((6*nN+1):(7*nN)) ] <- ParadisToNous(ectR)
	}
	if (calOp$convergence != 0) 
	{
		ValMuRec[indL, ((6*nN+1):(7*nN)) ] <- rep(NA, nN)
		ectR <- rep(NA, nN)
		ValSigmaRec[indL, ((6*nN+1):(7*nN)) ] <- rep(NA, nN)
	}
	ectRcomp <- c(ectRcomp, ectR)
	#
# GLSABM
#
	recGLSABM <- aceABM(phy=arbre, x=EtatFeuilles, CI=T)
	ValMuRec[indL, ((7*nN+1):(8*nN)) ] <- ParadisToNous(recGLSABM$ace)
	ectR <- (recGLSABM$CI95[,2] - recGLSABM$ace) / qnorm(0.975)
  	ValSigmaRec[indL, ((7*nN+1):(8*nN)) ] <- ParadisToNous(ectR)
	ectRcomp <- c(ectRcomp, ectR)
#
	if ( sum(is.na(ectRcomp)) == 0 )
	{
		k <- k+1
	}
#
	if ( sum(is.na(ectRcomp)) != 0 )
	{
	nbRec <- nbRec+1
		if (nbRec >= 10) 
		{
		print("A recommencé 10 fois")
		nbRec <- 1
		}
	}	




#nb
# Revoir le calcul de l'erreur, suffixe ABM
# fait formule avec lois NV
# pour l'erreur de la racine, c'est different
# 
# A MODIFIER si la racine n'est pas constante
#
	if (root == "constant")
	{
	muTheo <- rep(c(z0,muCond), nMe)
	sigmaTheo <- c(0,diag(sigmaCond))
	sigmaTheo <- rep(sigmaTheo, nMe)
	}	
#
#
	if (root == "gaussian")
	{
	muTheo <- rep(muCond, nMe)
	sigmaTheo <- diag(sigmaCond)
	sigmaTheo <- rep(sigmaTheo, nMe)
	}
#
# GLSOU, GLSABM
	lois <- c("Stud", "Norm", "Norm", "Stud", "Norm", "Norm", "Norm", "Norm")
	lois <- rep(lois, each=nN)
	#
	ErreursOUnrj[indL, 1:(nMe*nN)] <- NRJ(muR=ValMuRec[indL, (1:(nMe*nN)) ], muT=muTheo, sigma2R=ValSigmaRec[indL, (1:(nMe*nN))]^2, sigma2T=sigmaTheo, loiR=lois, ddl=nN)
	ErreursOUnrj[indL, (nMe*nN+1) ] <- aM[ia]
	ErreursOUnrj[indL, (nMe*nN+2) ] <- bM[ib]
  	ErreursOUnrj[indL, (nMe*nN+3) ] <- sigmaM[j]
	#
	ErreursOUnrjM[indL, 1:(nMe*nN)] <- 2 * sqrt(sigmaTheo) * (sqrt(2) - 1) /sqrt(pi)
  	ErreursOUnrjM[indL, (nMe*nN+1) ] <- aM[ia]
  	ErreursOUnrjM[indL, (nMe*nN+2) ] <- bM[ib]
  	ErreursOUnrjM[indL, (nMe*nN+3) ] <- sigmaM[j]
	#
	ErreursOUks[indL, 1:(nMe*nN)] <- KS(muR=ValMuRec[indL, (1:(nMe*nN)) ], muT=muTheo, sigma2R=ValSigmaRec[indL, (1:(nMe*nN))]^2, sigma2T=sigmaTheo, loiR=lois, ddl=nN)
	ErreursOUks[indL, (nMe*nN+1) ] <- aM[ia]
  	ErreursOUks[indL, (nMe*nN+2) ] <- bM[ib]
  	ErreursOUks[indL, (nMe*nN+3) ] <- sigmaM[j]
	#
	ErreursOUmse[indL, 1:(nMe*nN)] <- sigmaTheo + (muTheo-ValMuRec[indL, (1:(nMe*nN)) ])^2
  	ErreursOUmse[indL, (nMe*nN+1) ] <- aM[ia]
  	ErreursOUmse[indL, (nMe*nN+2) ] <- bM[ib]
  	ErreursOUmse[indL, (nMe*nN+3) ] <- sigmaM[j]
	#
	ErreursOUbiaisG[indL, 1:(nMe*nN)] <- muTheo - ValMuRec[indL, (1:(nMe*nN)) ]
  	ErreursOUbiaisG[indL, (nMe*nN+1) ] <- aM[ia]
  	ErreursOUbiaisG[indL, (nMe*nN+2) ] <- bM[ib]
  	ErreursOUbiaisG[indL, (nMe*nN+3) ] <- sigmaM[j]
	#
	ErreursOUbiaisAbsG[indL, 1:(nMe*nN)] <- NRJ(muR=ValMuRec[indL, (1:(nMe*nN)) ], muT=muTheo, sigma2R=rep(0,(nMe*nN)), sigma2T=sigmaTheo, loiR=rep("Norm", (nMe*nN)), ddl=nN)
  	ErreursOUbiaisAbsG[indL, (nMe*nN+1) ] <- aM[ia]
  	ErreursOUbiaisAbsG[indL, (nMe*nN+2) ] <- bM[ib]
  	ErreursOUbiaisAbsG[indL, (nMe*nN+3) ] <- sigmaM[j]


}
}
}
#
}
#
#
write.csv2(ErreursOUnrj, "NRJmodOU.csv", quote=F, row.names=F)
write.csv2(ErreursOUks, "KSmodOU.csv", quote=F, row.names=F)
write.csv2(ErreursOUmse, "MSEmodOU.csv", quote=F, row.names=F)
write.csv2(ErreursOUbiaisG, "BGmodOU.csv", quote=F, row.names=F)
write.csv2(ErreursOUbiaisAbsG, "BABSGmodOU.csv", quote=F, row.names=F)
