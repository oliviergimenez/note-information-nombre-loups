########################################################################
# fonction pour eviter explosion de log(0)
########################################################################
logprot <- function(v){
  eps <- 2.2204e-016
  u <- log(eps) * (1+vector(length=length(v)))
  index <- (v>eps)
  u[index] <- log(v[index])
  u
}

############################################################
# fonction pour eviter backward sequence
# CC from wrapr package
############################################################
seqi <- function (a, b) 
{
  a <- ceiling(a)
  b <- floor(b)
  if (a > b) {
    return(integer(0))
  }
  seq(a, b, by = 1L)
}

#################################################################################
# pi, phi(a2), p(mix), r(mix), psi 
#################################################################################
devCJS1 <- function(b,data,eff,e,lc,garb,nh,km1){
  # b parametres
  # data est le bloc de CH
  # eff les effectifs (on peut imaginer grouper les individus qui ont des histoires communes)
  # e le vecteur des dates des premières capture
  # lc le vecteur des dates de dernières captures (LECA vs. PASS et ANTA)
  # garb le vecteur des états de départ
  # nh nb individus
  # km1 le nombre d'occasions de capture - 1
  
  #---------------------------------------------
  # lien logit pour les survies (transience)
  phi1 <- 1/(1+exp(-b[1]))
  phi2 <- 1/(1+exp(-b[2]))
  # lien logit pour les captures (hétérogénéité)
  pp1 <- 1/(1+exp(-b[3]))
  pp2 <- 1/(1+exp(-b[4]))
  # lien logit pour les reprises (hétérogénéité)
  lambda1 <- 1/(1+exp(-b[5]))
  lambda2 <- 1/(1+exp(-b[6]))
  # prop
  prop <- 1/(1+exp(-b[7]))
  # trans entre classes d'hétérogénéité
  psi12 <- 1/(1+exp(-b[8]))
  psi21 <- 1/(1+exp(-b[9]))
  #---------------------------------------------
  
  
  #---------------------------------------------
  # Probabilités des événements (lignes) conditionnellement aux états (colonnes)
  B <- matrix(c(1-pp1,1-pp2,1-lambda1,1-lambda2,1,pp1,pp2,0,0,0,0,0,lambda1,lambda2,0),
              nrow=3,ncol=5,byrow=T)
  
  # Output initiaux
  BE <- matrix(c(0,0,1,1,1,1,1,0,0,0,0,0,0,0,0),nrow=3,ncol=5,byrow=T) 
  
  # Matrice stochastique de transitions entre états (processus Markovien)
  # age 1
  PHI1 <- matrix(c(phi1,0,1-phi1,0,0,
                   0,phi1,0,1-phi1,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  # age 2
  PHI2 <- matrix(c(phi2,0,1-phi2,0,0,
                   0,phi2,0,1-phi2,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  PSI <- matrix(c((1-psi12),psi12,0,0,0,
                  psi21,(1-psi21),0,0,0,
                  0,0,1,0,0,
                  0,0,0,1,0,
                  0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  A1 <- PHI1 %*% PSI
  A2 <- PHI2 %*% PSI
  
  # Distribution des états initiaux
  PI <- c(prop,1-prop,0,0,0)
  #---------------------------------------------
  
  #---------------------------------------------
  # on forme la vraisemblance, ou plutôt la log(vraisemblance) 
  # on linéarise car c'est plus stable numériquement de maximiser une somme qu'un produit
  # en fait, on forme -log(vraisemblance) plutôt, car R permet seulement de minimiser des fonctions
  # voir Pradel 2005, 2009
  
  l <- 0
  for (i in 1:nh) # boucles sur histoires de capture
  {
    ei <- e[i] # date marquage
    oe <- garb[i] + 1 # événement initial
    evennt <- data[,i] + 1 # les non-détections deviennent des 1, et les détections des 2
    ALPHA <- PI*BE[oe,]
    for (j in seqi(ei+1,lc[i])) # on conditionne par rapport à la première capture
    {
      if (j==(ei+1)) ALPHA <- (ALPHA %*% A1)*B[evennt[j],]
      if (j>(ei+1)) ALPHA <- (ALPHA %*% A2)*B[evennt[j],]
    }
    l <- l + logprot(sum(ALPHA))*eff[i]
  }
  l <- -l
  l
  #---------------------------------------------
  
} # fin fn vrais




#################################################################################
# pi, phi(a2), p, r(mix), psi ##########################
#################################################################################
devCJS2 <- function(b,data,eff,e,lc,garb,nh,km1){
  
  #---------------------------------------------
  # lien logit pour les survies (transience)
  phi1 <- 1/(1+exp(-b[1]))
  phi2 <- 1/(1+exp(-b[2]))
  # lien logit pour les captures (homogénéité)
  pp <- 1/(1+exp(-b[3]))
  # lien logit pour les reprises (hétérogénéité)
  lambda1 <- 1/(1+exp(-b[4]))
  lambda2 <- 1/(1+exp(-b[5]))
  # prop
  prop <- 1/(1+exp(-b[6]))
  # trans entre classes d'hétérogénéité
  psi12 <- 1/(1+exp(-b[7]))
  psi21 <- 1/(1+exp(-b[8]))
  #---------------------------------------------
  
  
  #---------------------------------------------
  # Probabilités des événements (lignes) conditionnellement aux états (colonnes)
  B <- matrix(c(1-pp,1-pp,1-lambda1,1-lambda2,1,pp,pp,0,0,0,0,0,lambda1,lambda2,0),
              nrow=3,ncol=5,byrow=T)
  
  # Output initiaux
  BE <- matrix(c(0,0,1,1,1,1,1,0,0,0,0,0,0,0,0),nrow=3,ncol=5,byrow=T) 
  
  # Matrice stochastique de transitions entre états (processus Markovien)
  # age 1
  PHI1 <- matrix(c(phi1,0,1-phi1,0,0,
                   0,phi1,0,1-phi1,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  # age 2
  PHI2 <- matrix(c(phi2,0,1-phi2,0,0,
                   0,phi2,0,1-phi2,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  PSI <- matrix(c((1-psi12),psi12,0,0,0,
                  psi21,(1-psi21),0,0,0,
                  0,0,1,0,0,
                  0,0,0,1,0,
                  0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  A1 <- PHI1 %*% PSI
  A2 <- PHI2 %*% PSI
  
  # Distribution des états initiaux
  PI <- c(prop,1-prop,0,0,0)
  #---------------------------------------------
  
  #---------------------------------------------
  # on forme la vraisemblance, ou plutôt la log(vraisemblance) 
  # on linéarise car c'est plus stable numériquement de maximiser une somme qu'un produit
  # en fait, on forme -log(vraisemblance) plutôt, car R permet seulement de minimiser des fonctions
  # voir Pradel 2005, 2009
  
  l <- 0
  for (i in 1:nh) # boucles sur histoires de capture
  {
    ei <- e[i] # date marquage
    oe <- garb[i] + 1 # événement initial
    evennt <- data[,i] + 1 # les non-détections deviennent des 1, et les détections des 2
    ALPHA <- PI*BE[oe,]
    for (j in seqi(ei+1,lc[i])) # on conditionne par rapport à la première capture
    {
      if (j==(ei+1)) ALPHA <- (ALPHA %*% A1)*B[evennt[j],]
      if (j>(ei+1)) ALPHA <- (ALPHA %*% A2)*B[evennt[j],]
    }
    l <- l + logprot(sum(ALPHA))*eff[i]
  }
  l <- -l
  l
  #---------------------------------------------
  
} # fin fn vrais


#################################################################################
# pi, phi(a2), p(mix), r, psi ##########################
#################################################################################
devCJS3 <- function(b,data,eff,e,lc, garb,nh,km1){
  
  #---------------------------------------------
  # lien logit pour les survies (transience)
  phi1 <- 1/(1+exp(-b[1]))
  phi2 <- 1/(1+exp(-b[2]))
  # lien logit pour les captures (hétérogénéité)
  pp1 <- 1/(1+exp(-b[3]))
  pp2 <- 1/(1+exp(-b[4]))
  # lien logit pour les reprises (hétérogénéité)
  lambda <- 1/(1+exp(-b[5]))
  # prop
  prop <- 1/(1+exp(-b[6]))
  # trans entre classes d'hétérogénéité
  psi12 <- 1/(1+exp(-b[7]))
  psi21 <- 1/(1+exp(-b[8]))
  #---------------------------------------------
  
  
  #---------------------------------------------
  # Probabilités des événements (lignes) conditionnellement aux états (colonnes)
  B <- matrix(c(1-pp1,1-pp2,1-lambda,1-lambda,1,pp1,pp2,0,0,0,0,0,lambda,lambda,0),
              nrow=3,ncol=5,byrow=T)
  
  # Output initiaux
  BE <- matrix(c(0,0,1,1,1,1,1,0,0,0,0,0,0,0,0),nrow=3,ncol=5,byrow=T) 
  
  # Matrice stochastique de transitions entre états (processus Markovien)
  # age 1
  PHI1 <- matrix(c(phi1,0,1-phi1,0,0,
                   0,phi1,0,1-phi1,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  # age 2
  PHI2 <- matrix(c(phi2,0,1-phi2,0,0,
                   0,phi2,0,1-phi2,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  PSI <- matrix(c((1-psi12),psi12,0,0,0,
                  psi21,(1-psi21),0,0,0,
                  0,0,1,0,0,
                  0,0,0,1,0,
                  0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  A1 <- PHI1 %*% PSI
  A2 <- PHI2 %*% PSI
  
  # Distribution des états initiaux
  PI <- c(prop,1-prop,0,0,0)
  #---------------------------------------------
  
  #---------------------------------------------
  # on forme la vraisemblance, ou plutôt la log(vraisemblance) 
  # on linéarise car c'est plus stable numériquement de maximiser une somme qu'un produit
  # en fait, on forme -log(vraisemblance) plutôt, car R permet seulement de minimiser des fonctions
  # voir Pradel 2005, 2009
  
  l <- 0
  for (i in 1:nh) # boucles sur histoires de capture
  {
    ei <- e[i] # date marquage
    oe <- garb[i] + 1 # événement initial
    evennt <- data[,i] + 1 # les non-détections deviennent des 1, et les détections des 2
    ALPHA <- PI*BE[oe,]
    for (j in seqi(ei+1,lc[i])) # on conditionne par rapport à la première capture
    {
      if (j==(ei+1)) ALPHA <- (ALPHA %*% A1)*B[evennt[j],]
      if (j>(ei+1)) ALPHA <- (ALPHA %*% A2)*B[evennt[j],]
    }
    l <- l + logprot(sum(ALPHA))*eff[i]
  }
  l <- -l
  l
  #---------------------------------------------
  
} # fin fn vrais



#################################################################################
# pi, phi(a2), p, r, psi ##########################
#################################################################################
devCJS4 <- function(b,data,eff,e,lc,garb,nh,km1){
  
  #---------------------------------------------
  # lien logit pour les survies (transience)
  phi1 <- 1/(1+exp(-b[1]))
  phi2 <- 1/(1+exp(-b[2]))
  # lien logit pour les captures (hétérogénéité)
  pp <- 1/(1+exp(-b[3]))
  # lien logit pour les reprises (hétérogénéité)
  lambda <- 1/(1+exp(-b[4]))
  # prop
  prop <- 1/(1+exp(-b[5]))
  # trans entre classes d'hétérogénéité
  psi12 <- 1/(1+exp(-b[6]))
  psi21 <- 1/(1+exp(-b[7]))
  #---------------------------------------------
  
  
  #---------------------------------------------
  # Probabilités des événements (lignes) conditionnellement aux états (colonnes)
  B <- matrix(c(1-pp,1-pp,1-lambda,1-lambda,1,pp,pp,0,0,0,0,0,lambda,lambda,0),
              nrow=3,ncol=5,byrow=T)
  
  # Output initiaux
  BE <- matrix(c(0,0,1,1,1,1,1,0,0,0,0,0,0,0,0),nrow=3,ncol=5,byrow=T) 
  
  # Matrice stochastique de transitions entre états (processus Markovien)
  # age 1
  PHI1 <- matrix(c(phi1,0,1-phi1,0,0,
                   0,phi1,0,1-phi1,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  # age 2
  PHI2 <- matrix(c(phi2,0,1-phi2,0,0,
                   0,phi2,0,1-phi2,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  PSI <- matrix(c((1-psi12),psi12,0,0,0,
                  psi21,(1-psi21),0,0,0,
                  0,0,1,0,0,
                  0,0,0,1,0,
                  0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  A1 <- PHI1 %*% PSI
  A2 <- PHI2 %*% PSI
  
  # Distribution des états initiaux
  PI <- c(prop,1-prop,0,0,0)
  #---------------------------------------------
  
  #---------------------------------------------
  # on forme la vraisemblance, ou plutôt la log(vraisemblance) 
  # on linéarise car c'est plus stable numériquement de maximiser une somme qu'un produit
  # en fait, on forme -log(vraisemblance) plutôt, car R permet seulement de minimiser des fonctions
  # voir Pradel 2005, 2009
  
  l <- 0
  for (i in 1:nh) # boucles sur histoires de capture
  {
    ei <- e[i] # date marquage
    oe <- garb[i] + 1 # événement initial
    evennt <- data[,i] + 1 # les non-détections deviennent des 1, et les détections des 2
    ALPHA <- PI*BE[oe,]
    for (j in seqi(ei+1,lc[i])) # on conditionne par rapport à la première capture
    {
      if (j==(ei+1)) ALPHA <- (ALPHA %*% A1)*B[evennt[j],]
      if (j>(ei+1)) ALPHA <- (ALPHA %*% A2)*B[evennt[j],]
    }
    l <- l + logprot(sum(ALPHA))*eff[i]
  }
  l <- -l
  l
  #---------------------------------------------
  
} # fin fn vrais


#################################################################################
# pi, phi(a2*mix), p(mix), r(mix), psi ##########################
#################################################################################
devCJS5 <- function(b,data,eff,e,lc,garb,nh,km1){
  
  #---------------------------------------------
  # lien logit pour les survies (transience)
  phi1c1 <- 1/(1+exp(-b[1]))
  phi2c1 <- 1/(1+exp(-b[2]))
  phi1c2 <- 1/(1+exp(-b[3]))
  phi2c2 <- 1/(1+exp(-b[4]))
  # lien logit pour les captures (hétérogénéité)
  pp1 <- 1/(1+exp(-b[5]))
  pp2 <- 1/(1+exp(-b[6]))
  # lien logit pour les reprises (hétérogénéité)
  lambda1 <- 1/(1+exp(-b[7]))
  lambda2 <- 1/(1+exp(-b[8]))
  # prop
  prop <- 1/(1+exp(-b[9]))
  # trans entre classes d'hétérogénéité
  psi12 <- 1/(1+exp(-b[10]))
  psi21 <- 1/(1+exp(-b[11]))
  #---------------------------------------------
  
  
  #---------------------------------------------
  # Probabilités des événements (lignes) conditionnellement aux états (colonnes)
  B <- matrix(c(1-pp1,1-pp2,1-lambda1,1-lambda2,1,pp1,pp2,0,0,0,0,0,lambda1,lambda2,0),
              nrow=3,ncol=5,byrow=T)
  
  # Output initiaux
  BE <- matrix(c(0,0,1,1,1,1,1,0,0,0,0,0,0,0,0),nrow=3,ncol=5,byrow=T) 
  
  # Matrice stochastique de transitions entre états (processus Markovien)
  # age 1
  PHI1 <- matrix(c(phi1c1,0,1-phi1c1,0,0,
                   0,phi1c2,0,1-phi1c2,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  # age 2
  PHI2 <- matrix(c(phi2c1,0,1-phi2c1,0,0,
                   0,phi2c2,0,1-phi2c2,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  PSI <- matrix(c((1-psi12),psi12,0,0,0,
                  psi21,(1-psi21),0,0,0,
                  0,0,1,0,0,
                  0,0,0,1,0,
                  0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  A1 <- PHI1 %*% PSI
  A2 <- PHI2 %*% PSI
  
  # Distribution des états initiaux
  PI <- c(prop,1-prop,0,0,0)
  #---------------------------------------------
  
  #---------------------------------------------
  # on forme la vraisemblance, ou plutôt la log(vraisemblance) 
  # on linéarise car c'est plus stable numériquement de maximiser une somme qu'un produit
  # en fait, on forme -log(vraisemblance) plutôt, car R permet seulement de minimiser des fonctions
  # voir Pradel 2005, 2009
  
  l <- 0
  for (i in 1:nh) # boucles sur histoires de capture
  {
    ei <- e[i] # date marquage
    oe <- garb[i] + 1 # événement initial
    evennt <- data[,i] + 1 # les non-détections deviennent des 1, et les détections des 2
    ALPHA <- PI*BE[oe,]
    for (j in seqi(ei+1,lc[i])) # on conditionne par rapport à la première capture
    {
      if (j==(ei+1)) ALPHA <- (ALPHA %*% A1)*B[evennt[j],]
      if (j>(ei+1)) ALPHA <- (ALPHA %*% A2)*B[evennt[j],]
    }
    l <- l + logprot(sum(ALPHA))*eff[i]
  }
  l <- -l
  l
  #---------------------------------------------
  
} # fin fn vrais




#################################################################################
# pi, phi(a2*mix), p, r(mix), psi ##########################
#################################################################################
devCJS6 <- function(b,data,eff,e,lc,garb,nh,km1){
  
  #---------------------------------------------
  # lien logit pour les survies (transience)
  phi1c1 <- 1/(1+exp(-b[1]))
  phi2c1 <- 1/(1+exp(-b[2]))
  phi1c2 <- 1/(1+exp(-b[3]))
  phi2c2 <- 1/(1+exp(-b[4]))
  # lien logit pour les captures (homogénéité)
  pp <- 1/(1+exp(-b[5]))
  # lien logit pour les reprises (hétérogénéité)
  lambda1 <- 1/(1+exp(-b[6]))
  lambda2 <- 1/(1+exp(-b[7]))
  # prop
  prop <- 1/(1+exp(-b[8]))
  # trans entre classes d'hétérogénéité
  psi12 <- 1/(1+exp(-b[9]))
  psi21 <- 1/(1+exp(-b[10]))
  #---------------------------------------------
  
  
  #---------------------------------------------
  # Probabilités des événements (lignes) conditionnellement aux états (colonnes)
  B <- matrix(c(1-pp,1-pp,1-lambda1,1-lambda2,1,pp,pp,0,0,0,0,0,lambda1,lambda2,0),
              nrow=3,ncol=5,byrow=T)
  
  # Output initiaux
  BE <- matrix(c(0,0,1,1,1,1,1,0,0,0,0,0,0,0,0),nrow=3,ncol=5,byrow=T) 
  
  # Matrice stochastique de transitions entre états (processus Markovien)
  # age 1
  PHI1 <- matrix(c(phi1c1,0,1-phi1c1,0,0,
                   0,phi1c2,0,1-phi1c2,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  # age 2
  PHI2 <- matrix(c(phi2c1,0,1-phi2c1,0,0,
                   0,phi2c2,0,1-phi2c2,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  PSI <- matrix(c((1-psi12),psi12,0,0,0,
                  psi21,(1-psi21),0,0,0,
                  0,0,1,0,0,
                  0,0,0,1,0,
                  0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  A1 <- PHI1 %*% PSI
  A2 <- PHI2 %*% PSI
  
  # Distribution des états initiaux
  PI <- c(prop,1-prop,0,0,0)
  #---------------------------------------------
  
  #---------------------------------------------
  # on forme la vraisemblance, ou plutôt la log(vraisemblance) 
  # on linéarise car c'est plus stable numériquement de maximiser une somme qu'un produit
  # en fait, on forme -log(vraisemblance) plutôt, car R permet seulement de minimiser des fonctions
  # voir Pradel 2005, 2009
  
  l <- 0
  for (i in 1:nh) # boucles sur histoires de capture
  {
    ei <- e[i] # date marquage
    oe <- garb[i] + 1 # événement initial
    evennt <- data[,i] + 1 # les non-détections deviennent des 1, et les détections des 2
    ALPHA <- PI*BE[oe,]
    for (j in seqi(ei+1,lc[i])) # on conditionne par rapport à la première capture
    {
      if (j==(ei+1)) ALPHA <- (ALPHA %*% A1)*B[evennt[j],]
      if (j>(ei+1)) ALPHA <- (ALPHA %*% A2)*B[evennt[j],]
    }
    l <- l + logprot(sum(ALPHA))*eff[i]
  }
  l <- -l
  l
  #---------------------------------------------
  
} # fin fn vrais


#################################################################################
# pi, phi(a2*mix), p(mix), r, psi ##########################
#################################################################################
devCJS7 <- function(b,data,eff,e,lc,garb,nh,km1){
  
  #---------------------------------------------
  # lien logit pour les survies (transience)
  phi1c1 <- 1/(1+exp(-b[1]))
  phi2c1 <- 1/(1+exp(-b[2]))
  phi1c2 <- 1/(1+exp(-b[3]))
  phi2c2 <- 1/(1+exp(-b[4]))
  # lien logit pour les captures (hétérogénéité)
  pp1 <- 1/(1+exp(-b[5]))
  pp2 <- 1/(1+exp(-b[6]))
  # lien logit pour les reprises (hétérogénéité)
  lambda <- 1/(1+exp(-b[7]))
  # prop
  prop <- 1/(1+exp(-b[8]))
  # trans entre classes d'hétérogénéité
  psi12 <- 1/(1+exp(-b[9]))
  psi21 <- 1/(1+exp(-b[10]))
  #---------------------------------------------
  
  
  #---------------------------------------------
  # Probabilités des événements (lignes) conditionnellement aux états (colonnes)
  B <- matrix(c(1-pp1,1-pp2,1-lambda,1-lambda,1,pp1,pp2,0,0,0,0,0,lambda,lambda,0),
              nrow=3,ncol=5,byrow=T)
  
  # Output initiaux
  BE <- matrix(c(0,0,1,1,1,1,1,0,0,0,0,0,0,0,0),nrow=3,ncol=5,byrow=T) 
  
  # Matrice stochastique de transitions entre états (processus Markovien)
  # age 1
  PHI1 <- matrix(c(phi1c1,0,1-phi1c1,0,0,
                   0,phi1c2,0,1-phi1c2,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  # age 2
  PHI2 <- matrix(c(phi2c1,0,1-phi2c1,0,0,
                   0,phi2c2,0,1-phi2c2,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  PSI <- matrix(c((1-psi12),psi12,0,0,0,
                  psi21,(1-psi21),0,0,0,
                  0,0,1,0,0,
                  0,0,0,1,0,
                  0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  A1 <- PHI1 %*% PSI
  A2 <- PHI2 %*% PSI
  
  # Distribution des états initiaux
  PI <- c(prop,1-prop,0,0,0)
  #---------------------------------------------
  
  #---------------------------------------------
  # on forme la vraisemblance, ou plutôt la log(vraisemblance) 
  # on linéarise car c'est plus stable numériquement de maximiser une somme qu'un produit
  # en fait, on forme -log(vraisemblance) plutôt, car R permet seulement de minimiser des fonctions
  # voir Pradel 2005, 2009
  
  l <- 0
  for (i in 1:nh) # boucles sur histoires de capture
  {
    ei <- e[i] # date marquage
    oe <- garb[i] + 1 # événement initial
    evennt <- data[,i] + 1 # les non-détections deviennent des 1, et les détections des 2
    ALPHA <- PI*BE[oe,]
    for (j in seqi(ei+1,lc[i])) # on conditionne par rapport à la première capture
    {
      if (j==(ei+1)) ALPHA <- (ALPHA %*% A1)*B[evennt[j],]
      if (j>(ei+1)) ALPHA <- (ALPHA %*% A2)*B[evennt[j],]
    }
    l <- l + logprot(sum(ALPHA))*eff[i]
  }
  l <- -l
  l
  #---------------------------------------------
  
} # fin fn vrais



#################################################################################
# pi, phi(a2*mix), p, r, psi ##########################
#################################################################################
devCJS8 <- function(b,data,eff,e,lc,garb,nh,km1){
  
  #---------------------------------------------
  # lien logit pour les survies (transience)
  phi1c1 <- 1/(1+exp(-b[1]))
  phi2c1 <- 1/(1+exp(-b[2]))
  phi1c2 <- 1/(1+exp(-b[3]))
  phi2c2 <- 1/(1+exp(-b[4]))
  # lien logit pour les captures (hétérogénéité)
  pp <- 1/(1+exp(-b[5]))
  # lien logit pour les reprises (hétérogénéité)
  lambda <- 1/(1+exp(-b[6]))
  # prop
  prop <- 1/(1+exp(-b[7]))
  # trans entre classes d'hétérogénéité
  psi12 <- 1/(1+exp(-b[8]))
  psi21 <- 1/(1+exp(-b[9]))
  #---------------------------------------------
  
  
  #---------------------------------------------
  # Probabilités des événements (lignes) conditionnellement aux états (colonnes)
  B <- matrix(c(1-pp,1-pp,1-lambda,1-lambda,1,pp,pp,0,0,0,0,0,lambda,lambda,0),
              nrow=3,ncol=5,byrow=T)
  
  # Output initiaux
  BE <- matrix(c(0,0,1,1,1,1,1,0,0,0,0,0,0,0,0),nrow=3,ncol=5,byrow=T) 
  
  # Matrice stochastique de transitions entre états (processus Markovien)
  # age 1
  PHI1 <- matrix(c(phi1c1,0,1-phi1c1,0,0,
                   0,phi1c2,0,1-phi1c2,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  # age 2
  PHI2 <- matrix(c(phi2c1,0,1-phi2c1,0,0,
                   0,phi2c2,0,1-phi2c2,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  PSI <- matrix(c((1-psi12),psi12,0,0,0,
                  psi21,(1-psi21),0,0,0,
                  0,0,1,0,0,
                  0,0,0,1,0,
                  0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  A1 <- PHI1 %*% PSI
  A2 <- PHI2 %*% PSI
  
  # Distribution des états initiaux
  PI <- c(prop,1-prop,0,0,0)
  #---------------------------------------------
  
  #---------------------------------------------
  # on forme la vraisemblance, ou plutôt la log(vraisemblance) 
  # on linéarise car c'est plus stable numériquement de maximiser une somme qu'un produit
  # en fait, on forme -log(vraisemblance) plutôt, car R permet seulement de minimiser des fonctions
  # voir Pradel 2005, 2009
  
  l <- 0
  for (i in 1:nh) # boucles sur histoires de capture
  {
    ei <- e[i] # date marquage
    oe <- garb[i] + 1 # événement initial
    evennt <- data[,i] + 1 # les non-détections deviennent des 1, et les détections des 2
    ALPHA <- PI*BE[oe,]
    for (j in seqi(ei+1,lc[i])) # on conditionne par rapport à la première capture
    {
      if (j==(ei+1)) ALPHA <- (ALPHA %*% A1)*B[evennt[j],]
      if (j>(ei+1)) ALPHA <- (ALPHA %*% A2)*B[evennt[j],]
    }
    l <- l + logprot(sum(ALPHA))*eff[i]
  }
  l <- -l
  l
  #---------------------------------------------
  
} # fin fn vrais



#################################################################################
# pi, phi(a2), p(mix), r(mix) ##########################
#################################################################################
devCJS9 <- function(b,data,eff,e,lc,garb,nh,km1){
  
  #---------------------------------------------
  # lien logit pour les survies (transience)
  phi1 <- 1/(1+exp(-b[1]))
  phi2 <- 1/(1+exp(-b[2]))
  # lien logit pour les captures (hétérogénéité)
  pp1 <- 1/(1+exp(-b[3]))
  pp2 <- 1/(1+exp(-b[4]))
  # lien logit pour les reprises (hétérogénéité)
  lambda1 <- 1/(1+exp(-b[5]))
  lambda2 <- 1/(1+exp(-b[6]))
  # prop
  prop <- 1/(1+exp(-b[7]))
  #---------------------------------------------
  
  
  #---------------------------------------------
  # Probabilités des événements (lignes) conditionnellement aux états (colonnes)
  B <- matrix(c(1-pp1,1-pp2,1-lambda1,1-lambda2,1,pp1,pp2,0,0,0,0,0,lambda1,lambda2,0),
              nrow=3,ncol=5,byrow=T)
  
  # Output initiaux
  BE <- matrix(c(0,0,1,1,1,1,1,0,0,0,0,0,0,0,0),nrow=3,ncol=5,byrow=T) 
  
  # Matrice stochastique de transitions entre états (processus Markovien)
  # age 1
  PHI1 <- matrix(c(phi1,0,1-phi1,0,0,
                   0,phi1,0,1-phi1,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  # age 2
  PHI2 <- matrix(c(phi2,0,1-phi2,0,0,
                   0,phi2,0,1-phi2,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  PSI <- matrix(c(1,0,0,0,0,
                  0,1,0,0,0,
                  0,0,1,0,0,
                  0,0,0,1,0,
                  0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  A1 <- PHI1 %*% PSI
  A2 <- PHI2 %*% PSI
  
  # Distribution des états initiaux
  PI <- c(prop,1-prop,0,0,0)
  #---------------------------------------------
  
  #---------------------------------------------
  # on forme la vraisemblance, ou plutôt la log(vraisemblance) 
  # on linéarise car c'est plus stable numériquement de maximiser une somme qu'un produit
  # en fait, on forme -log(vraisemblance) plutôt, car R permet seulement de minimiser des fonctions
  # voir Pradel 2005, 2009
  
  l <- 0
  for (i in 1:nh) # boucles sur histoires de capture
  {
    ei <- e[i] # date marquage
    oe <- garb[i] + 1 # événement initial
    evennt <- data[,i] + 1 # les non-détections deviennent des 1, et les détections des 2
    ALPHA <- PI*BE[oe,]
    for (j in seqi(ei+1,lc[i])) # on conditionne par rapport à la première capture
    {
      if (j==(ei+1)) ALPHA <- (ALPHA %*% A1)*B[evennt[j],]
      if (j>(ei+1)) ALPHA <- (ALPHA %*% A2)*B[evennt[j],]
    }
    l <- l + logprot(sum(ALPHA))*eff[i]
  }
  l <- -l
  l
  #---------------------------------------------
  
} # fin fn vrais




#################################################################################
# pi, phi(a2), p, r(mix) ##########################
#################################################################################
devCJS10 <- function(b,data,eff,e,lc,garb,nh,km1){
  
  #---------------------------------------------
  # lien logit pour les survies (transience)
  phi1 <- 1/(1+exp(-b[1]))
  phi2 <- 1/(1+exp(-b[2]))
  # lien logit pour les captures (homogénéité)
  pp <- 1/(1+exp(-b[3]))
  # lien logit pour les reprises (hétérogénéité)
  lambda1 <- 1/(1+exp(-b[4]))
  lambda2 <- 1/(1+exp(-b[5]))
  # prop
  prop <- 1/(1+exp(-b[6]))
  #---------------------------------------------
  
  
  #---------------------------------------------
  # Probabilités des événements (lignes) conditionnellement aux états (colonnes)
  B <- matrix(c(1-pp,1-pp,1-lambda1,1-lambda2,1,pp,pp,0,0,0,0,0,lambda1,lambda2,0),
              nrow=3,ncol=5,byrow=T)
  
  # Output initiaux
  BE <- matrix(c(0,0,1,1,1,1,1,0,0,0,0,0,0,0,0),nrow=3,ncol=5,byrow=T) 
  
  # Matrice stochastique de transitions entre états (processus Markovien)
  # age 1
  PHI1 <- matrix(c(phi1,0,1-phi1,0,0,
                   0,phi1,0,1-phi1,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  # age 2
  PHI2 <- matrix(c(phi2,0,1-phi2,0,0,
                   0,phi2,0,1-phi2,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  PSI <- matrix(c(1,0,0,0,0,
                  0,1,0,0,0,
                  0,0,1,0,0,
                  0,0,0,1,0,
                  0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  A1 <- PHI1 %*% PSI
  A2 <- PHI2 %*% PSI
  
  # Distribution des états initiaux
  PI <- c(prop,1-prop,0,0,0)
  #---------------------------------------------
  
  #---------------------------------------------
  # on forme la vraisemblance, ou plutôt la log(vraisemblance) 
  # on linéarise car c'est plus stable numériquement de maximiser une somme qu'un produit
  # en fait, on forme -log(vraisemblance) plutôt, car R permet seulement de minimiser des fonctions
  # voir Pradel 2005, 2009
  
  l <- 0
  for (i in 1:nh) # boucles sur histoires de capture
  {
    ei <- e[i] # date marquage
    oe <- garb[i] + 1 # événement initial
    evennt <- data[,i] + 1 # les non-détections deviennent des 1, et les détections des 2
    ALPHA <- PI*BE[oe,]
    for (j in seqi(ei+1,lc[i])) # on conditionne par rapport à la première capture
    {
      if (j==(ei+1)) ALPHA <- (ALPHA %*% A1)*B[evennt[j],]
      if (j>(ei+1)) ALPHA <- (ALPHA %*% A2)*B[evennt[j],]
    }
    l <- l + logprot(sum(ALPHA))*eff[i]
  }
  l <- -l
  l
  #---------------------------------------------
  
} # fin fn vrais


#################################################################################
# pi, phi(a2), p(mix), r ##########################
#################################################################################
devCJS11 <- function(b,data,eff,e,lc,garb,nh,km1){
  
  #---------------------------------------------
  # lien logit pour les survies (transience)
  phi1 <- 1/(1+exp(-b[1]))
  phi2 <- 1/(1+exp(-b[2]))
  # lien logit pour les captures (hétérogénéité)
  pp1 <- 1/(1+exp(-b[3]))
  pp2 <- 1/(1+exp(-b[4]))
  # lien logit pour les reprises (hétérogénéité)
  lambda <- 1/(1+exp(-b[5]))
  # prop
  prop <- 1/(1+exp(-b[6]))
  #---------------------------------------------
  
  
  #---------------------------------------------
  # Probabilités des événements (lignes) conditionnellement aux états (colonnes)
  B <- matrix(c(1-pp1,1-pp2,1-lambda,1-lambda,1,pp1,pp2,0,0,0,0,0,lambda,lambda,0),
              nrow=3,ncol=5,byrow=T)
  
  # Output initiaux
  BE <- matrix(c(0,0,1,1,1,1,1,0,0,0,0,0,0,0,0),nrow=3,ncol=5,byrow=T) 
  
  # Matrice stochastique de transitions entre états (processus Markovien)
  # age 1
  PHI1 <- matrix(c(phi1,0,1-phi1,0,0,
                   0,phi1,0,1-phi1,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  # age 2
  PHI2 <- matrix(c(phi2,0,1-phi2,0,0,
                   0,phi2,0,1-phi2,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  PSI <- matrix(c(1,0,0,0,0,
                  0,1,0,0,0,
                  0,0,1,0,0,
                  0,0,0,1,0,
                  0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  A1 <- PHI1 %*% PSI
  A2 <- PHI2 %*% PSI
  
  # Distribution des états initiaux
  PI <- c(prop,1-prop,0,0,0)
  #---------------------------------------------
  
  #---------------------------------------------
  # on forme la vraisemblance, ou plutôt la log(vraisemblance) 
  # on linéarise car c'est plus stable numériquement de maximiser une somme qu'un produit
  # en fait, on forme -log(vraisemblance) plutôt, car R permet seulement de minimiser des fonctions
  # voir Pradel 2005, 2009
  
  l <- 0
  for (i in 1:nh) # boucles sur histoires de capture
  {
    ei <- e[i] # date marquage
    oe <- garb[i] + 1 # événement initial
    evennt <- data[,i] + 1 # les non-détections deviennent des 1, et les détections des 2
    ALPHA <- PI*BE[oe,]
    for (j in seqi(ei+1,lc[i])) # on conditionne par rapport à la première capture
    {
      if (j==(ei+1)) ALPHA <- (ALPHA %*% A1)*B[evennt[j],]
      if (j>(ei+1)) ALPHA <- (ALPHA %*% A2)*B[evennt[j],]
    }
    l <- l + logprot(sum(ALPHA))*eff[i]
  }
  l <- -l
  l
  #---------------------------------------------
  
} # fin fn vrais



#################################################################################
# pi, phi(a2), p, r ##########################
#################################################################################
devCJS12 <- function(b,data,eff,e,lc,garb,nh,km1){
  
  #---------------------------------------------
  # lien logit pour les survies (transience)
  phi1 <- 1/(1+exp(-b[1]))
  phi2 <- 1/(1+exp(-b[2]))
  # lien logit pour les captures (hétérogénéité)
  pp <- 1/(1+exp(-b[3]))
  # lien logit pour les reprises (hétérogénéité)
  lambda <- 1/(1+exp(-b[4]))
  # prop
  prop <- 1/(1+exp(-b[5]))
  #---------------------------------------------
  
  
  #---------------------------------------------
  # Probabilités des événements (lignes) conditionnellement aux états (colonnes)
  B <- matrix(c(1-pp,1-pp,1-lambda,1-lambda,1,pp,pp,0,0,0,0,0,lambda,lambda,0),
              nrow=3,ncol=5,byrow=T)
  
  # Output initiaux
  BE <- matrix(c(0,0,1,1,1,1,1,0,0,0,0,0,0,0,0),nrow=3,ncol=5,byrow=T) 
  
  # Matrice stochastique de transitions entre états (processus Markovien)
  # age 1
  PHI1 <- matrix(c(phi1,0,1-phi1,0,0,
                   0,phi1,0,1-phi1,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  # age 2
  PHI2 <- matrix(c(phi2,0,1-phi2,0,0,
                   0,phi2,0,1-phi2,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  PSI <- matrix(c(1,0,0,0,0,
                  0,1,0,0,0,
                  0,0,1,0,0,
                  0,0,0,1,0,
                  0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  A1 <- PHI1 %*% PSI
  A2 <- PHI2 %*% PSI
  
  # Distribution des états initiaux
  PI <- c(prop,1-prop,0,0,0)
  #---------------------------------------------
  
  #---------------------------------------------
  # on forme la vraisemblance, ou plutôt la log(vraisemblance) 
  # on linéarise car c'est plus stable numériquement de maximiser une somme qu'un produit
  # en fait, on forme -log(vraisemblance) plutôt, car R permet seulement de minimiser des fonctions
  # voir Pradel 2005, 2009
  
  l <- 0
  for (i in 1:nh) # boucles sur histoires de capture
  {
    ei <- e[i] # date marquage
    oe <- garb[i] + 1 # événement initial
    evennt <- data[,i] + 1 # les non-détections deviennent des 1, et les détections des 2
    ALPHA <- PI*BE[oe,]
    for (j in seqi(ei+1,lc[i])) # on conditionne par rapport à la première capture
    {
      if (j==(ei+1)) ALPHA <- (ALPHA %*% A1)*B[evennt[j],]
      if (j>(ei+1)) ALPHA <- (ALPHA %*% A2)*B[evennt[j],]
    }
    l <- l + logprot(sum(ALPHA))*eff[i]
  }
  l <- -l
  l
  #---------------------------------------------
  
} # fin fn vrais



#################################################################################
# pi, phi(a2*mix), p(mix), r(mix) ##########################
#################################################################################
devCJS13 <- function(b,data,eff,e,lc,garb,nh,km1){
  
  #---------------------------------------------
  # lien logit pour les survies (transience)
  phi1c1 <- 1/(1+exp(-b[1]))
  phi2c1 <- 1/(1+exp(-b[2]))
  phi1c2 <- 1/(1+exp(-b[3]))
  phi2c2 <- 1/(1+exp(-b[4]))
  # lien logit pour les captures (hétérogénéité)
  pp1 <- 1/(1+exp(-b[5]))
  pp2 <- 1/(1+exp(-b[6]))
  # lien logit pour les reprises (hétérogénéité)
  lambda1 <- 1/(1+exp(-b[7]))
  lambda2 <- 1/(1+exp(-b[8]))
  # prop
  prop <- 1/(1+exp(-b[9]))
  #---------------------------------------------
  
  
  #---------------------------------------------
  # Probabilités des événements (lignes) conditionnellement aux états (colonnes)
  B <- matrix(c(1-pp1,1-pp2,1-lambda1,1-lambda2,1,pp1,pp2,0,0,0,0,0,lambda1,lambda2,0),
              nrow=3,ncol=5,byrow=T)
  
  # Output initiaux
  BE <- matrix(c(0,0,1,1,1,1,1,0,0,0,0,0,0,0,0),nrow=3,ncol=5,byrow=T) 
  
  # Matrice stochastique de transitions entre états (processus Markovien)
  # age 1
  PHI1 <- matrix(c(phi1c1,0,1-phi1c1,0,0,
                   0,phi1c2,0,1-phi1c2,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  # age 2
  PHI2 <- matrix(c(phi2c1,0,1-phi2c1,0,0,
                   0,phi2c2,0,1-phi2c2,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  PSI <- matrix(c(1,0,0,0,0,
                  0,1,0,0,0,
                  0,0,1,0,0,
                  0,0,0,1,0,
                  0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  A1 <- PHI1 %*% PSI
  A2 <- PHI2 %*% PSI
  
  # Distribution des états initiaux
  PI <- c(prop,1-prop,0,0,0)
  #---------------------------------------------
  
  #---------------------------------------------
  # on forme la vraisemblance, ou plutôt la log(vraisemblance) 
  # on linéarise car c'est plus stable numériquement de maximiser une somme qu'un produit
  # en fait, on forme -log(vraisemblance) plutôt, car R permet seulement de minimiser des fonctions
  # voir Pradel 2005, 2009
  
  l <- 0
  for (i in 1:nh) # boucles sur histoires de capture
  {
    ei <- e[i] # date marquage
    oe <- garb[i] + 1 # événement initial
    evennt <- data[,i] + 1 # les non-détections deviennent des 1, et les détections des 2
    ALPHA <- PI*BE[oe,]
    for (j in seqi(ei+1,lc[i])) # on conditionne par rapport à la première capture
    {
      if (j==(ei+1)) ALPHA <- (ALPHA %*% A1)*B[evennt[j],]
      if (j>(ei+1)) ALPHA <- (ALPHA %*% A2)*B[evennt[j],]
    }
    l <- l + logprot(sum(ALPHA))*eff[i]
  }
  l <- -l
  l
  #---------------------------------------------
  
} # fin fn vrais

#################################################################################
# pi, phi(a2*mix), p, r(mix) ##########################
#################################################################################
devCJS14 <- function(b,data,eff,e,lc,garb,nh,km1){
  
  #---------------------------------------------
  # lien logit pour les survies (transience)
  phi1c1 <- 1/(1+exp(-b[1]))
  phi2c1 <- 1/(1+exp(-b[2]))
  phi1c2 <- 1/(1+exp(-b[3]))
  phi2c2 <- 1/(1+exp(-b[4]))
  # lien logit pour les captures (homogénéité)
  pp <- 1/(1+exp(-b[5]))
  # lien logit pour les reprises (hétérogénéité)
  lambda1 <- 1/(1+exp(-b[6]))
  lambda2 <- 1/(1+exp(-b[7]))
  # prop
  prop <- 1/(1+exp(-b[8]))
  #---------------------------------------------
  
  
  #---------------------------------------------
  # Probabilités des événements (lignes) conditionnellement aux états (colonnes)
  B <- matrix(c(1-pp,1-pp,1-lambda1,1-lambda2,1,pp,pp,0,0,0,0,0,lambda1,lambda2,0),
              nrow=3,ncol=5,byrow=T)
  
  # Output initiaux
  BE <- matrix(c(0,0,1,1,1,1,1,0,0,0,0,0,0,0,0),nrow=3,ncol=5,byrow=T) 
  
  # Matrice stochastique de transitions entre états (processus Markovien)
  # age 1
  PHI1 <- matrix(c(phi1c1,0,1-phi1c1,0,0,
                   0,phi1c2,0,1-phi1c2,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  # age 2
  PHI2 <- matrix(c(phi2c1,0,1-phi2c1,0,0,
                   0,phi2c2,0,1-phi2c2,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  PSI <- matrix(c(1,0,0,0,0,
                  0,1,0,0,0,
                  0,0,1,0,0,
                  0,0,0,1,0,
                  0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  A1 <- PHI1 %*% PSI
  A2 <- PHI2 %*% PSI
  
  # Distribution des états initiaux
  PI <- c(prop,1-prop,0,0,0)
  #---------------------------------------------
  
  #---------------------------------------------
  # on forme la vraisemblance, ou plutôt la log(vraisemblance) 
  # on linéarise car c'est plus stable numériquement de maximiser une somme qu'un produit
  # en fait, on forme -log(vraisemblance) plutôt, car R permet seulement de minimiser des fonctions
  # voir Pradel 2005, 2009
  
  l <- 0
  for (i in 1:nh) # boucles sur histoires de capture
  {
    ei <- e[i] # date marquage
    oe <- garb[i] + 1 # événement initial
    evennt <- data[,i] + 1 # les non-détections deviennent des 1, et les détections des 2
    ALPHA <- PI*BE[oe,]
    for (j in seqi(ei+1,lc[i])) # on conditionne par rapport à la première capture
    {
      if (j==(ei+1)) ALPHA <- (ALPHA %*% A1)*B[evennt[j],]
      if (j>(ei+1)) ALPHA <- (ALPHA %*% A2)*B[evennt[j],]
    }
    l <- l + logprot(sum(ALPHA))*eff[i]
  }
  l <- -l
  l
  #---------------------------------------------
  
} # fin fn vrais


#################################################################################
# pi, phi(a2*mix), p(mix), r ##########################
#################################################################################
devCJS15 <- function(b,data,eff,e,lc,garb,nh,km1){
  
  #---------------------------------------------
  # lien logit pour les survies (transience)
  phi1c1 <- 1/(1+exp(-b[1]))
  phi2c1 <- 1/(1+exp(-b[2]))
  phi1c2 <- 1/(1+exp(-b[3]))
  phi2c2 <- 1/(1+exp(-b[4]))
  # lien logit pour les captures (hétérogénéité)
  pp1 <- 1/(1+exp(-b[5]))
  pp2 <- 1/(1+exp(-b[6]))
  # lien logit pour les reprises (hétérogénéité)
  lambda <- 1/(1+exp(-b[7]))
  # prop
  prop <- 1/(1+exp(-b[8]))
  #---------------------------------------------
  
  
  #---------------------------------------------
  # Probabilités des événements (lignes) conditionnellement aux états (colonnes)
  B <- matrix(c(1-pp1,1-pp2,1-lambda,1-lambda,1,pp1,pp2,0,0,0,0,0,lambda,lambda,0),
              nrow=3,ncol=5,byrow=T)
  
  # Output initiaux
  BE <- matrix(c(0,0,1,1,1,1,1,0,0,0,0,0,0,0,0),nrow=3,ncol=5,byrow=T) 
  
  # Matrice stochastique de transitions entre états (processus Markovien)
  # age 1
  PHI1 <- matrix(c(phi1c1,0,1-phi1c1,0,0,
                   0,phi1c2,0,1-phi1c2,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  # age 2
  PHI2 <- matrix(c(phi2c1,0,1-phi2c1,0,0,
                   0,phi2c2,0,1-phi2c2,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  PSI <- matrix(c(1,0,0,0,0,
                  0,1,0,0,0,
                  0,0,1,0,0,
                  0,0,0,1,0,
                  0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  A1 <- PHI1 %*% PSI
  A2 <- PHI2 %*% PSI
  
  # Distribution des états initiaux
  PI <- c(prop,1-prop,0,0,0)
  #---------------------------------------------
  
  #---------------------------------------------
  # on forme la vraisemblance, ou plutôt la log(vraisemblance) 
  # on linéarise car c'est plus stable numériquement de maximiser une somme qu'un produit
  # en fait, on forme -log(vraisemblance) plutôt, car R permet seulement de minimiser des fonctions
  # voir Pradel 2005, 2009
  
  l <- 0
  for (i in 1:nh) # boucles sur histoires de capture
  {
    ei <- e[i] # date marquage
    oe <- garb[i] + 1 # événement initial
    evennt <- data[,i] + 1 # les non-détections deviennent des 1, et les détections des 2
    ALPHA <- PI*BE[oe,]
    for (j in seqi(ei+1,lc[i])) # on conditionne par rapport à la première capture
    {
      if (j==(ei+1)) ALPHA <- (ALPHA %*% A1)*B[evennt[j],]
      if (j>(ei+1)) ALPHA <- (ALPHA %*% A2)*B[evennt[j],]
    }
    l <- l + logprot(sum(ALPHA))*eff[i]
  }
  l <- -l
  l
  #---------------------------------------------
  
} # fin fn vrais


#################################################################################
# pi, phi(a2*mix), p, r ##########################
#################################################################################
devCJS16 <- function(b,data,eff,e,lc,garb,nh,km1){
  
  #---------------------------------------------
  # lien logit pour les survies (transience)
  phi1c1 <- 1/(1+exp(-b[1]))
  phi2c1 <- 1/(1+exp(-b[2]))
  phi1c2 <- 1/(1+exp(-b[3]))
  phi2c2 <- 1/(1+exp(-b[4]))
  # lien logit pour les captures (hétérogénéité)
  pp <- 1/(1+exp(-b[5]))
  # lien logit pour les reprises (hétérogénéité)
  lambda <- 1/(1+exp(-b[6]))
  # prop
  prop <- 1/(1+exp(-b[7]))
  #---------------------------------------------
  
  
  #---------------------------------------------
  # Probabilités des événements (lignes) conditionnellement aux états (colonnes)
  B <- matrix(c(1-pp,1-pp,1-lambda,1-lambda,1,pp,pp,0,0,0,0,0,lambda,lambda,0),
              nrow=3,ncol=5,byrow=T)
  
  # Output initiaux
  BE <- matrix(c(0,0,1,1,1,1,1,0,0,0,0,0,0,0,0),nrow=3,ncol=5,byrow=T) 
  
  # Matrice stochastique de transitions entre états (processus Markovien)
  # age 1
  PHI1 <- matrix(c(phi1c1,0,1-phi1c1,0,0,
                   0,phi1c2,0,1-phi1c2,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  # age 2
  PHI2 <- matrix(c(phi2c1,0,1-phi2c1,0,0,
                   0,phi2c2,0,1-phi2c2,0,
                   0,0,0,0,1,
                   0,0,0,0,1,
                   0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  PSI <- matrix(c(1,0,0,0,0,
                  0,1,0,0,0,
                  0,0,1,0,0,
                  0,0,0,1,0,
                  0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  A1 <- PHI1 %*% PSI
  A2 <- PHI2 %*% PSI
  
  # Distribution des états initiaux
  PI <- c(prop,1-prop,0,0,0)
  #---------------------------------------------
  
  #---------------------------------------------
  # on forme la vraisemblance, ou plutôt la log(vraisemblance) 
  # on linéarise car c'est plus stable numériquement de maximiser une somme qu'un produit
  # en fait, on forme -log(vraisemblance) plutôt, car R permet seulement de minimiser des fonctions
  # voir Pradel 2005, 2009
  
  l <- 0
  for (i in 1:nh) # boucles sur histoires de capture
  {
    ei <- e[i] # date marquage
    oe <- garb[i] + 1 # événement initial
    evennt <- data[,i] + 1 # les non-détections deviennent des 1, et les détections des 2
    ALPHA <- PI*BE[oe,]
    for (j in seqi(ei+1,lc[i])) # on conditionne par rapport à la première capture
    {
      if (j==(ei+1)) ALPHA <- (ALPHA %*% A1)*B[evennt[j],]
      if (j>(ei+1)) ALPHA <- (ALPHA %*% A2)*B[evennt[j],]
    }
    l <- l + logprot(sum(ALPHA))*eff[i]
  }
  l <- -l
  l
  #---------------------------------------------
  
} # fin fn vrais


#################################################################################
# pi, phi(a2*t), p(mix), r(mix), psi ##########################
#################################################################################
devCJS17 <- function(b,data,eff,e,lc,garb,nh,km1){
  
  #---------------------------------------------
  # lien logit pour les survies (transience)
  phi1 <- 1/(1+exp(-b[8:(7+km1)]))
  phi2 <- 1/(1+exp(-b[(8+km1):(6+2*km1)]))
  # lien logit pour les captures (hétérogénéité)
  pp1 <- 1/(1+exp(-b[1]))
  pp2 <- 1/(1+exp(-b[2]))
  # lien logit pour les reprises (hétérogénéité)
  lambda1 <- 1/(1+exp(-b[3]))
  lambda2 <- 1/(1+exp(-b[4]))
  # prop
  prop <- 1/(1+exp(-b[5]))
  # trans entre classes d'hétérogénéité
  psi12 <- 1/(1+exp(-b[6]))
  psi21 <- 1/(1+exp(-b[7]))
  #---------------------------------------------
  
  #---------------------------------------------
  # Probabilités des événements (lignes) conditionnellement aux états (colonnes)
  B <- matrix(c(1-pp1,1-pp2,1-lambda1,1-lambda2,1,pp1,pp2,0,0,0,0,0,lambda1,lambda2,0),
              nrow=3,ncol=5,byrow=T)
  
  # Output initiaux
  BE <- matrix(c(0,0,1,1,1,1,1,0,0,0,0,0,0,0,0),nrow=3,ncol=5,byrow=T) 
  
  # Matrice stochastique de transitions entre états (processus Markovien)
  
  PSI <- matrix(c((1-psi12),psi12,0,0,0,
                  psi21,(1-psi21),0,0,0,
                  0,0,1,0,0,
                  0,0,0,1,0,
                  0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  
  # age 1
  PHI1 <- array(NA, dim = c(5,5,km1))
  A1 <- PHI1
  for (kk in 1:km1){
    PHI1[,,kk] <- matrix(c(phi1[kk],0,1-phi1[kk],0,0,
                           0,phi1[kk],0,1-phi1[kk],0,
                           0,0,0,0,1,
                           0,0,0,0,1,
                           0,0,0,0,1),nrow=5,ncol=5,byrow=T)
    A1[,,kk] <- PHI1[,,kk] %*% PSI
  }
  # age 2
  PHI2 <- array(NA, dim = c(5,5,km1-1))
  A2 <- PHI2
  for (kk in 1:(km1-1)){
    PHI2[,,kk] <- matrix(c(phi2[kk],0,1-phi2[kk],0,0,
                           0,phi2[kk],0,1-phi2[kk],0,
                           0,0,0,0,1,
                           0,0,0,0,1,
                           0,0,0,0,1),nrow=5,ncol=5,byrow=T)
    A2[,,kk] <- PHI2[,,kk] %*% PSI
  }
  
  # Distribution des états initiaux
  PI <- c(prop,1-prop,0,0,0)
  #---------------------------------------------
  
  #---------------------------------------------
  # on forme la vraisemblance, ou plutôt la log(vraisemblance) 
  # on linéarise car c'est plus stable numériquement de maximiser une somme qu'un produit
  # en fait, on forme -log(vraisemblance) plutôt, car R permet seulement de minimiser des fonctions
  # voir Pradel 2005, 2009
  
  l <- 0
  for (i in 1:nh) # boucles sur histoires de capture
  {
    ei <- e[i] # date marquage
    oe <- garb[i] + 1 # événement initial
    evennt <- data[,i] + 1 # les non-détections deviennent des 1, et les détections des 2
    ALPHA <- PI*BE[oe,]
    for (j in seqi(ei+1,lc[i])) # on conditionne par rapport à la première capture
    {
      if (j==(ei+1)) ALPHA <- (ALPHA %*% A1[,,j-1])*B[evennt[j],]
      if (j>(ei+1)) ALPHA <- (ALPHA %*% A2[,,j-2])*B[evennt[j],]
    }
    l <- l + logprot(sum(ALPHA))*eff[i]
  }
  l <- -l
  l
  #---------------------------------------------
  
} # fin fn vrais


#################################################################################
# pi, phi(a2*mix*t), p(mix), r(mix), psi ##########################
#################################################################################
devCJS18 <- function(b,data,eff,e,lc,garb,nh,km1){
  
  #---------------------------------------------
  # lien logit pour les survies (transience)
  phi1c1 <- 1/(1+exp(-b[8:(7+km1)]))
  phi1c2 <- 1/(1+exp(-b[(8+km1):(7+2*km1)]))
  phi2c1 <- 1/(1+exp(-b[(8+2*km1):(6+3*km1)]))
  phi2c2 <- 1/(1+exp(-b[(7+3*km1):(5+4*km1)]))
  
  # lien logit pour les captures (hétérogénéité)
  pp1 <- 1/(1+exp(-b[1]))
  pp2 <- 1/(1+exp(-b[2]))
  # lien logit pour les reprises (hétérogénéité)
  lambda1 <- 1/(1+exp(-b[3]))
  lambda2 <- 1/(1+exp(-b[4]))
  # prop
  prop <- 1/(1+exp(-b[5]))
  # trans entre classes d'hétérogénéité
  psi12 <- 1/(1+exp(-b[6]))
  psi21 <- 1/(1+exp(-b[7]))
  #---------------------------------------------
  
  
  #---------------------------------------------
  # Probabilités des événements (lignes) conditionnellement aux états (colonnes)
  B <- matrix(c(1-pp1,1-pp2,1-lambda1,1-lambda2,1,pp1,pp2,0,0,0,0,0,lambda1,lambda2,0),
              nrow=3,ncol=5,byrow=T)
  
  # Output initiaux
  BE <- matrix(c(0,0,1,1,1,1,1,0,0,0,0,0,0,0,0),nrow=3,ncol=5,byrow=T) 
  
  # Matrice stochastique de transitions entre états (processus Markovien)
  
  PSI <- matrix(c((1-psi12),psi12,0,0,0,
                  psi21,(1-psi21),0,0,0,
                  0,0,1,0,0,
                  0,0,0,1,0,
                  0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  
  # age 1
  PHI1 <- array(NA, dim = c(5,5,km1))
  A1 <- PHI1
  for (kk in 1:km1){
    PHI1[,,kk] <- matrix(c(phi1c1[kk],0,1-phi1c1[kk],0,0,
                           0,phi1c2[kk],0,1-phi1c2[kk],0,
                           0,0,0,0,1,
                           0,0,0,0,1,
                           0,0,0,0,1),nrow=5,ncol=5,byrow=T)
    A1[,,kk] <- PHI1[,,kk] %*% PSI
  }
  
  # age 2
  PHI2 <- array(NA, dim = c(5,5,km1-1))
  A2 <- PHI2
  for (kk in 1:(km1-1)){
    PHI2[,,kk] <- matrix(c(phi2c1[kk],0,1-phi2c1[kk],0,0,
                           0,phi2c2[kk],0,1-phi2c2[kk],0,
                           0,0,0,0,1,
                           0,0,0,0,1,
                           0,0,0,0,1),nrow=5,ncol=5,byrow=T)
    A2[,,kk] <- PHI2[,,kk] %*% PSI
  }
  
  # Distribution des états initiaux
  PI <- c(prop,1-prop,0,0,0)
  #---------------------------------------------
  
  #---------------------------------------------
  # on forme la vraisemblance, ou plutôt la log(vraisemblance) 
  # on linéarise car c'est plus stable numériquement de maximiser une somme qu'un produit
  # en fait, on forme -log(vraisemblance) plutôt, car R permet seulement de minimiser des fonctions
  # voir Pradel 2005, 2009
  
  l <- 0
  for (i in 1:nh) # boucles sur histoires de capture
  {
    ei <- e[i] # date marquage
    oe <- garb[i] + 1 # événement initial
    evennt <- data[,i] + 1 # les non-détections deviennent des 1, et les détections des 2
    ALPHA <- PI*BE[oe,]
    for (j in seqi(ei+1,lc[i])) # on conditionne par rapport à la première capture
    {
      if (j==(ei+1)) ALPHA <- (ALPHA %*% A1[,,j-1])*B[evennt[j],]
      if (j>(ei+1)) ALPHA <- (ALPHA %*% A2[,,j-2])*B[evennt[j],]
    }
    l <- l + logprot(sum(ALPHA))*eff[i]
  }
  l <- -l
  l
  #---------------------------------------------
  
} # fin fn vrais





#################################################################################
# pi, phi(a2*2period), p(mix), r(mix), psi ##########################
#################################################################################
devCJS19 <- function(b,data,eff,e,lc,garb,nh,km1){
  
  #---------------------------------------------
  # lien logit pour les survies (transience)
  phi1 <- 1/(1+exp(-b[8:9]))
  phi2 <- 1/(1+exp(-b[10:11]))
  # lien logit pour les captures (hétérogénéité)
  pp1 <- 1/(1+exp(-b[1]))
  pp2 <- 1/(1+exp(-b[2]))
  # lien logit pour les reprises (hétérogénéité)
  lambda1 <- 1/(1+exp(-b[3]))
  lambda2 <- 1/(1+exp(-b[4]))
  # prop
  prop <- 1/(1+exp(-b[5]))
  # trans entre classes d'hétérogénéité
  psi12 <- 1/(1+exp(-b[6]))
  psi21 <- 1/(1+exp(-b[7]))
  #---------------------------------------------
  
  #---------------------------------------------
  # Probabilités des événements (lignes) conditionnellement aux états (colonnes)
  B <- matrix(c(1-pp1,1-pp2,1-lambda1,1-lambda2,1,pp1,pp2,0,0,0,0,0,lambda1,lambda2,0),
              nrow=3,ncol=5,byrow=T)
  
  # Output initiaux
  BE <- matrix(c(0,0,1,1,1,1,1,0,0,0,0,0,0,0,0),nrow=3,ncol=5,byrow=T) 
  
  # Matrice stochastique de transitions entre états (processus Markovien)
  
  PSI <- matrix(c((1-psi12),psi12,0,0,0,
                  psi21,(1-psi21),0,0,0,
                  0,0,1,0,0,
                  0,0,0,1,0,
                  0,0,0,0,1),nrow=5,ncol=5,byrow=T)
  
  # age 1
  PHI1 <- array(NA, dim = c(5,5,km1))
  A1 <- PHI1
  for (kk in 1:18){
    PHI1[,,kk] <- matrix(c(phi1[1],0,1-phi1[1],0,0,
                           0,phi1[1],0,1-phi1[1],0,
                           0,0,0,0,1,
                           0,0,0,0,1,
                           0,0,0,0,1),nrow=5,ncol=5,byrow=T)
    A1[,,kk] <- PHI1[,,kk] %*% PSI
  }
  for (kk in 19:km1){
    PHI1[,,kk] <- matrix(c(phi1[2],0,1-phi1[2],0,0,
                           0,phi1[2],0,1-phi1[2],0,
                           0,0,0,0,1,
                           0,0,0,0,1,
                           0,0,0,0,1),nrow=5,ncol=5,byrow=T)
    A1[,,kk] <- PHI1[,,kk] %*% PSI
  }
  # age 2
  PHI2 <- array(NA, dim = c(5,5,km1-1))
  A2 <- PHI2
  for (kk in 1:18){
    PHI2[,,kk] <- matrix(c(phi2[1],0,1-phi2[1],0,0,
                           0,phi2[1],0,1-phi2[1],0,
                           0,0,0,0,1,
                           0,0,0,0,1,
                           0,0,0,0,1),nrow=5,ncol=5,byrow=T)
    A2[,,kk] <- PHI2[,,kk] %*% PSI
  }
  for (kk in 19:(km1-1)){
    PHI2[,,kk] <- matrix(c(phi2[2],0,1-phi2[2],0,0,
                           0,phi2[2],0,1-phi2[2],0,
                           0,0,0,0,1,
                           0,0,0,0,1,
                           0,0,0,0,1),nrow=5,ncol=5,byrow=T)
    A2[,,kk] <- PHI2[,,kk] %*% PSI
  }
  
  # Distribution des états initiaux
  PI <- c(prop,1-prop,0,0,0)
  #---------------------------------------------
  
  #---------------------------------------------
  # on forme la vraisemblance, ou plutôt la log(vraisemblance) 
  # on linéarise car c'est plus stable numériquement de maximiser une somme qu'un produit
  # en fait, on forme -log(vraisemblance) plutôt, car R permet seulement de minimiser des fonctions
  # voir Pradel 2005, 2009
  
  l <- 0
  for (i in 1:nh) # boucles sur histoires de capture
  {
    ei <- e[i] # date marquage
    oe <- garb[i] + 1 # événement initial
    evennt <- data[,i] + 1 # les non-détections deviennent des 1, et les détections des 2
    ALPHA <- PI*BE[oe,]
    for (j in seqi(ei+1,lc[i])) # on conditionne par rapport à la première capture
    {
      if (j==(ei+1)) ALPHA <- (ALPHA %*% A1[,,j-1])*B[evennt[j],]
      if (j>(ei+1)) ALPHA <- (ALPHA %*% A2[,,j-2])*B[evennt[j],]
    }
    l <- l + logprot(sum(ALPHA))*eff[i]
  }
  l <- -l
  l
  #---------------------------------------------
  
} # fin fn vrais



