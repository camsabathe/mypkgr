
# attention on a bien une matrice d'entrée
# de la forme p lignes x n colonnes
# les observations sont en colonne'

# x une matrice, à n colonnes (les observations) et p
#
# lignes
#
# mean un vecteur de moyennes
#
# varcovM une matrice de variance-covariance
#
# Log un paramètre logique valant TRUE par défaut

# mvnpdf<-function(x, mean, varcovM, Log=TRUE){
#   x_transp<-t(x)
#
#   n<-nrow(x_transp)
#   p<-ncol(x_transp)
#   det_sigma<-det(varcovM)
#   inv_sigma<-solve(varcovM)
#
#   f<-c()
#   for (i in 1:n){
#     prod1<-1/( (2*pi)**(p/2) * sqrt(det_sigma))
#     prod2<- exp( -1/2 * t(x_transp[i,]-mean) %*% inv_sigma %*% (x_transp[i,]-mean)   )
#     f[i]<-prod1*prod2
#   }
#
#
#
#   if (Log==TRUE){
#     f<-log(f)
#   }
#   return( list(x=x, y=(f)))
# }


#' fonction normale multivariee pour la densite
#'  description
#' @param x les points
#' @param mean les moyennes
#' @param varcovM la matrice de var cov sigma
#' @param Log si on veut les résultats en log, par défaut = TRUE
#'
#' @return les densites
#' @export
#'
#' @examples
#' mvnpdf(x=matrix(1.96), Log=FALSE)
#'dnorm(1.96)
mvnpdf <- function(x, mean =  rep(0, nrow(x)),
                   varcovM = diag(nrow(x)), Log = TRUE) {
  n <- ncol(x)
  p <- nrow(x)
  x0 <- x - mean
  Rinv <- solve(varcovM)
  LogDetvarcovM <- log(det(varcovM))

  y <- NULL
  for (j in 1:n) {
    yj <- - p/2 * log(2*pi) - 0.5 * LogDetvarcovM -
      0.5 * t(x0[, j]) %*% Rinv %*% x0[, j]
    y <- c(y, yj)
  }

  if (!Log) {
    y <- exp(y)
  }

  return( list(x=x, y=(y)))
}
#
#
# mvnpdf(x=t(x), mean=mean,varcovM = sigma, Log=FALSE)

