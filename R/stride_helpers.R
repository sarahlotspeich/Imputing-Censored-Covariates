Kern.FUN <- function(zz,zi,bw) ## returns an (n x nz) matrix
{
  out = (VTM(zz,length(zi))- zi)/bw
  norm.k = dnorm(out)/bw
  norm.k
}

VTM<-function(vc, dm){
  matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}

#' @export
sum.I <- function(yy,FUN,Yi,Vi=NULL)   ## sum I(yy FUN Yi)Vi
{
  if(FUN=="<"|FUN=="<=") { yy <- -yy; Yi <- -Yi}
  if(substring(FUN,2,2)=="=") yy <- yy + 1e-8 else yy <- yy - 1e-8
  pos <- rank(c(yy,Yi))[1:length(yy)] - rank(yy)
  if(is.null(Vi)){return(pos)}else{
    Vi <- cumsum2(as.matrix(Vi)[order(Yi),,drop=F])
    out <- matrix(0, nrow=length(yy), ncol=dim(as.matrix(Vi))[2])
    out[pos!=0,] <- Vi[pos,]
    if(is.null(dim(Vi))) out <- c(out)
    return(out) ## n.y x p
  }
}

cumsum2 <- function(mydat)     #cumsum by row, col remains the same
{
  if(is.null(dim(mydat))) return(cumsum(mydat))
  else{
    out <- matrix(cumsum(mydat), nrow=nrow(mydat))
    out <- out - VTM(c(0, out[nrow(mydat), -ncol(mydat)]), nrow(mydat))
    return(out)
  }
}
