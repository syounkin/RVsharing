gene.drop.fn <- function(g1,g2){
  if( is.na(g1) | is.na(g2) ){
    goff <- NA
  }else if( g1 == 0 && g2 == 1 ){
    goff <- sample( x = 0:1, size = 1, prob = rep(0.5,2) )
  }else if( g1 == 1 && g2 == 0 ){
    goff <- sample( x = 0:1, size = 1, prob = rep(0.5,2) )
  }else if( g1 == 1 && g2 == 1 ){
    goff <- sample( x = 0:2, size = 1, prob = c(0.25,0.5,0.25) )
  }else if( g1 == 0 && g2 == 0 ){
    goff <- 0
  }else if( g1 == 2 && g2 == 0 ){
    goff <- 1
  }else if( g1 == 0 && g2 == 2 ){
    goff <- 1
  }else if( g1 == 1 && g2 == 2 ){
    goff <- sample( x = 1:2, size = 1, prob = rep(0.5,2) )
  }else if( g1 == 2 && g2 == 1 ){
    goff <- sample( x = 1:2, size = 1, prob = rep(0.5,2) )
  }else if( g1 == 2 && g2 == 2 ){
    goff <- 2
  }else{
    goff <- "error"
  }
  return(goff)
}

GeneDropSim.fn <- function(trio, geno.vec, dt.vec, fd.indices, n = 1e3, k = 10, nf = 1){
  n.bail <- k*n; i <- 1;
  share.vec <- logical(n.bail); occur.vec <- logical(n.bail);
  while( sum(occur.vec) < n & i <= n.bail ){
      founder <- sample(fd.indices,nf,replace = FALSE)
      geno.vec[founder] <- 1
      geno.vec.sim <- GeneDrop(trio, geno.vec)
      if( any(is.na(geno.vec.sim))) stop("GeneDrop returned NA genotpye.")
      share.vec[i] <- all( geno.vec.sim[dt.vec]==1 )
      occur.vec[i] <- any( geno.vec.sim[dt.vec]==1 )
      geno.vec[founder] <- 0
      i <- i + 1
    }
  return(sum(share.vec)/sum(occur.vec)) # note that these vectors have length n.bail
}
