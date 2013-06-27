setMethod("initialize", "Trio", function(.Object, ... ){
  .Object <- callNextMethod()
  .Object
})

setMethod("GeneDrop", signature( trio = "Trio", geno.vec = "numeric"), function( trio, geno.vec ){
  g1 <- geno.vec[trio@id]
  g2 <- geno.vec[trio@spouse]
  for (i in 1:length(trio@offspring)){
    if( is.character(trio@offspring[[i]]) ){
      if( is.na(geno.vec[ trio@offspring[[i]]])){
        goff <- gene.drop.fn(g1,g2)
        geno.vec[ trio@offspring[[i]] ] <- goff
      }
    }else{
      # this is the case where the offspring is a trio object
      if( is.na(geno.vec[ trio@offspring[[i]]@id])){
        goff <- gene.drop.fn(g1,g2)
        geno.vec[ trio@offspring[[i]]@id ] <- goff
        geno.vec <- GeneDrop(trio@offspring[[i]], geno.vec)
      }
    }
  }
  return( geno.vec )
})

GeneDropSim.fn <- function(trio, geno.vec, dt.vec, n = 1e3, k = 10, nf = 1){
  n.bail <- k*n; i <- 1;
  share.vec <- logical(n.bail); occur.vec <- logical(n.bail);
  while( sum(occur.vec) < n & i <= n.bail ){
      founder <- sample(c(1,2,3,6),nf,replace = FALSE)
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
