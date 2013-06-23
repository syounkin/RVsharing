setMethod("initialize", "Trio", function(.Object, ... ){
#  .Object@pedigree <- pedigree
  .Object <- callNextMethod()
  .Object
})

setMethod("GeneDrop", signature( trio = "Trio", geno.vec = "numeric"), function( trio, geno.vec ){
  g1 <- geno.vec[trio@id]
  g2 <- geno.vec[trio@spouse]
  if( is.na(g1) | is.na(g2) ) stop
  for (i in 1:length(trio@offspring)){
    if( is.character(trio@offspring[[i]]) ){
      goff <- gene.drop.fn(g1,g2)
      geno.vec[ trio@offspring[[i]] ] <- goff
    }else{
      # this is the case where the offspring is a trio object
      goff <- gene.drop.fn(g1,g2)
      geno.vec[ trio@offspring[[i]]@id ] <- goff
      geno.vec <- GeneDrop(trio@offspring[[i]], geno.vec)
    }
  }
  return( geno.vec )
})
