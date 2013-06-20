setMethod("RVsharing",  signature(ped="pedigree", id = "missing",dad.id = "missing",mom.id = "missing"), function(ped){

  fa.index <- ifelse(ped$findex !=0, ped$findex, NA)
  fa.id <- ped$id[fa.index]
  ma.index <- ifelse(ped$mindex !=0, ped$mindex, NA)
  ma.id <- ped$id[ma.index]
  
    RVsharing.fn( id = ped$id, dad.id = fa.id, mom.id = ma.id )
})
setMethod("RVsharing",  signature(ped = "missing", id="character",dad.id="character",mom.id="character"), function(id, dad.id, mom.id){
    RVsharing.fn( id = id, dad.id = dad.id, mom.id = mom.id )
})
