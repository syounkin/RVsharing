setMethod("RVsharing",  signature(ped="pedigree", id = "missing",dad.id = "missing",mom.id = "missing"), function(ped){
    RVsharing.fn( id = ped$id, dad.id = ifelse(ped$findex!=0,ped$id[ped$findex],0), mom.id = ifelse(ped$mindex!=0,ped$id[ped$mindex],0) )
})
setMethod("RVsharing",  signature(ped = "missing", id="character",dad.id="character",mom.id="character"), function(id, dad.id, mom.id){
    RVsharing.fn( id = id, dad.id = dad.id, mom.id = mom.id )
})
