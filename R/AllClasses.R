setClass("Trio", representation = list( id = "character", spouse = "character", offspring = "list"))
setClass("RVsharing.prob", representation = list( pshare="numeric",iancestors="character",desfounders="list",id="character",dad.id="character",mom.id="character")

setClassUnion("TrioOrChar", c("Trio","character"))


