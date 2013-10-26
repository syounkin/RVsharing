setClass("Trio", representation = list( id = "character", spouse = "character", offspring = "list"))
setClass("RVsharing_prob", representation = list( pshare="numeric",iancestors="character",desfounders="list",id="character",dad.id="character",mom.id="character"))

setClassUnion("TrioOrChar", c("Trio","character"))


