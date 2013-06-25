setClass("Trio", representation = list( id = "character", spouse = "character", offspring = "list"))

setClassUnion("TrioOrChar", c("Trio","character"))


