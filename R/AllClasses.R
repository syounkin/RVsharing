setClass("Trio", slots = c( id = "character", spouse = "character", offspring = "list"))

setClassUnion("TrioOrChar", c("Trio","character"))


