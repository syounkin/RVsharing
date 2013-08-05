### R code from vignette source 'ped2trio.asciidoc'

###################################################
### code chunk number 1: junk
###################################################
par(mar=rep(0,4))


###################################################
### code chunk number 2: options
###################################################
  options(width=75, continue = " ")
  library("Bureau")


###################################################
### code chunk number 3: RVsharing.toy2
###################################################
id <- paste0("sub", 1:8)
fa.id <- c(NA,NA,NA,"sub1","sub1", NA,"sub3","sub5")
ma.id <- c(NA,NA,NA,"sub2","sub2", NA,"sub4","sub6")
test.ped.first.cousins <- pedigree(id = id, dadid = fa.id, momid = ma.id, sex = c(1,2,1,2,1,2,2,2))


###################################################
### code chunk number 4: plotped2
###################################################
plot(test.ped.first.cousins)

###################################################
### code chunk number 5: ped2trio
###################################################
trio.obj <- RVsharing:::ped2trio(test.ped.first.cousins)

geno.vec <- c(1,0,0,NA,NA,0,NA,NA)
names(geno.vec) <- paste0("sub",1:8)

RVsharing:::GeneDropSim.fn( trio = trio.obj, geno.vec = geno.vec, dt.vec = c("sub7","sub8"), fd.indices = c(1,2,3,6), n = 1e2)

## trios <- rep(list(NA),length(ped2trio.list))
## names(trios) <- names(ped2trio.list)
## trio.ids <- ped2trio.list[[names(trios)[1]]] # first trio
## trios[paste0("trio",trio.ids$id)] <- new("Trio", id = trio.ids$id, spouse = trio.ids$spouse, offspring = list(trio.ids$offspring))

## trios


## trio.obj.list <- list()
## k <- length(ped2trio.list)

## for( i in 1:(k-1) ){

##     trio.obj.list[[i]] <- new("Trio", id = ped2trio.list[[i]]$id, spouse = ped2trio.list[[i]]$spouse, offspring = ped2trio.list[[i]]$offspring)

## }

## trio.obj <- new("Trio", id = ped2trio.list[[k]]$id, spouse = ped2trio.list[[k]]$spouse, offspring = trio.obj.list )
