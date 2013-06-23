gene.drop.fn <- function(g1,g2){
  if( is.na(g1) || is.na(g2) ){
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

## trio5 <- new("Trio", id = "sub8", spouse = "sub11", offspring = list("sub12") )
## trio4 <- new("Trio", id = "sub7", spouse = "sub9", offspring = list("sub10") )
## trio3 <- new("Trio", id = "sub4", spouse = "sub6", offspring = list(trio5) )
## trio2 <- new("Trio", id = "sub3", spouse = "sub5", offspring = list(trio4) )
## trio1 <- new("Trio", id = "sub1", spouse = "sub2", offspring = list(trio2, trio3) )
## geno.vec <- c(1,0,NA,NA,0,0,NA,NA,0,NA,0,NA)
## names(geno.vec) <- paste0("sub",1:12)
## GeneDrop(trio1, geno.vec)
