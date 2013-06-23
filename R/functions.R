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
