gene.drop.fn <- function(g1,g2){
  if( is.na(g1) || is.na(g2) ){
    goff <- NA
  }else if( g1 == 0 && g2 == 1 ){
    goff <- sample( x = 0:1, size = 1, prob = rep(0.5,2) )
  }else if( g1 == 1 && g2 == 0 ){
    goff <- sample( x = 0:1, size = 1, prob = rep(0.5,2) )
  ## }else if( g1 == 0 && g2 == 1 ){
  ##   goff <- sample( x = 0:1, size = 1, prob = rep(0.5,2) )
  ## }else if( g1 == 0 && g2 == 1 ){
  ##   goff <- sample( x = 0:1, size = 1, prob = rep(0.5,2) )
  }else{
    goff <- NA
  }
  return(goff)
}
