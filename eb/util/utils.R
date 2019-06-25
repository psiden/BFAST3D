### Some helpful functions for multigrid
#
# Author: Per Siden
# Created 2017-08-28
###


get.lattice.points <- function(x,y)
{
  sz = c(length(x),length(y))
  X=matrix(rep(x,sz[2]),nrow=sz[1])
  Y=matrix(rep(y,each=sz[1]),nrow=sz[1])
  dim(X) = c(prod(sz),1)
  dim(Y) = c(prod(sz),1)
  points = cbind(X,Y)
  return(points)
}

