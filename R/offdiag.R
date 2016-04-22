
offdiag <- function(m){
  # returns the unique elements of a symmetric matrix that are off the diagonals. 
  vm = (upper.tri(m))*m
  bools = unlist(upper.tri(m))
  values = unlist(vm)
  values = values[bools]
  return(values)
}