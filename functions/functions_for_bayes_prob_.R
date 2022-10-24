###############################################
#
#  Cao et al. 2020
#  "A Bayesian method for synthesizing multiple diagnostic outcomes of COVID-19"
#
###############################################

# functions_for_bayes_prob

get.posterior.prob = function(test.info.mat, COVID.prev){
  this.numerator = get.numerator(test.info.mat = test.info.mat, COVID.prev = COVID.prev)
  this.denominator = get.denominator(test.info.mat = test.info.mat, COVID.prev = COVID.prev)
  #
  posterior.prob = this.numerator / this.denominator
  return(posterior.prob)
}


get.numerator = function(test.info.mat, COVID.prev){
  numerator.component.array.when.T = get.component.array.when.T(test.info.mat)
  numerator = prod(numerator.component.array.when.T)*COVID.prev
  # 
  return(numerator)
}


get.denominator = function(test.info.mat, COVID.prev){
  denominator.component.array.when.T = get.component.array.when.T(test.info.mat)
  denominator.component.array.when.F = get.component.array.when.F(test.info.mat)
  denominator = prod(denominator.component.array.when.T)*COVID.prev + prod(denominator.component.array.when.F)*(1-COVID.prev)
  # 
  return(denominator)
}


get.component.array.when.T = function(test.info.mat){
  test.info.mat = as.matrix(test.info.mat, nrow = 3)
  this.pred.rate.array = ifelse(test.info.mat[3,], test.info.mat[1,], 1-test.info.mat[1,])
  this.pred.rate.array = unlist(this.pred.rate.array)
  names(this.pred.rate.array) = NULL
  # this.pred.rate.array
  return(this.pred.rate.array)
}


get.component.array.when.F = function(test.info.mat){
  test.info.mat = as.matrix(test.info.mat, nrow = 3)
  this.pred.rate.array = ifelse(test.info.mat[3,], 1-test.info.mat[2,], test.info.mat[2,])
  this.pred.rate.array = unlist(this.pred.rate.array)
  names(this.pred.rate.array) = NULL
  # this.pred.rate.array
  return(this.pred.rate.array)
}
