
#########################################
##                                     ##
##   MULTIVARIATE DENSITY ESTIMATION   ##
##     GLOBAL AND LOCAL BANDWIDTHS     ##
##        CONDITIONAL DENSITIES        ##
##                                     ##
#########################################

#' Estimate a multivariate density function using the LGDE.
#' @param data The data matrix, one row per observation.
#' @param bandwidths A matrix of bandwidths. Must be in the same format as is
#' produced by the HLocal-function.
#' @param gsize If grid is not provided, the density will be estimated at
#' gsize points going diagonally through the data.
#' @param grid The grid at which to estimate the density. Must be a matrix with
#' the same number of columns as the data matrix.
#' @return A list with 6 elements:
#'     $f.est are the density estimates at the grid points,
#'     $d is the dimension,
#'     $loc.cor are the local correlations,
#'     $h is the matrix of bandwidths that was used,
#'     $pairs is atwo-column matrix indicating all pairs of variables,
#'     $grid is the grid at which the density has been estimated.
#' @examples
#' data <- cbind(rt(100, df = 10), rt(100, df = 10))
#' multiLocal(data)

multiLocal <- function(data,
                       bandwidths = NULL,
                       gsize = 15,
                       grid = apply(data, 2, function(x) seq(quantile(x,.001),quantile(x,.999), length.out=gsize))) {
    
    # Sample size and number of variables
    n <- dim(data)[1]
    d <- dim(data)[2]

    # Check that the dimensions agree
    if(d != dim(grid)[2]) stop("Grid and data dimensions do not match!")

    # Transform the data and grid to standard normality
    transformed <- transLocal(data = data,
                              grid = grid,
                              return.normalizing.constants = TRUE)

    transformed.data <- transformed$transformed.data
    z.grid <- transformed$transformed.grid
    normalizing.constants <- transformed$normalizing.constants

    # Calculate the bandwidths if they have not been provided by the user
    if(!is.null(bandwidths)) {
        h <- bandwidths
    } else {
        h <- HLocal(transformed.data)
    }
            
    pairs <- h[,1:2]
    if(length(pairs) == 2) {pairs <- t(as.matrix(pairs))} 
    if(length(pairs) == 2) {n.pairs <- 1} else {n.pairs <- dim(pairs)[1]}

    
    # Function to calculate local correlation
    joint.densities = function(pair.number) {
        x <- transformed.data[,pairs[pair.number,]]
        bw <- h[pair.number, c('h1', 'h2')]; 

        rho <- biLocal(data = x,
                       grid = z.grid[,h[pair.number,1:2]],
                       h = bw)$par.est                                     
    }

    # Estmate the local correlations
    joints <- lapply(X = as.list(1:n.pairs), FUN = joint.densities)

    # Calculate the density estimates
    f.estimate = function(parameters) {
          x <- parameters[1:d]
          rho.vec <- parameters[(d+1):length(parameters)]
          sigma <- diag(d)
          for(j in 1:n.pairs) {
              var1 <- pairs[j,1]
              var2 <- pairs[j,2]
              rho <- rho.vec[j]
              sigma[var1, var2] <- rho
              sigma[var2, var1] <- sigma[var1, var2]
          }
         mvtnorm::dmvnorm(as.numeric(x), mean = rep(0, d), sigma = sigma)
     }
            
    if(length(grid[,1]) == 1) {
        parameter.matrix <- t(c(do.call(cbind, z.grid),
                                do.call(cbind, joints))) 
    } else {
        parameter.matrix <- cbind(z.grid, do.call(cbind, joints))
    }

    parameter.matrix <- split(parameter.matrix, row(parameter.matrix))
            
    f.est <- do.call(rbind, lapply(X = parameter.matrix, FUN = f.estimate))*apply(normalizing.constants, 1, prod)
            
    # Return results
    ret <- list(loc.cor = joints, f.est = f.est, h = h, dimension = d, pairs = pairs, grid = grid)
    return(ret)
        
}

#' Estimate a multivariate density function using the LGDE and kNN bandwidths.
#' @param data The data matrix, one row per observation.
#' @param k.mat A matrix of k. Must be in the same format as is
#' produced by the kLocal-function.
#' @param gsize If grid is not provided, the density will be estimated at
#' gsize points going diagonally through the data.
#' @param grid The grid at which to estimate the density. Must be a matrix with
#' the same number of columns as the data matrix.
#' @return A list with 6 elements:
#'     $f.est are the density estimates at the grid points,
#'     $d is the dimension,
#'     $loc.cor are the local correlations,
#'     $k.mat is the matrix of k's that was used,
#'     $pairs is a two-column matrix indicating all pairs of variables,
#'     $grid is the grid at which the density has been estimated,
#'     $test is a vector of k's to test for cross-validation.
#' @examples
#' data <- cbind(rt(100, df = 10), rt(100, df = 10))
#' multiLocal.knn(data)

multiLocal.knn <- function(data,
                           k.mat = NULL,
                           gsize = 15,
                           grid = apply(data, 2, function(x)
                               seq(quantile(x,.001),quantile(x,.999), length.out=gsize)),
                           test = seq(20, 100, by = 3)) {
    
    # Sample size and number of variables
    n <- dim(data)[1]
    d <- dim(data)[2]

    # Transform the data and grid to standard normality
    transformed <- transLocal(data = data,
                              grid = grid,
                              return.normalizing.constants = TRUE)

    transformed.data <- transformed$transformed.data
    z.grid <- transformed$transformed.grid
    normalizing.constants <- transformed$normalizing.constants

    # Calculate the bandwidths if they have not been provided by the user
    if(!is.null(k.mat)) {
        h <- k.mat
    } else {
        h <- kLocal(transformed.data, test = test)
    }
            
    pairs <- h[,1:2]
    if(length(pairs) == 2) {pairs <- t(as.matrix(pairs))} 
    if(length(pairs) == 2) {n.pairs <- 1} else {n.pairs <- dim(pairs)[1]}

    
    # Function to calculate local correlation
    joint.densities = function(pair.number) {
        x <- transformed.data[,pairs[pair.number,]]
        k <- h[pair.number, c('k')];
        bw <- knnBi(x,
                    z.grid[,pairs[pair.number,]],
                    k)

        rho <- biLocal.knn(data = x,
                           grid = z.grid[,h[pair.number,1:2]],
                           h = bw)$par.est                                     
    }

    # Estmate the local correlations
    joints <- lapply(X = as.list(1:n.pairs), FUN = joint.densities)

    # Calculate the density estimates
    f.estimate = function(parameters) {
          x <- parameters[1:d]
          rho.vec <- parameters[(d+1):length(parameters)]
          sigma <- diag(d)
          for(j in 1:n.pairs) {
              var1 <- pairs[j,1]
              var2 <- pairs[j,2]
              rho <- rho.vec[j]
              sigma[var1, var2] <- rho
              sigma[var2, var1] <- sigma[var1, var2]
          }
         mvtnorm::dmvnorm(as.numeric(x), mean = rep(0, d), sigma = sigma)
     }
            
    if(length(grid[,1]) == 1) {
        parameter.matrix <- t(c(do.call(cbind, z.grid),
                                do.call(cbind, joints))) 
    } else {
        parameter.matrix <- cbind(z.grid, do.call(cbind, joints))
    }

    parameter.matrix <- split(parameter.matrix, row(parameter.matrix))
            
    f.est <- do.call(rbind, lapply(X = parameter.matrix, FUN = f.estimate))*
                        apply(normalizing.constants, 1, prod)
            
    # Return results
    ret <- list(joints = joints,
                f.est = f.est,
                h = h,
                dimension = d,
                pairs = pairs,
                grid = grid,
                k.mat = k.mat)
    return(ret)
        
}

#' Estimate conditional densities using the LGDE
#' @param data The data matrix, one row per observation.
#' @param cond A vector of values to be conditioned on.
#' @param bandwidths A matrix of bandwidths. Must be in the same format as is
#' produced by the HLocal-function.
#' @param gsize If grid is not provided, the density will be estimated at
#' gsize points going diagonally through the data.
#' @param grid The grid at which to estimate the density. Must be a matrix with
#' the same number of columns as the data matrix.
#' @return A list with 1 element:
#'     $f.est.cond is the estimated conditional density.
#' @examples
#' # Three dimensional example. The conditioning variables must be the last columns in the data matrix.
#' Estimating the density of X1|X2 = 0, X3 = 0:
#' data <- cbind(rt(100, df = 10), rt(100, df = 10), rt(100, df = 10))
#' condLocal(data, cond = c(0,0))

condLocal <- function(data,
                      cond,
                      bandwidths = NULL,
                      gsize = 15,
                      grid = apply(as.matrix(data[, 1:(dim(data)[2] - length(cond))]),
                                   2, function(x) seq(quantile(x,.001),quantile(x,.999),
                                                      length.out = gsize))) {

    # Check if the dimensions of the grid, data and conditioning vector match
    if(dim(data)[2] != (dim(grid)[2] + length(cond))) {
        stop('The data matrix, grid matrix and conditioning vector do not match.')
    }

    # Sample size and number of variables
    n <- dim(data)[1]
    d <- dim(data)[2]

    # The grid where we must estimate the joint density. On the x-scale
    x.grid <- cbind(grid, matrix(rep(cond, dim(grid)[1]), ncol = length(cond), byrow = TRUE))

    # Transform the data and grid to standard normality
    transformed <- transLocal(data = data,
                              grid = x.grid,
                              return.normalizing.constants = TRUE)
    transformed.data <- transformed$transformed.data
    z.grid <- transformed$transformed.grid
    normalizing.constants <- transformed$normalizing.constants

    # Calculate the bandwidths if they have not been provided by the user
    if(!is.null(bandwidths)) {
        h <- bandwidths
    } else {
        h <- HLocal(transformed.data)
    }
            
    pairs <- h[,1:2]
    if(length(pairs) == 2) {pairs <- t(as.matrix(pairs))} 
    if(length(pairs) == 2) {n.pairs <- 1} else {n.pairs <- dim(pairs)[1]}

    # Function to calculate local correlation
    joint.densities = function(pair.number) {
        x <- transformed.data[,pairs[pair.number,]]
        bw <- h[pair.number, c('h1', 'h2')]; 

        rho <- biLocal(data = x,
                       grid = z.grid[,h[pair.number,1:2]],
                       h = bw)$par.est                                     
    }

    # Estmate the local correlations
    joints <- lapply(X = as.list(1:n.pairs), FUN = joint.densities)

    # Collect the grid and the local correlations        
    if(length(grid[,1]) == 1) {
        parameter.matrix <- t(c(do.call(cbind, z.grid[,1:(d-length(cond))]),
                                do.call(cbind, joints))) 
    } else {
        parameter.matrix <- cbind(z.grid[,1:(d-length(cond))], do.call(cbind, joints))
    }

    parameter.matrix <- split(parameter.matrix, row(parameter.matrix))
            
    
    # Calculate the *conditional* density estimates
    f.estimate = function(parameters) {
         x <- parameters[1:(d-length(cond))]
         rho.vec <- parameters[(d-length(cond)+1):length(parameters)]
         sigma <- diag(d)
         for(j in 1:n.pairs) {
              var1 <- pairs[j,1]
              var2 <- pairs[j,2]
              rho <- rho.vec[j]
              sigma[var1, var2] <- rho
              sigma[var2, var1] <- sigma[var1, var2]
         }
         # Partition sigma
         S11 <- as.matrix(sigma[1:(d - length(cond)), 1:(d-length(cond))])
         S22 <- as.matrix(sigma[(d - length(cond) +1):d, (d-length(cond) + 1):d])
         S12 <- as.matrix(sigma[1:(d - length(cond)),  (d-length(cond) + 1):d])
         dim(S12) <- c(d - length(cond), length(cond))
         S21 <- as.matrix(sigma[(d-length(cond) + 1):d, 1:(d - length(cond))])
         dim(S21) <- rev(dim(S12))
         mvtnorm::dmvnorm(as.numeric(x),
                          mean = S12%*%solve(S22)%*%as.matrix(z.grid[1,(d-length(cond)+1):d]),
                          sigma = S11 - S12%*%solve(S22)%*%S21)
     }

    f.est.cond <- do.call(rbind,
                     lapply(X = parameter.matrix, FUN = f.estimate))*
                         apply(as.matrix(normalizing.constants[,1:(d-length(cond))]), 1, prod)

    ret <- list(f.est.cond = f.est.cond)
    return(ret)
}    
    

