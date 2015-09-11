
#########################################
##                                     ##
##   MULTIVARIATE DENSITY ESTIMATION   ##
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
#'     $f.est are the density estimates at the grid points
#'     $d is the dimension
#'     $loc.cor are the local correlations
#'     $h is the matrix of bandwidths that was used
#'     $pairs is atwo-column matrix indicating all pairs of variables
#'     $grid is the grid at which the density has been estimated
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
