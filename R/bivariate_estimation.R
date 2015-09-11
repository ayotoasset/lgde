
##############################################################
##                                                          ##
##   ESTIMATE THE LOCAL CORRELATION BETWEEN TWO VARIABLES   ##
##           WHEN THE MARGINALS ARE STANDARD NORMAL         ##
##                                                          ##
##             CV AND KNN BANDWIDTH SELECTORS               ##
##                                                          ##
##############################################################

#' Estimate the local correlation between two variables when the marginals are standard normal.
#' @param data The data matrix, one row per observation. The marginals are assumed to be standard normal. They can be transformed using the transLocal-function.
#' @param grid The grid where the local correlation will be estimated.
#' @param h The vector of bandwidths. Length must be 2.
#' @return A list with 4 elements: $data is the data, $grid is the grid, $par.est is a vector of the estimated local correlations, and $f.est are the local likelihood density estimates there.
#' @examples
#' data <- cbind(rt(100, df = 10), rt(100, df = 10))
#' grid <- matrix(c(-1, -1, 0, 0, 1, 1), ncol = 2, byrow = TRUE)
#' trans <- transLocal(data = data, grid = grid)
#'
#' estimates <- biLocal(data = trans$transformed.data,
#'                      grid = trans$transformed.grid,
#'                      h = c(1, 1))

biLocal = function(data, 
                   grid,
                   h) {
    # Tests
    if(length(h) != 2) stop("The bandwidth must have length 2")
    if(dim(data)[2] != 2) stop("The data must have dimension 2")
    if(dim(grid)[2] != 2) stop("The grid must have dimension 2")
    
    # The data
    x <- data[,1]
    y <- data[,2]

    # The bandwidths
    h1 <- h[1]
    h2 <- h[2]
    
    maximize.likelihood = function(grid.point) {
      
        x0 <- grid.point[1]
        y0 <- grid.point[2]
            
        # We need weights and some empirical moments in this grid point
        W <- dnorm(x, mean = x0, sd = h1)*dnorm(y, mean = y0, sd = h2)

        m1 <- mean(W)
        m2 <- mean(W*x^2)
        m3 <- mean(W*y^2)
        m4 <- mean(W*x*y)
            
        # The likelihood function
        lik <- function(rho) {
            - log(2*pi*sqrt(1-rho^2))*m1 - m2/(2*(1-rho^2)) - m3/(2*(1-rho^2)) + rho*m4/(1-rho^2) -
            1/2*exp(-1/2*(y0^2*h1^2+x0^2+x0^2*h2^2-2*x0*rho*y0+y0^2)/(-rho^2+h2^2+1+h1^2+h1^2*h2^2))/
            (pi*(-rho^2+h2^2+1+h1^2+h1^2*h2^2)^(1/2))
        }
            
        # Return the maximum of the likelihood and the density estimate
        opt <- try(optimise(lik,
                            lower = -1,
                            upper=1,
                            maximum = TRUE,
                            tol = .Machine$double.eps^0.25/10^4),
                   silent=TRUE) 
        if(class(opt)!="try-error") {
            return(c(opt$maximum, mvtnorm::dmvnorm(c(x0, y0), mean = c(0,0), sigma = matrix(c(1, opt$maximum, opt$maximum, 1), 2))))
        } else {
            return(c(NA, NA))
        }
      
    }

    ## Send the grid points to 'maximize.likelihood'
    est <- cbind(do.call(rbind, lapply(X = split(grid, row(grid)), FUN = maximize.likelihood)))
    par.est <- matrix(est[,1])
    colnames(par.est) = c('rho')
     
    return(list(data = data,
                grid = grid,
                par.est = par.est,
                f.est = est[,2]))
}    
                            
#' Estimate the local correlation between two variables when the marginals are standard normal. Uses local bandwidths, one bandwidth per grid point.
#' @param data The data matrix, one row per observation. The marginals are assumed to be standard normal. They can be transformed using the transLocal-function.
#' @param grid The grid where the local correlation will be estimated.
#' @param h The vector of bandwidths. Length must be equal to the number of grid points.
#' @return A list with 4 elements: $data is the data, $grid is the grid, $par.est is a vector of the estimated local correlations, and $f.est are the local likelihood density estimates there.
#' @examples
#' data <- cbind(rt(100, df = 10), rt(100, df = 10))
#' grid <- matrix(c(-1, -1, 0, 0, 1, 1), ncol = 2, byrow = TRUE)
#' trans <- transLocal(data = data, grid = grid)
#'
#' estimates <- biLocal(data = trans$transformed.data,
#'                      grid = trans$transformed.grid,
#'                      h = c(1, 1))


biLocal.knn = function(data, 
                       grid,
                       h) {

    # Check for matching dimensions
    if(dim(grid)[1] != length(h)) {
        stop(" One bandwidth for each grid point! Does not match!")
    }
    
    # The data
    x <- data[,1]
    y <- data[,2]

    maximize.likelihood = function(grid.point.index) {
      
        x0 <- grid[grid.point.index, 1]
        y0 <- grid[grid.point.index, 2]

        bw <- h[grid.point.index]
          
        # We need weights and some empirical moments in this grid point
        W <- dnorm(x, mean = x0, sd = bw)*dnorm(y, mean = y0, sd = bw)

        m1 <- mean(W)
        m2 <- mean(W*x^2)
        m3 <- mean(W*y^2)
        m4 <- mean(W*x*y)
            
        # The likelihood function
        lik <- function(rho) {
            - log(2*pi*sqrt(1-rho^2))*m1 - m2/(2*(1-rho^2)) - m3/(2*(1-rho^2)) + rho*m4/(1-rho^2) -
            1/2*exp(-1/2*(y0^2*bw^2+x0^2+x0^2*bw^2-2*x0*rho*y0+y0^2)/(-rho^2+bw^2+1+bw^2+bw^2*bw^2))/
            (pi*(-rho^2+bw^2+1+bw^2+bw^2*bw^2)^(1/2))
        }
            
        # Return the maximum of the likelihood and the density estimate
        opt <- try(optimise(lik,
                            lower = -1,
                            upper = 1,
                            maximum = TRUE,
                            tol = .Machine$double.eps^0.25/10^4),
                   silent=TRUE) 
        if(class(opt)!="try-error") {
            return(c(opt$maximum, dmvnorm(c(x0, y0), mean = c(0,0), sigma = matrix(c(1, opt$maximum, opt$maximum, 1), 2))))
        } else {
            return(c(NA, NA))
        }
      
    }

    ## Send the grid points to 'maximize.likelihood'
    est <- t(apply(matrix(1:dim(grid)[1], ncol = 1), 1, maximize.likelihood))
    par.est <- matrix(est[,1])
    colnames(par.est) = c('rho')
     
    return(list(data = data,
                grid = grid,
                par.est = par.est,
                f.est = est[,2]))
}    
                            

