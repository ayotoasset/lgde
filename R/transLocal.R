
#######################################################
##                                                   ##
##   TRANSFORM DATA AND GRID TO STANDARD NORMALITY   ##
##        AND RETURN NORMALIZING CONSTANTS           ##
##                                                   ##
#######################################################

##      Author:        Håkon Otneim
##      Date:          01.12.2014
##      Description:   Takes a data matrix (one observation per row) and transforms each
##                     observation vector to standard normality using the logspline desnity
##                     estimator by Stone, Hansen and Kooperberg (1997).
##
##      Arguments:     data:                 The data matrix
##                     grid:                 Provide the grid if it should be transformed as well
##                     return.normalizing.constants:
##                                           Logical. If true, return the normalizing constants as well.
##
##      Value:         transformed.data:          The transformed data matrix.
##                     transformed.grid:          The transformed grid.
##                     normalizing.constants:     The normalizing constants for use in 'multiLocal'.
##
##      Necessary packages:
##                     logspline
##
##      References:    Stone, Charles J., et al. "Polynomial splines and their tensor products in extended linear modeling: 1994 Wald memorial lecture."
##                         The Annals of Statistics 25.4 (1997): 1371-1470.

transLocal <- function(data,
                       grid = NULL,
                       return.normalizing.constants = FALSE) {

    # The sample size and the dimension
    n <- dim(data)[1]
    d <- dim(data)[2]
                       
    # Estimate the marginals using the logspline
    estimate.marginal <- function(i) logspline(data[,i])
    marginal.estimates <- lapply(X = as.list(1:d), FUN = estimate.marginal)
            
    transform.data <- function(i) qnorm(plogspline(data[,i], marginal.estimates[[i]]))
    transformed.data <- matrix(unlist(lapply(X = as.list(1:d), FUN = transform.data)), ncol = d)

    ret = list(transformed.data = transformed.data)
    
    # Transform the grid points
    if(!is.null(grid)) {
        transform.grid <- function(i) qnorm(plogspline(grid[,i], marginal.estimates[[i]]))
        transformed.grid <- matrix(unlist(lapply(X = as.list(1:d), FUN = transform.grid)), ncol = d)
        ret$transformed.grid <- transformed.grid
    }

    # Calculate the normalizing constants
    if(return.normalizing.constants) {
        calculate.normalizing.constants <- function(i) dlogspline(grid[,i], marginal.estimates[[i]])/dnorm(qnorm(plogspline(grid[,i], marginal.estimates[[i]])))
        normalizing.constants <- matrix(unlist(lapply(X = as.list(1:d), FUN = calculate.normalizing.constants)), ncol = d)
        ret$normalizing.constants <- normalizing.constants
    }         

    return(ret)
}
