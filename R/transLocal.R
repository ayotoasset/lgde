
#######################################################
##                                                   ##
##   TRANSFORM DATA AND GRID TO STANDARD NORMALITY   ##
##        AND RETURN NORMALIZING CONSTANTS           ##
##                                                   ##
#######################################################

#' Transform the marginals of a multivariate data set to standard normality based
#' on the logspline density estimator.
#' @param data The data matrix, one row per observation.
#' @param grid The grid if it should be transformed as well.
#' @param return.normalizing.constants TRUE if the normalizing constants to be
#' used in the LGDE should be returned as well.
#' @return A list containing the transformed data ($transformed.data), the transformed
#' grid if provided ($transformed.grid), and the normalizing constants if chosen
#' ($normalizing.constants).
#' @examples
#' data <- cbind(rt(100, df = 10), rt(100, df = 10))
#' grid <- matrix(c(-1, -1, 0, 0, 1, 1), ncol = 2, byrow = TRUE)
#' 
#' transLocal(data, grid, return.normalizing.constants = TRUE)

transLocal <- function(data,
                       grid = NULL,
                       return.normalizing.constants = FALSE) {

    # The sample size and the dimension
    n <- dim(data)[1]
    d <- dim(data)[2]
                       
    # Estimate the marginals using the logspline
    estimate.marginal <- function(i) logspline::logspline(data[,i])
    marginal.estimates <- lapply(X = as.list(1:d), FUN = estimate.marginal)

    # Transform each observation vector
    transform.data <- function(i) qnorm(logspline::plogspline(data[,i],
                                                              marginal.estimates[[i]]))
    transformed.data <- matrix(unlist(lapply(X = as.list(1:d),
                                             FUN = transform.data)), ncol = d)

    # The list to be returned
    ret = list(transformed.data = transformed.data)
    
    # Transform the grid points
    if(!is.null(grid)) {
        transform.grid <- function(i)
            qnorm(logspline::plogspline(grid[,i], marginal.estimates[[i]]))
        transformed.grid <- matrix(unlist(lapply(X = as.list(1:d),
                                                 FUN = transform.grid)), ncol = d)
        ret$transformed.grid <- transformed.grid
    }

    # Calculate the normalizing constants
    if(return.normalizing.constants) {
        calculate.normalizing.constants <- function(i) {
            logspline::dlogspline(grid[,i], marginal.estimates[[i]])/
                dnorm(qnorm(plogspline(grid[,i], marginal.estimates[[i]])))
        }
        normalizing.constants <- matrix(unlist(lapply(X = as.list(1:d),
                                                      FUN = calculate.normalizing.constants)),
                                        ncol = d)
        ret$normalizing.constants <- normalizing.constants
    }         

    return(ret)
}
