################################################
##|                                          |##
##|          BANDWIDTH SELECTION FOR:        |##
##| MULTIVARIATE LOCAL LIKELIHOOD ESTIMATION |##
##|                                          |##
##|     GLOBAL CROSS-VALIDATED BANDWIDTHS    |##
##|                                          |##
################################################

#' Calculate the global cross-validated bandwidths for the LGDE.
#' @param data The data matrix, one row per observation.
#' @return A matrix. The first two columns contain each variable pair,
#' the last two columns contain the two bandwidths for the pair
#' in question. Suitable for multiLocal-function.
#' @examples
#' data <- cbind(rnorm(100), rnorm(100))
#' HLocal(data)

HLocal <- function(data) {

    n <- dim(data)[1]                    # Number of observations
    d <- dim(data)[2]                    # Dimension of the problem

    # Initialize the matrix
    joint.bandwidths <- cbind(t(combn(c(1:d), 2)), 0*t(combn(c(1:d), 2)))
    colnames(joint.bandwidths) <- c('x1', 'x2', 'h1', 'h2')

    # Calculate the bandwisth for each pair of variables
    for(i in 1:dim(joint.bandwidths)[1]) {
        biv.data <- cbind(data[,joint.bandwidths[i,'x1']], data[,joint.bandwidths[i,'x2']])
        # CV function to optimize
        CV.bivariate <- function(h) {
            obj <- function(j) {
            log(biLocal(data = biv.data[-j,],
                        grid = t(as.matrix(biv.data[j,])),
                        h = h)$f.est)
            }
            -mean(do.call(rbind, lapply(X = as.list(1:n), FUN = obj)), na.rm = TRUE)
        }
        joint.bandwidths[i, 3:4] <- optim(c(1, 1),
                                         CV.bivariate)$par
      }
      return(joint.bandwidths)
}

