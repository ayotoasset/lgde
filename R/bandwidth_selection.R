################################################
##|                                          |##
##|          BANDWIDTH SELECTION FOR:        |##
##| MULTIVARIATE LOCAL LIKELIHOOD ESTIMATION |##
##|                                          |##
##|     GLOBAL CROSS-VALIDATED BANDWIDTHS    |##
##|              PLUGIN BANDWIDTHS           |##
##|            LOCAL kNN-BANDWIDTHS          |##
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

#' Create a bandwidth object using the plugin method h = 1.75n^{-1/6}
#' @param n The number of observations.
#' @param nvar The number of variables
#' @return A matrix with bandwidths to be used in the multiLocal- or condLocal-functions.
#' @examples
#' data <- cbind(rnorm(100), rnorm(100))
#' pluginLocal(data)

pluginLocal <- function(n, nvar) {
    joint.bandwidths <- cbind(t(combn(c(1:nvar), 2)), 0 * t(combn(c(1:nvar), 
        2)))
    colnames(joint.bandwidths) <- c("x1", "x2", "h1", "h2")
    joint.bandwidths[,c("h1", "h2")] <- 1.75*n^(-1/6)
    joint.bandwidths
}

#' Calculate the distance to the k'th nearest neighbour for a given bivariate data set, on a given grid.
#' @param data The data set, a nx2-matrix.
#' @param grid The grid, a dx2-matrix.
#' @param k The number of neighbours.
#' @return A vector of length d containing the distance to the k'th nearest neighbour for each grid point.
#' @examples
#' data <- cbind(rnorm(100), rnorm(100))
#' grid <- cbind(c(-1, 0, 1), c(-1, 0, 1))
#' k <- 10
#' knnBi(data, grid, k)

knnBi <- function(data, grid, k) {
    # Calculate the Eucledian distance between all pairs of data, and observations.
    # Deal with infinite grid values by putting them far out, and therefore at the bottom of the list everywhere.
    grid[grid == Inf] <- 100
    euclid <- fields::rdist(grid, data)
    knn <- function(grid.point.index) {
        euclid[grid.point.index, sort.int(euclid[grid.point.index,], index.return = TRUE)$ix[k]]
    }
    # Return the knn-value for each grid points 
    apply(matrix(1:dim(grid)[1], ncol = 1), 1, knn)
}

#' Uses cross validation to determine the k for kNN bandwidth selection.
#' @param data The data matrix.
#' @param test Vector of integers that shall be tested. Chooses the k that has the highest likelihood.
#' @return Matrix with 4 columns, one for each pair of variables. The first two columns indicate the
#' pair in question, the third column givs the selected k, and the last column returns a 1 if the selected
#' k is an endpoint for the "test"-vector. Perhaps one should change the search area?
#' @examples
#' data <- cbind(rnorm(100), rnorm(100))
#' kLocal(data, test = seq(20, 80, 3))

kLocal <- function(data, test) {

    n <- dim(data)[1]                    # Number of observations
    d <- dim(data)[2]                    # Dimension of the problem

    # Initialize the matrix
    k <- cbind(t(combn(c(1:d), 2)), 0*t(combn(c(1:d), 2)))
    colnames(k) <- c('x1', 'x2', 'k', 'endpoint')

    # Calculate k for each pair of observations
    for(i in 1:dim(k)[1]) {
        biv.data <- cbind(data[, k[i,'x1']], data[, k[i,'x2']])
        # CV function to optimize
        CV.bivariate <- function(k) {
            k = k + 1 # Correct for evaluating in the grid points
            h <- knnBi(biv.data, biv.data, k)
            obj <- function(j) {
                log(biLocal.knn(data = biv.data[-j,],
                                grid = t(as.matrix(biv.data[j,])),
                                h = h[j])$f.est)
            }
            -mean(do.call(rbind, lapply(X = as.list(1:n), FUN = obj)), na.rm = TRUE)
        }
        CV <- apply(matrix(test, ncol = 1), 1, CV.bivariate)
        ind <- which.min(CV)
        k[i, 3] <- test[ind]
        if((ind == 1) | (ind == length(test)))
            k[i, 4] = 1
      }
      return(k)
}




