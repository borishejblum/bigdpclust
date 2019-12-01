#' Gaussian Dirichlet Process Mixture CLustering For Tall Data
#'
#' @param data n by p matrix wuith n observations in rows and p dimensions in columns.
#'
#' @param coresets a list with 3 components
#'
#' @param variance logical flag indicating whether the covariance matrix of each cluster is
#' constrained as diagonal, or unconstrained full matrix.
#' Default is \code{FALSE} (unconstrained covariance).
#'
#' @param burnin an integer giving the number of MCMC iterations to burn. Ddefaults is half)
#'
#' @param plotevery_nit an integer indicating the interval between plotted iterations
#' when \code{doPlot} is \code{TRUE}. Default is \code{Nmcmc/10}
#'
#'@param doPlot logical flag indicating whether to plot MCMC iteration or not.
#'Default is \code{FALSE}.
#'
#'@param verbose logical flag indicating whether partition info is messaged over
#'at each MCMC iteration. Default is \code{FALSE}.
#'
#' @author Boris Hejblum, Paul Kirk
#' @importFrom NPflow DPMGibbsN
#' @importFrom stats kmeans
#'
#' @export
#'
#' @examples
#' n1 <- 50000
#' n2 <- 500
#' mydata <- rbind(cbind(rnorm(n1), rnorm(n = n1)),
#'                 cbind(rnorm(n2, m=10), rnorm(n = n2, m=10)))
#' #plot(mydata)
#' #coresets <- stats::kmeans(mydata, centers = 100)[c("cluster", "centers", "size")]
#'
#' res <- bigdpclust(mydata, nclumps=200)
#' table(res$cluster[1:n1])
#' table(res$cluster[n1 + 1:n2])

bigdpclust <- function(data, coresets = NULL, clumping_fn = stats::kmeans,
                       nclumps = min(500, nrow(data)/10),
                       hyperG0 = NULL,
                       Ninit = 50, Nmcmc = 1000,
                       burnin = Nmcmc/5, thin = 2, loss_fn = "MBinderN",
                       diagVar = FALSE,
                       plotevery_nit = Nmcmc/10,
                       doPlot = FALSE,
                       verbose = FALSE){


    if(!is.null(coresets)){
        stopifnot(is.list(coresets))
        stopifnot(names(coresets) == c("cluster", "centers", "size"))
        stopifnot(is.matrix(coresets$centers))
        stopifnot(is.vector(coresets$size))
        stopifnot(nrow(coresets$centers) == length(coresets$size))
    }

    if(!is.null(data)){
        stopifnot(is.matrix(data))
    }


    if(is.null(data)){
        if(is.null(coresets)){
            stop("either 'coresets' or 'data' must be specified")
        }
        else{
            d <- ncol(coresets$centers)
        }

    }
    else{
        d <- ncol(data)
    }


    if(is.null(coresets)){
        if(is.null(clumping_fn)){
            stop("both 'coresets' and 'clumping_fn' arguments are both NULL")
        }
        else{
            coresets <- suppressWarnings(clumping_fn(data, nclumps)[c("cluster", "centers", "size")])
        }
    }

    if(is.null(hyperG0)){
        # setting default hyperpriors
        hyperG0 <- list()
        hyperG0[["mu"]] <- colMeans(coresets$centers) # Empirical Bayes
        hyperG0[["kappa"]] <- 0.001 # Weakly informative
        hyperG0[["nu"]] <- d + 2 # Weakly informative
        hyperG0[["lambda"]] <- diag(d)/10
    }

    Ninit <- min(Ninit, nrow(coresets$centers))

    res_mcmc <- NPflow::DPMGibbsN(z = t(coresets$centers), obs_weights = coresets$size,
                             hyperG0 = hyperG0, nbclust_init = Ninit, N = Nmcmc,
                             doPlot = doPlot, diagVar = FALSE, plotevery = plotevery_nit,
                             verbose = verbose)
    s <- summary(res_mcmc, burnin = burnin, thin = thin, lossFn = loss_fn)
    coresets_clust <- s$point_estim$c_est
    names(coresets_clust) <- rownames(coresets$centers)

    output <- list("cluster" = as.numeric(as.factor(as.numeric(coresets_clust[coresets$cluster]))),
                   "mcmc_output" = res_mcmc)
    return(output)
}
