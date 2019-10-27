#' Gaussian Dirichlet Process Mixture CLustering For Tall Data
#'
#' @param data n by p matrix wuith n observations in rows and p dimensions in columns.
#' @param coresets a list with 3 components
#'
#' @author Boris Hejblum, Paul Kirk
#' @importFrom NPflow DPMGibbsN
#' @importFrom stats kmeans
#'
#' @export
#'
#' @examples
#' mydata <- rbind(cbind(rnorm(1000), rnorm(n = 1000)),
#'                 cbind(rnorm(500, m=10), rnorm(n = 500, m=10)))
#' coresets <- kmeans(mydata, 100)[c("cluster", "centers", "size")]
#'
#' res <- bigdpclust(mydata)
#' table(res$cluster[1:1000])
#' table(res$cluster[1001:1500])

bigdpclust <- function(data, coresets = NULL, clumping_fn = stats::kmeans,
                       nclumps = min(1000, nrow(data)/10),
                       hyperG0 = NULL,
                       Ninit = 50, Nmcmc = 1000,
                       burnin = Nmcmc/5, thin = 2, loss_fn = "MBinderN"){


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
            coresets <- clumping_fn(data, nclumps)[c("centers", "size")]
        }
    }

    if(is.null(hyperG0)){
        # default hyperpriors
        hyperG0 <- list()
        hyperG0[["mu"]] <- colMeans(coresets$centers)
        hyperG0[["kappa"]] <- 0.001
        hyperG0[["nu"]] <- d + 2
        hyperG0[["lambda"]] <- diag(d)/10
    }

    Ninit <- min(Ninit, nrow(coresets$centers))

    res_mcmc <- NPflow::DPMGibbsN(z = t(coresets$centers), obs_weights = coresets$size,
                             hyperG0 = hyperG0, nbclust_init = Ninit, N = Nmcmc,
                             doPlot = TRUE)
    s <- summary(res_mcmc, burnin = burnin, thin = thin, lossFn = loss_fn)
    coresets_clust <- s$point_estim$c_est
    names(coresets_clust) <- rownames(coresets$centers)

    output <- list("cluster" = as.numeric(as.factor(as.numeric(coresets_clust[coresets$cluster]))),
                   "mcmc_output" = res_mcmc)
    return(output)
}
