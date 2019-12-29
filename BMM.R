# TODO: use log probabilities.

# Paramaters.
# x: contains the values of the variable. It is assumed to be a mixture of independent Beta distributions.
# k: the number of desired components in the Beta Mixture Model.
# n.iter: the number of iterations of the Expectacion Maximization algorithm. [default: 1e3]

# Value:
# y: identity of the component to which each observation most likely belongs. It is an integer in the range [1, k].
# params: estimated `a` and `b` paremeters of each Beta component.
# prior: mixture coefficients of the BMM (i.e. the prior probabilities of each Beta component).

BMM <- function(x, k, n.iter = 1e3) {
    # Scale observations to the Beta domain.
    x <- (x - min(x) + .Machine$double.eps) / (max(x) - min(x) + 2 * .Machine$double.eps)
    
    # Sort observations to make component initialization easier.
    o <- order(x); r <- rank(x); x <- x[o]
    
    # Initialize components by grouping observations by value into clusters of equal size.  
    y <- rep(1:k, each = length(x) / k)
    y <- c(y, rep(k, length(x) - length(y))) # Fill up the last component with the remaining (highest-valued) observations.
    
    # Initialize BMM parameters and posterior probabilities.
    params <- data.frame(a = rep(0, k), b = rep(0, k))
    probs <- matrix(nrow = length(x), ncol = k)

    for (i in 1:n.iter) {
        # E step.
        means <- tapply(x, y, mean)
        vars <- tapply(x, y, var)
        com <- means * (1 - means) / vars - 1 # Common part of the expression used to compute `a` and `b`.
        params$a <- means * com
        params$b <- (1 - means) * com

        # M step.
        prior <- (table(y) + 1) / (length(y) + k) # Laplace estimate to avoid 0- or 1-valued prior probabilities.
        for (p in 1:nrow(params)) probs[, p] <- dbeta(x, params$a[p], params$b[p]) * prior[p] # Unnormalized posterior probabilities.
        probs <- apply(probs, 2, '/', rowSums(probs)) # Normalize posterior probabilities.
        y <- max.col(probs) # Select most likely component.
    }
    list(y = y[r], params = params, prior = prior)
}
