BMM <- function(x, k, n.iter = 1e3) {
    x <- sort(x)
    x <- (x - min(x) + .Machine$double.eps) / (max(x) - min(x) + 2 * .Machine$double.eps)
    y <- rep(1:k, each = length(x) / k)
    y <- c(y, rep(k, length(x) - length(y)))
    params <- data.frame(a = rep(0, k), b = rep(0, k))
    probs <- matrix(nrow = length(x), ncol = k)

    for (i in 1:n.iter) {
        # E step.
        means <- tapply(x, y, mean)
        vars <- tapply(x, y, var)
        com <- means * (1 - means) / vars - 1
        params$a <- means * com
        params$b <- (1 - means) * com

        # M step.
        prior <- table(y) / length(y)
        for (p in 1:nrow(params)) probs[, p] <- dbeta(x, params$a[p], params$b[p]) * prior[p]
        probs <- apply(probs, 2, '/', rowSums(probs))
        y <- max.col(y)
    }
    list(y = y, params = params, prior = prior)
}