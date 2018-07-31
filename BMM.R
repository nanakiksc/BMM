BMM <- function(x, k, n.iter = 1e3) {
    y <- rep(1:k, each = length(x) / k)
    y <- c(y, rep(k, length(x) - length(y)))
    params <- data.frame(a = rep(0, k), b = rep(0, k))
    probs <- matrix(nrow = k, ncol = length(x))

    for (i in 1:n.iter) {
        # E step.
        means <- tapply(x, y, mean)
        vars <- tapply(x, y, var)
        com <- means * (1 - means) / vars - 1
        params$a <- means * com
        params$b <- (1 - means) * com

        # M step.
        prior <- table(y) / length(y)
        for (p in 1:nrow(params)) probs[p, ] <- dbeta(x, params$a[p], params$b[p]) * prior[p]
        probs <- t(apply(probs, 1, '/', colSums(probs)))
        y <- apply(probs, 2, which.max)
    }
    list(y = y, params = params, prior = prior)
}