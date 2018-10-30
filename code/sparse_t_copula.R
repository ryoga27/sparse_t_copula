source("sparse_t.R")

any2unif = function(x){
    d = ncol(x)
    n = nrow(x)
    u = matrix(NA, nrow = n, ncol = d)
    for(j in 1:d){
        cdf = ecdf(x[, j])
        u[, j] = cdf(x[, j])
    }
    return(u)
}

unif2t = function(u, nu){
    d = ncol(u)
    n = nrow(u)
    u[u == 1] = 0.99
    t = matrix(NA, nrow = n, ncol = d)
    for(j in 1:d){
        t = qt(u, df = nu)
    }
    return(t)
}

Sigma2rho = function(Sigma){
    d = nrow(Sigma)
    rho = matrix(NA, nrow = d, ncol = d)
    sigma = sqrt(diag(Sigma))
    for(i in 1:d){
        for(j in 1:d){
            rho[i, j] = Sigma[i, j]/(sigma[i]*sigma[j])
        }
    }
    return(rho)
}

sparse_t_copula = function(
    x, nu, lambda = 1,
    iter_max = 100, epsilon = 1e-3, step.size = 0.1
){
    n = nrow(x)
    d = ncol(x)
    u = any2unif(x = x)
    t = unif2t(u = u, nu = nu)
    fit = sparse_t(x = t, nu = nu, lambda = lambda, iter_max = iter_max, epsilon = epsilon, step.size = step.size)
    rho = Sigma2rho(fit$Sigma)
    fit$rho = rho
    return(fit)
}
