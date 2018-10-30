if(require(mvtnorm) == FALSE){
    install.packages("mvtnorm", dependencies = TRUE)
}
library(mvtnorm)

if(require(spcov) == FALSE){
    install.packages("spcov", dependencies = TRUE)
}
library(spcov)


# calculate mahalanobis distance
delta = function(x, mu, Sigma){
    d_vec = rep(NA, length = n)
    for(i in 1:n){
        d_vec[i] = t(x[i, ] - mu)%*%solve(Sigma)%*%(x[i, ] - mu)
    }
    return(d_vec)
}

# average each column of x
mean_col = function(x){
    d = ncol(x)
    mu = rep(NA, length = d)
    for(j in 1:d){
        mu[j] = mean(x[, j])
    }
    return(mu)
}

# calculate tau
calc_tau = function(x, nu, mu, Sigma){
    d = ncol(x)
    tau = (nu + d)/(nu + delta(x, mu, Sigma))
    return(tau)
}

# calculate mu
calc_mu = function(x, tau){
    d = ncol(x)
    mu = rep(NA, length = d)
    for(j in 1:d){
        mu[j] = sum(tau*x[, j])/sum(tau)
    }
    return(mu)
}

# calculate S_tau
calc_S_tau = function(x, tau, mu){
    S_tau = t(tau*x - mu)%*%(tau*x - mu)
    return(S_tau)
}

# calculate likelihood functiond
log_lilelihood = function(x, mu, Sigma, nu, lambda){
    log_likelihood = 0
    for(i in 1:n){
        log_likelihood = log_likelihood + log(dmvt(x[1, ], delta = mu, sigma = Sigma, log = FALSE))
    }
    penalty = lambda*(sum(solve(Sigma)) - sum(diag(solve(Sigma))))
    loss = log_likelihood + penalty
    return(loss)
}

sparse_t_copula = function(
    x, nu, lambda = 0.1,
    iter_max = 100, epsilon = 1e-3
){
    n = nrow(x)
    d = ncol(x)

    mu_array = array(NA, dim = c(d, iter_max + 1))
    Sigma_array = array(NA, dim = c(d, d, iter_max + 1))
    tau_array = array(NA, dim = c(n, iter_max))
    S_tau_array = array(NA, dim = c(d, d, iter_max))
    l_vec = rep(NA, length = iter_max + 1)

    # set intial values
    mu_0 = mean_col(x)
    mu_array[, 1] = mu_0
    S_0 = cov(x)
    Sigma_array[, , 1] = solve(S_0)
    l_vec[1] = 0

    for(s in 1:iter_max){
        tau_array[, s] = calc_tau(x = x, nu = nu, mu = mu_array[, s], Sigma = Sigma_array[, , s])
        mu_array[, s + 1] = calc_mu(x = x, tau = tau_array[, s])
        S_tau_array[, , s] = calc_S_tau(x = x, tau = tau_array[, s], mu = mu_array[, s + 1])
        fit = spcov(Sigma = diag(d), S = S_tau_array[, , s], lambda = lambda, step.size = 0.1)
        Sigma_array[, , s + 1] = fit$Sigma
        l_vec[s + 1] = log_lilelihood(x, mu = mu_array[, s + 1], Sigma = Sigma_array[, , s + 1], nu = nu, lambda = lambda)
        convergence = (abs(l_vec[s + 1] - l_vec[s]) < epsilon)
        if(convergence){
            iter_num = s
            message("converged")
            break
        }
        cat("optimizing Sigma:", s, "\n")
    }

    args_list = list(
        x = x,
        mu = mu_array[, iter_num + 1],
        Sigma = Sigma_array[, , iter_num + 1],
        nu = nu,
        log_likelihood = l_vec[iter_num + 1],
        iter_num = iter_num
    )
    return(args_list)
}

# select nu (defree of freedom)
select_nu = function(x, nu, lambda){

    N = length(nu)

    fit = list()
    log_likelihood_vec = rep(NA, length = N)

    for(i in 1:N){
        cat("searching nu:", nu_vec[i], "\n")
        fit[[i]] = sparse_t_copula(x = x, nu = nu_vec[i], lambda = lambda)
        log_likelihood_vec[i] = fit[[i]]$log_likelihood
    }

    opt_index = (1:N)[log_likelihood_vec == max(log_likelihood_vec)]

    args_list = list(
        mu = fit[[opt_index]]$mu,
        Sigma = fit[[opt_index]]$Sigma,
        nu = nu[opt_index],
        lambda = lamdba,
        log_likelihood = log_likelihood_vec
    )

    return(args_list)
}
