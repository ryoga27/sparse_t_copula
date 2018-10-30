rm(list = ls(all = TRUE))
source("sparse_t.R")

data_path = "../data/stock_us/stock_us_20170606.csv"
data_set = read.csv(data_path)

data_set_label = data_set[, 1]

remove_na = function(data_set){
    d = ncol(data_set)
    na_index = rep(NA, length = d)
    for(j in 1:d){
        na_index[j] = ifelse(all(!is.na(data_set[, j])), TRUE, FALSE)
    }
    return(data_set[, na_index])
}

X = remove_na(data_set[, -1])
head(X[, 1:10])

n = nrow(X)
d = ncol(X)
u = matrix(NA, nrow = n, ncol = d)
for(j in 1:d){
    cdf = ecdf(X[, j])
    u[, j] = cdf(X[, j])
}

nu_vec = seq(from = 0, to = 10, by = 1)
fit = select_nu(x = u[, 1:100], nu = nu_vec, lambda = 0.1)
