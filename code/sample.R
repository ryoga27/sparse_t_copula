rm(list = ls(all = TRUE))
source("sparse_t_copula.R")

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

x = remove_na(data_set[, -1])
n = nrow(x)
y = log(x[-1, ]) - log(x[-n, ])
fit = sparse_t_copula(x = y[, 1:10], nu = 3, lambda = 0.5)
