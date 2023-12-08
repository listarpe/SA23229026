## -----------------------------------------------------------------------------
adj_matrix <- matrix(c(1, 0, 0, 0, 
                      0, 1, 0, 0, 
                      0, 0, 1, 0, 
                      1, 1, 1, 1), 
                    nrow = 4, byrow=TRUE)

## -----------------------------------------------------------------------------
library(SA23229026)
pagerankr(adj_matrix)
pagerankc(adj_matrix)

## -----------------------------------------------------------------------------
library(microbenchmark)
ts <- microbenchmark(pagerankR=pagerankr(adj_matrix), pagerankc=pagerankc(adj_matrix))
summary(ts)[,c(1,3,5,6)]

## -----------------------------------------------------------------------------
set.seed(123)
data <- matrix(c(rnorm(30, 0, 0.2), rnorm(30, 2, 0.2), rnorm(30, 2, 0.2), 
                 rnorm(30, 0, 0.2), rnorm(30, 0, 0.2), rnorm(30, 2, 0.2)), ncol = 2)
plot(data)

## -----------------------------------------------------------------------------
result <- kmeans(data, 3)
result

