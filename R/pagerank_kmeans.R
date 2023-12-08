#' @import Rcpp
#' @import DAAG
#' @import boot
#' @import bootstrap
#' @import coda
#' @import microbenchmark
#' @useDynLib SA23229026
NULL


#' @title Pagerank using R
#' @description Pagerank using R
#' @param adj_matrix adjacency matrix
#' @param alpha Correction probability of transition probability matrix
#' @param max_iter Maximum Number Of Iterations
#' @return Pagerank values for each node
#' @examples
#' \dontrun{
#' r <- pagerankr(adj_matrix, alpha=0.9, max_iter=5000)
#' }
#' @export
pagerankr <- function(adj_matrix, alpha=0.85, max_iter=10000){
  divide_by_row_sum <- function(row) {
    return(row / sum(row))
  }
  p <- t(apply(adj_matrix, 1, divide_by_row_sum))
  n <- nrow(p)
  p <- alpha * p + (1-alpha)/n * matrix(1, nrow = n, ncol = n)
  
  r <- matrix(1/n, nrow = 1, ncol = n)
  err <- 1
  i <- 0
  while (err > 1e-6 & i < max_iter) {
    i <-  i + 1
    r0 <- r
    r <- r0 %*% adj_matrix
    r <- r / sum(r)
    err <- norm(r-r0, type = "F")/norm(r, type = "F")
  }
  return(r)
}


#' @title kmeans using R
#' @description kmeans using R
#' @param data Sample points
#' @param k Number of clusters
#' @param max_iter Maximum Number Of Iterations
#' @return Center and cluster of each sample
#' @examples
#' \dontrun{
#' kmeans(data, k=3, max_iter=5000)
#' }
#' @export
kmeans <- function(data, k, max_iter=10000) {
  centroids <- data[sample(1:nrow(data), k), ]
  distances <- matrix(0, nrow = nrow(data), ncol = k)
  i <- 0
  while (i < max_iter) {
    for (i in 1:k) {
      distances[, i] <- sqrt(rowSums((t(t(data)-centroids[i, ]))^2))
    }
    clusters <- apply(distances, 1, which.min)
    new_centroids <- matrix(0, nrow = k, ncol = 2)
    for (i in 1:k) {
      new_centroids[i, ] <- colMeans(data[clusters == i, ])
    }
    if (all(centroids == new_centroids)) {
      break
    }
    # 更新聚类中心并增加迭代次数
    centroids <- new_centroids
    i <- i + 1
  }
  list(centroids=centroids, clusters=clusters)
}