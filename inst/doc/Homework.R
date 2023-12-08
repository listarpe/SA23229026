## -----------------------------------------------------------------------------
data <- iris
summary(data)

## ----plot---------------------------------------------------------------------
# par(pin=c(3,2))
plot(data$Sepal.Length, data$Sepal.Width, col=data$Species, xlab = "Sepal Length(cm)", ylab = "Sepal Width(cm)", main="Three Classes of Iris")

## ----echo=FALSE---------------------------------------------------------------
height <- c(151, 174, 138, 186, 128, 136, 179, 163, 152, 131)
weight <- c(63, 81, 56, 91, 47, 57, 76, 72, 62, 48)
data = data.frame(Height = height, Weight = weight)
knitr::kable(data, format = "markdown", align = "c")

## ----echo=FALSE---------------------------------------------------------------
regression <- lm(height~weight)
plot(weight,height, main = "Height & Weight Regression", xlab = "Weight(Kg)",ylab = "Height(cm)")
abline(regression, col = "blue")

## ----echo=FALSE---------------------------------------------------------------
set.seed(123)
data <- data.frame(
  x = c(rnorm(50, mean = 0), rnorm(50, mean = 4)),
  y = c(rnorm(50, mean = 0), rnorm(50, mean = 4))
)
plot(data, main = "100 Samples")

## ----echo=FALSE---------------------------------------------------------------
kmeans_result <- kmeans(data, centers = 2)
plot(data, col = kmeans_result$cluster, main = "Result of K-means")
points(kmeans_result$centers, col = 1:3, pch = 8, cex = 2)

## -----------------------------------------------------------------------------
Mysample <- function(x, size, prob=NULL){
  n <- length(x) # 获取区间x中的元素个数
  if (is.null(prob)) {
    p <- rep(1/n, n) # 如果prob为空，所有样本概率相同
  }
  cp <- cumsum(prob)
  U <- runif(size)
  r <- x[findInterval(U, cp)+1]
  return(r)
}

## -----------------------------------------------------------------------------
Mysample(0:1, 10)

## -----------------------------------------------------------------------------
r <- Mysample(c(2, 5, 8, 9), 1000, prob = c(.2, .3, .4, .1))
barplot(table(r))

## -----------------------------------------------------------------------------
n <- 1e4
u <- runif(n)
r = sign(0.5-u)*log(1-abs(2*u-1))
hist(r, prob=TRUE, ylim=c(0, 0.55), main=expression(f(x)==frac(1,2)*e^-abs(x)))
y <- seq(-5, 5, .01)
lines(y, 0.5*exp(-abs(y)))

## -----------------------------------------------------------------------------
myBeta <- function(a, b, n){
  t = (a-1)/(a+b-2)
  c = 1/beta(a, b)*t**(a-1)*(1-t)**(b-1)
  k <- 0
  j <- 0
  y <- numeric(n)
  while(k<n){
    u1 <- runif(1)
    u2 <- runif(1)
    j <- j + 1
    if (dbeta(u2, a, b) > u1){
      k <- k+1
      y[k] <- u2
    }
  }
  out <- list(y = y, j = j, c = c)
  return(out)
}

## -----------------------------------------------------------------------------
out = myBeta(3, 2, 1000)
cat("iterations: ", out$j)

## -----------------------------------------------------------------------------
x <- out$y
hist(x, prob=TRUE, ylim=c(0, out$c), main="Beta(2,3)")
y <- seq(0, 1, .001)
lines(y, dbeta(y, 3, 2))

## -----------------------------------------------------------------------------
myfe <- function(n){
  k <- 0
  y <- numeric(n)
  while(k<n){
    u1 <- runif(1, -1, 1)
    u2 <- runif(1, -1, 1)
    u3 <- runif(1, -1, 1)
    if ((abs(u3)>abs(u2) || abs(u3)== abs(u2)) && (abs(u3)>abs(u1) || abs(u3)== abs(u1))){
      k <- k + 1
      y[k] <- u2
    } else {
      k <- k + 1
      y[k] <- u3
    }
  }
  return(y)
}

## -----------------------------------------------------------------------------
y <- myfe(1e5)
hist(y, prob = TRUE, main=expression(f[e](x)==frac(3,4)(1-x^2)))
x<- seq(-1, 1, .01)
lines(x, 0.75 * (1 - x**2))

## -----------------------------------------------------------------------------
n <- 10^6 # 试验次数
K <- 100 # 模拟次数
rhos <- c(0.5, 0.8, 1) # rho的取值
for (rho in rhos){
  m <- rbinom(K, n, 2*rho/pi) # 以2*rho/pi的概率生成随机数
  pi_hat <- 2*rho*n/m # 计算pi_hat
  D = var(pi_hat) # 计算方差
  cat("rho =", rho,"    variance =", D, "\n")
}

## -----------------------------------------------------------------------------
K = 10^6 # 模拟次数
n = 100 # 采样次数
MC1 = numeric(K)
MC2 = numeric(K)
for (i in 1:K){
  u = runif(n/2)
  v = runif(n/2)
  MC1[i] = mean(exp(c(u,v)))
  MC2[i] = mean(exp(c(u, 1-u)))
}
theta1 = mean(MC1)
theta2 = mean(MC2)
var1 = var(MC1)
var2 = var(MC2)
cat("theta1:",theta1,"\ntheta2:",theta2,"\nvariance reduction:", 1-var2/var1)

## -----------------------------------------------------------------------------
m <- 1e4
theta.hat <- var.hat <- numeric(2)
g <- function(x){
  x^2/sqrt(2*pi)*exp(-x^2/2)*(x>1)
}

x <- rnorm(m) # f1
fg <- g(x) / dnorm(x)
theta.hat[1] <- mean(fg)
var.hat[1] <- var(fg)

x <- rexp(m, 1) # f2
fg <- g(x) * exp(x)
theta.hat[2] <- mean(fg)
var.hat[2] <- var(fg)

rbind(theta.hat, var.hat)

## -----------------------------------------------------------------------------
m <- 1e4
g <- function(x){
  x^2/sqrt(2*pi)*exp(-x^2/2)*(x>1)
}

x <- rnorm(m, mean=1.5) # 设置正太分布的均值为1.5
fg <- g(x) / dnorm(x, mean=1.5)
theta.hat <- mean(fg)
se <- sd(fg)

rbind(theta.hat, se)

## -----------------------------------------------------------------------------
M <- 10000 # 总采样次数
k <- 5 # 层数
E <- numeric(k) # 保存每次中每层的均值
N = 50 # 估计次数
r <- M/k/N
est = numeric(N) # 保存每次估计的均值

g<-function(x){
  exp(-x)/(1+x^2)*(x>0)*(x<1)
}

for(i in 1:N){
  for(j in 1:k){
    E[j]<-mean(g(runif(r,(j-1)/k,j/k)))
  }
  est[i] = mean(E)
}
theta.hat = mean(est)
se = sd(est)

rbind(theta.hat, se)

## -----------------------------------------------------------------------------
n <- 20
alpha <- .05
calcCI <- function(n, alpha){
  y <- rchisq(n, df = 2)
  return(c(mean(y)-qt(1-alpha/2, df=n-1)*sd(y)/sqrt(n), mean(y)+qt(1-alpha/2, df=n-1)*sd(y)/sqrt(n)))
}

CL <- replicate(1000, expr = calcCI(n=n, alpha=alpha))
rate = mean((CL[1,]<2)*(CL[2,]>2))
print(rate)

## -----------------------------------------------------------------------------
n <- 20
alpha <- .05
mu0 <- 1

m <- 10000
p <- matrix(0, 3, m)

for (j in 1:m) {
  x <- rchisq(n, df=1) # X^2(1)
  ttest <- t.test(x, alternative = "two.sided", mu = mu0)
  p[1, ] <- ttest$p.value

  x <- runif(n, 0, 2) # U(0, 2)
  ttest <- t.test(x, alternative = "two.sided", mu = mu0)
  p[2, ] <- ttest$p.value

  x <- rexp(n, 1) # Exponential(1)
  ttest <- t.test(x, alternative = "two.sided", mu = mu0)
  p[3, ] <- ttest$p.value
}

p.hat <- apply(p<alpha, 1, mean)
se.hat <- sqrt(p.hat * (1 - p.hat) / m)
table = data.frame(
  p.hat = p.hat,
  se.hat = se.hat,
  row.names=c("X^2(1)", "U(0, 2)", "Exponential(1)")
)
print(table)

## -----------------------------------------------------------------------------
M = 1000
m = 1000
alpha = 0.1

result_bfn <- result_bh <- matrix(0,M,3)

for (i in 1:M){
  set.seed(1000+i)
  # 对立假设成立为1， 原假设成立为0
  # 拒绝原假设为1，不拒绝原假设为0
  label_true <- c(rep(FALSE, 0.95*m), rep(TRUE, 0.05*m))
  p = c(runif(0.95*m, 0, 1), rbeta(0.05*m, 0.1, 1))

  # 校正
  p.bfn = p.adjust(p,method='bonferroni')
  p.bh = p.adjust(p,method='BH')
  # 获取检验结果
  label_bfn <- p.bfn < alpha
  label_bh <- p.bh < alpha
  # 构造混淆矩阵
  confusion_bfn <- table(Predicted = label_bfn, Actual = label_true)
  confusion_bh <- table(Predicted = label_bh, Actual = label_true)
  result_bfn[i, ] <- c(confusion_bfn[2,1]>0, confusion_bfn[2,1]/sum(confusion_bfn[2, ]), confusion_bh[2,2]/sum(confusion_bfn[, 2]))
  result_bh[i, ] <- c(confusion_bh[2,1]>0, confusion_bh[2,1]/sum(confusion_bh[2, ]), confusion_bh[2,2]/sum(confusion_bh[, 2]))
}

table = data.frame(
  Bonferroni = colMeans(result_bfn),
  B.H = colMeans(result_bh),
  row.names = c("FWER", "FDR", "TPR")
)

print(t(table))

## -----------------------------------------------------------------------------
library(boot)
sample_size = c(5, 10, 20) # 采样个数
B = 1000 # bootstrap重复次数
m = 1000 # 模拟次数

# 计算lambda的估计值
b.lambda <- function(x, i){
  1/mean(x[i])
}

total = matrix(0, 3, 2) # 记录bootstrap估计值
for (i in 1:3){
  n = sample_size[i]
  matrix_temp = matrix(0, 1000, 2)
  for (m in 1:1000){
    x = rexp(n, 2)
    obj <- boot(data=x, statistic = b.lambda, R=1000)
    matrix_temp[m, ] <- c(mean(obj$t)-obj$t0, sd(obj$t))
  }
  total[i, ] <- colMeans(matrix_temp)
}

table = data.frame(
  bias.bootstrap = total[, 1],
  bias.theoretical = 2/(sample_size-1),
  se.bootstrap = total[, 2],
  se.theoretical = 2*sample_size/((sample_size-1)*sqrt(sample_size-2)),
  row.names = sample_size
)
print(table)

## -----------------------------------------------------------------------------
library(bootstrap)
library(boot)

boot.t.ci <- function(x, B = 500, R = 100, level = .95, statistic){
  # 计算bootstrap的t置信区间
  x <- as.matrix(x); n <- nrow(x)
  stat <- numeric(B); se <- numeric(B)
  boot.se <- function(x, R, f) {
    #计算bootstrap的估计标准差，f是统计量
    x <- as.matrix(x); m <- nrow(x)
    th <- replicate(R, expr = {
    i <- sample(1:m, size = m, replace = TRUE)
    f(x[i, ])
    })
    return(sd(th))
  }
  for (b in 1:B) {
    j <- sample(1:n, size = n, replace = TRUE)
    y <- x[j, ]
    stat[b] <- statistic(y)
    se[b] <- boot.se(y, R = R, f = statistic)
  }
  stat0 <- statistic(x)
  t.stats <- (stat - stat0) / se
  se0 <- sd(stat)
  alpha <- 1 - level
  Qt <- quantile(t.stats, c(alpha/2, 1-alpha/2), type = 1)
  names(Qt) <- rev(names(Qt))
  CI <- rev(stat0 - Qt * se0)
}

dat <- law # 导入数据

stat <- function(dat) {
  # 计算相关性的统计量
  cor(dat[,1], dat[,2])
}

ci <- boot.t.ci(dat, statistic = stat, B = 2000, R = 200)

print(ci)

## -----------------------------------------------------------------------------
library(boot)
x = c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487)
b.mean <- function(x,i){
  mean(x[i])
}

obj <- boot(data=x, statistic = b.mean, R=1000)
ci <- boot.ci(obj, type=c("norm","basic","perc","bca"), conf=0.95)
ci

## -----------------------------------------------------------------------------
library(bootstrap)
x <- scor

cov_matrix <- cov(scor, scor)
lambda <- eigen(cov_matrix)$values
theta <- lambda[1]/sum(lambda)

n <- nrow(x)
theta.hat <- numeric(n)
for (i in 1:n){
  cov.hat <- cov(x[-i,], x[-i,])
  lambda.hat <- eigen(cov.hat)$values
  theta.hat[i] <- lambda.hat[1]/sum(lambda.hat)
}

bias.jack = (n-1)*(mean(theta.hat)-theta)
se.jack = sqrt((n-1)^2/n)*sd(theta.hat)

c(bias.jack = bias.jack, se.jack = se.jack)

## -----------------------------------------------------------------------------
library(DAAG); attach(ironslag)

n <- length(magnetic) #in DAAG ironslag
e1 <- e2 <- e3 <- e4 <- numeric(n*(n-1)/2)

# for n-fold cross validation
# fit models on leave-one-out samples
k = 0
for (i in 1:n-1) {
  for (j in (i+1):n) {
    y <- magnetic[-c(i, j)]
    x <- chemical[-c(i, j)]

    J1 <- lm(y ~ x)
    yhat1 <- J1$coef[1] + J1$coef[2] * chemical[c(i, j)]
    e1[k] <- mean((magnetic[c(i, j)] - yhat1)^2)

    J2 <- lm(y ~ x + I(x^2))
    yhat2 <- J2$coef[1] + J2$coef[2] * chemical[c(i, j)] +
    J2$coef[3] * chemical[c(i, j)]^2
    e2[k] <- mean((magnetic[c(i, j)] - yhat2)^2)

    J3 <- lm(log(y) ~ x)
    logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[c(i, j)]
    yhat3 <- exp(logyhat3)
    e3[k] <- mean((magnetic[c(i, j)] - yhat3)^2)

    J4 <- lm(log(y) ~ log(x))
    logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[c(i, j)])
    yhat4 <- exp(logyhat4)
    e4[k] <- mean((magnetic[c(i, j)] - yhat4)^2)

    k = k + 1
  }
}

table = data.frame(
  J1 = mean(e1),
  J2 = mean(e2),
  J3 = mean(e3),
  J4 = mean(e4),
  row.names = "MSE"
)

print(table)

## ----warning=FALSE------------------------------------------------------------
rm(list = ls())
attach(chickwts)
x <- sort(as.vector(weight[feed == "soybean"]))
y <- sort(as.vector(weight[feed == "linseed"]))
detach(chickwts)

w <- function(x, y){
  ecdf1 <- ecdf(x)
  ecdf2 <- ecdf(y)
  n <- length(x)
  m <- length(y)
  sum1 = sum((ecdf1(x)-ecdf2(x))^2)
  sum2 = sum((ecdf1(y)-ecdf2(y))^2)
  return(m*n/(m+n)^2*(sum1+sum2))
}

R <- 999 #number of replicates
z <- c(x, y) #pooled sample
K <- 1:26
D <- numeric(R) #storage for replicates
D0 <- w(x, y)
for (i in 1:R) {
  #generate indices k for the first sample
  k <- sample(K, size = 14, replace = FALSE)
  x1 <- z[k]
  y1 <- z[-k] #complement of x1
  D[i] <- w(x1, y1)
}
p <- mean(c(D0, D) >= D0)
p

## -----------------------------------------------------------------------------
n1 <- 20
n2 <- 30
mu1 <- mu2 <- 0
sigma1 <- sigma2 <- 1
x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)

count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  return(max(c(outx, outy)))
}

R <- 999 #number of replicates
z <- c(x, y) #pooled sample
K <- 1:(n1+n2)
D <- numeric(R) #storage for replicates
D0 <- count5test(x, y)
for (i in 1:R) {
  #generate indices k for the first sample
  k <- sample(K, size = n1, replace = FALSE)
  x1 <- z[k]
  y1 <- z[-k] #complement of x1
  D[i] <- count5test(x1, y1)
}
p <- mean(c(D0, D) >= D0)
p

## -----------------------------------------------------------------------------
logistica <- function(N, b1, b2, b3, f0){ 
  x1 <- rpois(N, 1) 
  x2 <- rexp(N, 1)
  x3 <- rbeta(N, 1, 0.5)
  gg <- function(alpha){ # 内置函数
    t <- exp(-alpha-b1*x1-b2*x2-b3*x3)
    p <- 1/(1+t)
    mean(p)-f0
  }
  result <- uniroot(gg, c(-20, 0)) # 使用uniroot求方程的根
  return(result)
}

## -----------------------------------------------------------------------------
N = 1e6
b1 = 0
b2 = 1
b3 = -1
f0 = c(0.1, 0.01, 0.001, 0.0001)
a = numeric(4) # 记录a的值
for (i in 1:4){
  a[i] = logistica(N, b1, b2, b3, f0[i])$root
}

table = data.frame(f0=f0, a=a)
print(round(table, 4))

## -----------------------------------------------------------------------------
plot(-log(f0), a, type="l")

## -----------------------------------------------------------------------------
# laplace分布概率密度函数
laplace <- function(x){
  0.5*exp(-abs(x))
}

# random walk metropolis 采样
rw.Metropolis <- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0  # markov chain的第一个值
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (laplace(y) / laplace(x[i-1]))){ # 接受
      x[i] <- y
    } else { # 拒绝
      x[i] <- x[i-1]
      k <- k + 1
    }
  }
  return(list(x=x, k=k))
}

N <- 2000 # 迭代步数
sigma <- c(.05, .5, 2, 16) # 方差
x0 <- 20 # markov chain的第一个值
rw1 <- rw.Metropolis(sigma[1], x0, N)
rw2 <- rw.Metropolis(sigma[2], x0, N)
rw3 <- rw.Metropolis(sigma[3], x0, N)
rw4 <- rw.Metropolis(sigma[4], x0, N)

index = 1:2000
par(mar = c(3, 3, 3, 3))
par(mfrow=c(2,2))
plot(index, rw1$x, type="l", xlab="sigma=.05", ylab="x")
plot(index, rw2$x, type="l", xlab="sigma=.5", ylab="x")
plot(index, rw3$x, type="l", xlab="sigma=2", ylab="x")
plot(index, rw4$x, type="l", xlab="sigma=16", ylab="x")

## -----------------------------------------------------------------------------
accptance = numeric(4)
accptance[1] <- 1-rw1$k/N
accptance[2] <- 1-rw2$k/N
accptance[3] <- 1-rw3$k/N
accptance[4] <- 1-rw4$k/N
table.accptance = data.frame(sigma = sigma, apptance.rate = accptance)
table.accptance

## -----------------------------------------------------------------------------
# 初始化参数
N <- 5000 # markov chain的长度
burn <- 1000 # burn-in 长度
X <- matrix(0, N, 2) # 记录markov chain
rho <- 0.9 # 相关系数
mu1 <- 0 # 均值
mu2 <- 0
sigma1 <- 1 # 方差
sigma2 <- 1
s1 <- sqrt(1-rho^2)*sigma1 # 边缘分布的方差
s2 <- sqrt(1-rho^2)*sigma2

X[1, ] <- c(mu1, mu2) # 初始化
for (i in 2:N) {
  x2 <- X[i-1, 2]
  m1 <- mu1 + rho * (x2 - mu2) * sigma1/sigma2
  X[i, 1] <- rnorm(1, m1, s1)
  x1 <- X[i, 1]
  m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1
  X[i, 2] <- rnorm(1, m2, s2)
}
b <- burn + 1 
x <- X[b:N, ] # 去除burn-in部分的markov chain

## -----------------------------------------------------------------------------
plot(x, main="", cex=.5, xlab=bquote(X), ylab=bquote(Y), ylim=range(x[,2]))

## -----------------------------------------------------------------------------
model <- lm(x[, 2]~x[, 1])
summary(model)

## -----------------------------------------------------------------------------
shapiro.test(model$residuals)

## -----------------------------------------------------------------------------
se <- sd(model$residuals)
se

## -----------------------------------------------------------------------------
# 目标分布
f <- function(x, sigma) {
  if (any(x < 0)) { # 定义域x大于0
    return (0)
  }
  stopifnot(sigma > 0) # sigma需要大于0
  return((x / sigma^2) * exp(-x^2 / (2*sigma^2)))
}

# Gelman.Rubin方法
Gelman.Rubin <- function(psi) {
  psi <- as.matrix(psi) # psi是markov chain的统计量
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi) #行均值
  B <- n * var(psi.means) #between variance est.
  psi.w <- apply(psi, 1, "var") #within variances
  W <- mean(psi.w) #within est.
  v.hat <- W*(n-1)/n + (B/n) #方差的上界
  r.hat <- v.hat / W #G-R 统计量
  return(r.hat)
}

# 参数
sigma <- .2 # 提议分布的参数
k <- 4 #chain的数量
b <- 1000 #burn-in长度

#生成初始值
set.seed(123)
x0 <- rchisq(k, df=1)

X <- matrix(x0) # 记录markov chain
psi <- matrix(x0) # 记录统计量
R.hat <- 2 # 初始化R.hat，使其大于1.2

i = 2 # chain的长度
while (R.hat >= 1.2) { # 循环直到R.hat大于1.2，实现了控制
  xt <- X[, i-1]
  y <- rchisq(k, df = xt)
  r1 <- f(y, sigma) * dchisq(xt, df = y)
  r2 <- f(xt, sigma) * dchisq(y, df = xt)
  flag <- runif(4) <= r1/r2 # 这里flag用来表示时候接受提议分布
  X <- cbind(X, flag*y+(1-flag)*xt) # 将拒绝或接受的值加入
  psi <- cbind(psi, psi[, i-1]+(X[, i]-psi[, i-1])/i) # 更新统计量
  R.hat <- Gelman.Rubin(psi)
  i = i + 1
}

print(R.hat)

n = ncol(psi)

#画chain
par(mar = c(3, 3, 3, 3))
par(mfrow=c(2,2))
for (i in 1:k)
  plot(psi[i, (b+1):n], type="l", xlab=i, ylab=bquote(psi))
par(mfrow=c(1,1)) #restore default

#画R-hat statistics
rhat <- rep(0, n)
for (j in (b+1):n)
  rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
abline(h=1.1, lty=2)


## -----------------------------------------------------------------------------
library(coda)

mcmclist <- mcmc.list(mcmc(X[1, ]), mcmc(X[2, ]), mcmc(X[3, ]), mcmc(X[4, ]))
gelman.diag(mcmclist)

## -----------------------------------------------------------------------------
gelman.plot(mcmclist)

## -----------------------------------------------------------------------------
u <- c(11, 8, 27, 13, 16, 0, 23, 10, 24, 2)
v <- c(12, 9, 28, 14, 17, 1, 24, 11, 25, 3)
mlogL <- function(lambda) {
    # minus log-likelihood
    return(sum(log(exp(-lambda*u)-exp(-lambda*v))))
}

res <- optimize(mlogL, lower=0, upper=10, maximum=TRUE)
res$maximum

## -----------------------------------------------------------------------------
u <- c(11, 8, 27, 13, 16, 0, 23, 10, 24, 2)
v <- c(12, 9, 28, 14, 17, 1, 24, 11, 25, 3)

# EM算法的迭代函数
f <- function(lambda) {
    1/mean(((-v-1/lambda)*exp(-lambda*v)+(u+1/lambda)*exp(-lambda*u))/(exp(-lambda*u)-exp(-lambda*v)))
}

loss <- 1 # 前后两次迭代的差
t <- 1 # 迭代次数

lambda <- 0.1
while (t<1000 & loss>1e-8) { # 迭代次数大于1000或误差小于1e-8时停止迭代
    lambda1 <- lambda
    lambda <- f(lambda)
    loss <- abs(lambda - lambda1)
    t = t + 1
}
lambda

## -----------------------------------------------------------------------------
solve.game <- function(A) {
    #solve the two player zero-sum game by simplex method
    #optimize for player 1, then player 2
    #maximize v subject to ...
    #let x strategies 1:m, and put v as extra variable
    #A1, the <= constraints
    #
    min.A <- min(A)
    A <- A - min.A #so that v >= 0
    max.A <- max(A)
    A <- A / max(A)
    m <- nrow(A)
    n <- ncol(A)
    it <- n^3
    a <- c(rep(0, m), 1) #objective function
    A1 <- -cbind(t(A), rep(-1, n)) #constraints <=
    b1 <- rep(0, n)
    A3 <- t(as.matrix(c(rep(1, m), 0))) #constraints sum(x)=1
    b3 <- 1
    sx <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
        maxi=TRUE, n.iter=it)
    #the ’solution’ is [x1,x2,...,xm | value of game]
    #
    #minimize v subject to ...
    #let y strategies 1:n, with v as extra variable
    a <- c(rep(0, n), 1) #objective function
    A1 <- cbind(A, rep(-1, m)) #constraints <=
    b1 <- rep(0, m)
    A3 <- t(as.matrix(c(rep(1, n), 0))) #constraints sum(y)=1
    b3 <- 1
    sy <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
        maxi=FALSE, n.iter=it)
    soln <- list("A" = A * max.A + min.A,
        "x" = sx$soln[1:m],
        "y" = sy$soln[1:n],
        "v" = sx$soln[m+1] * max.A + min.A)
    soln
}

## -----------------------------------------------------------------------------
#enter the payoff matrix
A <- matrix(c(  0,-2,-2,3,0,0,4,0,0,
                2,0,0,0,-3,-3,4,0,0,
                2,0,0,3,0,0,0,-4,-4,
                -3,0,-3,0,4,0,0,5,0,
                0,3,0,-4,0,-4,0,5,0,
                0,3,0,0,4,0,-5,0,-5,
                -4,-4,0,0,0,5,0,0,6,
                0,0,4,-5,-5,0,0,0,6,
                0,0,4,0,0,5,-6,-6,0), 9, 9)
library(boot) #needed for simplex function
s <- solve.game(A)

## -----------------------------------------------------------------------------
round(cbind(s$x, s$y), 7)

## -----------------------------------------------------------------------------
B <- A + 2
s1 <- solve.game(B)
round(cbind(s1$x, s1$y), 7)

## -----------------------------------------------------------------------------
cat("A: ", s$v)
cat("B: ", s1$v)

## -----------------------------------------------------------------------------
x <- list(a = 1, b = list(c = 2, d = 3))
unlist(x)
as.vector(x)

## -----------------------------------------------------------------------------
a <- c(1, 2, 3)
dim(a)

## -----------------------------------------------------------------------------
x <- matrix(data = 1:6, nrow = 2, ncol = 3)
is.matrix(x)
is.array(x)

## -----------------------------------------------------------------------------
x <- data.frame(
    numeric_col = c(1, 2, 3),
    character_col = c("a", "b", "c"),
    logical_col = c(TRUE, FALSE, TRUE)
)

as.matrix(x)

## -----------------------------------------------------------------------------
x1 <- data.frame(matrix(nrow = 0, ncol = 3))
x1

## -----------------------------------------------------------------------------
x2 <- data.frame(matrix(nrow = 3, ncol = 0))
x2

## -----------------------------------------------------------------------------
scale01 <- function(x) {
    rng <- range(x, na.rm = TRUE)
    (x - rng[1]) / (rng[2] - rng[1])
}

x <- data.frame(
    a = 1:3,
    b = 7:9,
    c = 4:6
)

lapply(x, scale01)

## -----------------------------------------------------------------------------
y <- data.frame(
    a = 1:3,
    b = c('q', 'w', 'e'),
    c = 4:6
)

lapply(y, function(x) {if (is.numeric(x)) scale01(x) else x})

## -----------------------------------------------------------------------------
x <- data.frame(
    a = 1:3,
    b = 7:9,
    c = 4:6
)

vapply(x, sd, numeric(1))

## -----------------------------------------------------------------------------
y <- data.frame(
    a = 1:3,
    b = c('q', 'w', 'e'),
    c = 4:6
)

vapply(y[vapply(y, is.numeric, logical(1))], sd, numeric(1))

## -----------------------------------------------------------------------------
#initialize constants and parameters
N <- 50000 #length of chain
burn <- 1000 #burn-in length

n <- 10
a <- 2
b <- 3

gibbsR <- function(N, burn, n, a, b){
    x1 <- 3
    y1 <- 0.5
    X <- matrix(0, N, 2) #the chain, a bivariate sample
    X[1, ] <- c(x1, y1)

    for (i in 2:N) {
        yt <- X[i-1, 2]
        X[i, 1] <- rbinom(1, n, yt)
        xt <- X[i, 1]
        X[i, 2] <- rbeta(1, xt + a, n - xt + b)
    }
    X
}

## -----------------------------------------------------------------------------
library(Rcpp)
library(microbenchmark)

# dir_cpp <- "./"
# sourceCpp(paste0(dir_cpp, "gibbsC.cpp"))

cppFunction('NumericMatrix gibbsC(int N, int burn, int n, int a, int b) {
    int x1 = 3;
    double y1 = 0.5;
    NumericMatrix X(N, 2);
    X(0, 0) = x1;
    X(0, 1) = y1;

    int xt = x1;
    double yt = y1;
    for (int i = 1; i < N; i++) {
        yt = X(i-1, 1);
        X(i, 0) = R::rbinom(n, yt);
        xt = X(i, 0);
        X(i, 1) = R::rbeta(xt + a, n - xt + b);
    }

    return X;
}')

ts <- microbenchmark(gibbsR=gibbsR(N, burn, n, a, b), gibbsC=gibbsC(N, burn, n, a, b))
summary(ts)[,c(1,3,5,6)]

