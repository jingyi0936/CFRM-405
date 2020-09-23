f <- function(x)
  c(2*x[1]*x[5] + 2*x[1]*x[6] + 1,
    2*x[2]*x[5] + 4*x[2]*x[6] - 2,
    6*x[3]*x[6] + 3,
    -2*x[4]*x[5] - 4,
    x[1]^2 + x[2]^2 - x[4]^2 - 1,
    x[1]^2 + 2*x[2]^2 + 3*x[3]^2 - 6)

df <- function(x){
  grad <- matrix(0, 6, 6)
  grad[1,] <- c(2*x[5] + 2*x[6], 0, 0, 0, 2*x[1], 2*x[1])
  grad[2,] <- c(0, 2*x[5] + 4*x[6], 0, 0, 2*x[2], 4*x[2])
  grad[3,] <- c(0, 0, 6*x[6], 0, 0, 6*x[3])
  grad[4,] <- c(0, 0, 0, -2*x[5], -2*x[4], 0)
  grad[5,] <- c(2*x[1], 2*x[2], 0, -2*x[4], 0, 0)
  grad[6,] <- c(2*x[1], 4*x[2], 6*x[3], 0, 0, 0)
  grad
}

x <- c(1, -1, 1, -1, 1, -1)
u <- rep(1, 6)

while(sqrt(sum(u^2)) / sqrt(sum(x^2)) > 1e-6){
  u <- solve(df(x), f(x))
             x <- (x - u)
}

round(df(x), digits = 3)

# section 5.7 1
# five year semiannual with 3.375% coupon rate and price 100 1/32
# t_cash_flow = cash flow dates
# v_vash_flow = cash flows
# B = 132/100
B = 132/100
t_cash_flow = c(6/12, 12/12, 18/12, 24/12, 30/12, 36/12, 42/12, 48/12, 54/12, 60/12)
v_cash_flow = c(rep(3.375/2, 9), 100+100*3.375/2)
yield <- function(old)
  (sum(v_cash_flow *exp(-old*t_cash_flow)) - B)/(sum((-1)*t_cash_flow*v_cash_flow*
                                                       exp(-old*t_cash_flow)))
old = 0.1
new = old - 1
while(abs(new)/abs(old) > 1e-6){
  old = new
  new <- old +  yield(old)  
}


# correlation matrix
cor <- matrix(c(1, -0.3, 0.4, 0.25, -0.2, 
                -0.3, 1, -0.1, -0.2, 0.15,
                0.4, -0.1, 1, 0.35, 0.25, 
                0.25, -0.2, 0.35, 1, -0.15,
                -0.2, 0.15, 0.25, -0.15, 1), nrow = 5, ncol = 5)
mu <-c(0.08, 0.10, 0.13, 0.15, 0.20)
sd <- c(0.14, 0.18, 0.23, 0.25, 0.35)

# covariance matrix
c11 = (sd[1])^2
c12 = cor[1,2]*sd[1]*sd[2]
c13 = cor[1,3]*sd[1]*sd[3]
c14 = cor[1,4]*sd[1]*sd[4]
c15 = cor[1,5]*sd[1]*sd[5]
c22 = (sd[2])^2
c23 = cor[2,3]*sd[2]*sd[3]
c24 = cor[2,4]*sd[2]*sd[4]
c25 = cor[2,5]*sd[2]*sd[5]
c33 = (sd[3])^2
c34 = cor[3,4]*sd[3]*sd[4]
c35 = cor[3,5]*sd[3]*sd[5]
c44 = (sd[4])^2
c45 = cor[4,5]*sd[4]*sd[5]
c55 = (sd[5])^2

cov <- matrix(c(c11, c12, c13, c14, c15,
                c12, c22, c23, c24, c25,
                c13, c23, c33, c34, c35,
                c14, c24, c34, c44, c45,
                c15, c25, c35, c45, c55), nrow = 5, ncol = 5)
# 3(i)
G1 <- function(x, mu, cov, mup){
  n <- length(mu)
  c(2*t(x[1:n])*cov + rep(x[n+1], n) + 2*x[n+2]*mu,
  sum(x[1:n]) - 1,
  t(mu) %*% x[1:n] - mup)
}

DG1 <- function(x, mu, cov, mup){
  n <- length(mu)
  grad <- matrix(0.0, n+2, n+2)
  grad[1:n, 1:n] <- 2*cov
  grad[1:n, n+1] <- 1
  grad[1:n, n+2] <- (mu)
  grad[n+1, 1:n] <- 1
  grad[n+2, 1:n] <- t(mu)
  grad
}
x1 <- c(rep(0.5, 5), 1, 1)
u1 <- rep(1, length(x1))
while(sqrt(sum(u1^2)) / sqrt(sum(x1^2)) > 1e-6){
  u1 <- solve(DG(x1, mu, cov, 0.15), G(x1, mu, cov, 0.15))
  x1 <- x1-u1
}


# 3(ii)
G <- function(x, mu, cov, sigmap2){
  n <- length(mu)
  c(mu + rep(x[n+1], n) + 2*x[n+2]*(cov %*% x[1:n]),
  sum(x[1:n]) - 1,
  t(x[1:n] %*% cov %*% x[1:n] - sigmap2))
}
DG <- function(x, mu, cov, sigmap2){
  n <- length(mu)
  grad <- matrix(0.0, n+2, n+2)
  grad[1:n, 1:n] <- 2*x[n+2]*cov
  grad[1:n, n+1] <- 1
  grad[1:n, n+2] <- 2*(cov %*% x[1:n])
  grad[n+1, 1:n] <- 1
  grad[n+2, 1:n] <- 2*t(x[1:n]) %*% cov
  grad
}
x <- c(rep(0.5, 5), 1, 1)
u <- rep(1, length(x))
while(sqrt(sum(u^2)) / sqrt(sum(x^2)) > 1e-6){
  u <- solve(DG(x, mu, cov, 0.25^2), G(x, mu, cov, 0.25^2))
  x <- x-u
}

# 4 R function
bsc <- function(S, Tim, t, K, r, s, q){
  d1 <- (log(S/K) + (r - q + 0.5*s^2)*(Tim - t))/(s*sqrt(Tim - t))
  d2 <- d1 - s*sqrt(Tim - t)
  S*exp(-q*(Tim - t))*pnorm(d1) - K*exp(-r*(Tim - t))*pnorm(d2)
}

bisection <- function(f, a, b, tol){
  while(b-a > tol){
    c <- (a+b)/2
    if(sign(f(c))==sign(f(a)))
      a <- c
    else
      b <- c
  }
  (a+b)/2 # bisection function
}

fsig <- function(sigma)
  bsc(30, 0.25, 0.0, 25, 0.06, sigma, 0.02) - 2.5

bisection(fsig, 0.0001, 0.1, 0.001)

dfsig <- function(sigma)
  bsc(30, 0.25, 0, 25, 0.06, sigma, 0.02)
fig <- function(sigma)
  bsc(30, 0.25, 0, 25, 0.06, sigma, 0.02) - 2.5

u<- 500
x<- 0.5 # initial guess
while(abs(u)/abs(x) > 1e-6){
  u <- fsig(x)/dfsig(x)
  x <- x-u
}
x