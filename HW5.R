# 5(A) R function
bsc <- function(S, Tim, t, K, r, s, q){
  d1 <- (log(S/K) + (r - q + 0.5*s^2)*(Tim - t))/(s*sqrt(Tim - t))
  d2 <- d1 - s*sqrt(Tim - t)
  S*exp(-q*(Tim - t))*pnorm(d1) - K*exp(-r*(Tim - t))*pnorm(d2)
}

# 5(B)
sigmas <- seq(0.05, 0.5, by = 0.01)
fsig <- bsc(50, 0.5, 0.0, 45, 0.06, sigmas, 0.02) - 8
plot(sigmas, fsig, type = "l", col = "red", 
     xlab = expression(sigma), ylab = expression(f(sigma)), lwd = 2)

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
  bsc(50, 0.5, 0.0, 45, 0.06, sigma, 0.02) - 8

bisection(fsig, 0.3, 0.4, 0.001)

# 5(C) check good estiamte or not 
bsc(50, 0.5, 0.0, 45, 0.06, 0.3433594, 0.02) - 8


# 6(B)
g <- function(x){
  x^3 - 5*x^2 -7*x # g(x) in the problem 
}
y <- seq(-2, 1, by = 0.01) # make y a sequence from -2 to 1, step is 0.01
f <- g(y) - 1 # at each point of y, get the value - 1 named as f
plot(y, f, type = "l", col = "red") # plot the graph 
abline(a = 0, b = 0) # draw the line g = 0

g1 <- function(x)
  g(x) - 1
bisection(g1, -1.5, -0.5, 0.00001) # use the bsc nethod get one solution

bisection(g1, -0.4, 0.0, 0.00001) # use bsc get another solution

# choice for the initial interval
# First, I draw the graph of y and g (g is the response variable)
# Second, I draw the line f = 0, and see in which interval f = 0
# The interval that f = 0 is what we need to find

# 7(b)
f <- function(x) # original function
  -x^5 - cos(x)
df <- function(x) # first derivative
  -5*x^4 + sin(x)

x = 1 # intial value of x
u = f(1)/df(1) # intial f(x)/f'(x)
while (abs(u)/abs(x) > 1e-4) { # tolerance is 1e-4 where we want to stop
  u <- f(x)/df(x) # f(x)/f'(x)
  x <- x - u # x_k+1 = x_k - f(x_k)/f'(x_k)
}
# show final answer
x

# 7(c)
# we can make the tolerance be smaller to get a more accurate answer
# after I run the code, I found the answer gets more accurate after the 
# sixth digit and more, even though the iteration number may increase

x = 1 # intial value of x
u = f(1)/df(1) # intial f(x)/f'(x)
while (abs(u)/abs(x) > 1e-6) { # tolerance is 1e-4 where we want to stop
  u <- f(x)/df(x) # f(x)/f'(x)
  x <- x - u # x_k+1 = x_k - f(x_k)/f'(x_k)
}
# show final answer
x
