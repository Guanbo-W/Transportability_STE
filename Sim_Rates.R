## one simulation replication
## target parameter is phi(a = 1, s = 1, ~x = x0 = 1)
replication <- function(n, err) {
  
  h <- 2.5
  X0 <-  rbinom(n, 1, 0.5) # covariate define \tilde{x}
  X1 <- runif(n) # another confounder
  S <- rbinom(n, 1, plogis(-0.2 + 1.2 * X1 + 0.5 * X0)) # 2 data sources (1 and 0)
  epsilon <- 0.10 # bounding treatment probabilities away from zero
  A <- rbinom(n, 1, pmin(pmax(epsilon, ifelse(S == 1,
                                              plogis(0.8 - 0.8 * X1 + 0.9 * X0),
                                              plogis(-0.8 + 0.8 * X1 - 0.9 * X0))),
                         1 - epsilon))
  Y <- rnorm(n, 5.2 + (1.2 - 0.6 * X0)*A + X0 - 1.2 * X1, 1)
  
  q.1 <- plogis(-0.2 + 1.2 * X1 + 0.5 * X0 + 1.3 * h * rnorm(n, 1/n^err, 1/n^err))
  eta.1 <- plogis(qlogis(plogis(-0.2 + 1.2 * X1 + 0.5 * X0) * plogis(0.8 - 0.8 * X1 + 0.9 * X0) +
    (1 - plogis(-0.2 + 1.2 * X1 + 0.5 * X0)) * plogis(-0.8 + 0.8 * X1 - 0.9 * X0)) +
                    1.3 * h * rnorm(n, 1/n^err, 1/n^err))
  mu.1 <- 5.2 + (1.2 - 0.6 * X0) + X0 - 1.2 * X1 + h * rnorm(n, 1/n^err, 1/n^err)
  
  ## plug-in estimator
  plug.l <- 1 * (S == 1 & X0 == 1) * mu.1 / mean(S == 1 & X0 == 1)
  
  ## IF-based estimator
  IF.l <- (1 * (S == 1 & X0 == 1) * mu.1 + 
             1 * (A == 1 & X0 == 1) * q.1 * (Y - mu.1) / eta.1) / 
    mean(S == 1 & X0 == 1)
  
  return(c(plug = mean(plug.l), IF = mean(IF.l)))
  
}

set.seed(590)
R <- 5000
n <- 1000
err <- seq(0.10, 0.5, by = 0.05)

res <- sapply(1:length(err), function(l) {
  err.l <- err[l]
  sapply(1:R, FUN = function(j) {
    replication(n, err.l)
  })
}, simplify = "array")
res <- aperm(res, c(3, 2, 1))

### FIND TRUTH
N <- 10000000
X0 <-  rbinom(N, 1, 0.5) # covariate define \tilde{x}
X1 <- runif(N) # another confounder
S <- rbinom(N, 1, plogis(-0.2 + 1.2 * X1 + 0.5 * X0)) # 2 data sources (1 and 0)
epsilon <- 0.10 # bounding treatment probabilities away from zero
A <- rbinom(N, 1, pmin(pmax(epsilon, ifelse(S == 1,
                                            plogis(0.8 - 0.8 * X1 + 0.9 * X0),
                                            plogis(-0.8 + 0.8 * X1 - 0.9 * X0))),
                       1 - epsilon))
Y <- rnorm(N, 5.2 + (1.2 - 0.6 * X0)*A + X0 - 1.2 * X1, 1)
truth <- mean(5.2 + (1.2 - 0.6) + 1 - 1.2 * X1[S == 1])

rmse <- function(res) {
  sqrt(colMeans((res - truth)^2))
}

res.rmse <- cbind.data.frame(err, t(apply(res, 1, rmse)))

par(mar = c(4,4,1,1), mfrow = c(1,1))
plot(res.rmse$err, res.rmse$plug, type = "b", xlim = c(0.10, 0.5), pch=20,
     ylim = c(0, max(c(res.rmse$IF, res.rmse$plug))), 
     lty = "dashed", col = "red", ylab = "RMSE", xlab = "Error Rate")
lines(res.rmse$err, res.rmse$IF, type = "b", pch = 20)
legend("topright", lty = c("dashed","solid"), pch = 20,
       col = c("red", "black"), bty = "n",
       legend = c("Plug-in", "IF-based"))
