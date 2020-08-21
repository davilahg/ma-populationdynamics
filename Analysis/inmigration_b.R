# migration estimating function (for population size) ## + c (1st year) ## type = 0
{
  MSE_0 <- function(parm) {
    dist.N <- c()
    n.0 <- which(s1.db$Age == 0)
    n.0.h <- log(s1.db[n.0,]$h2)
    n.0.v <- hist(n.0.h, breaks = e.pred, plot = FALSE)$counts
    n.0.v <- n.0.v*area_scale
    c <- parm
    init.n.a.v <- n.0.v+c*F5
    pred.N.total <- c()
    for (a in Age.pred) {
      n.a.v <- k.i.j.a[a,,]%*%init.n.a.v
      pred.N.total <- c(pred.N.total, sum(n.a.v))
      init.n.a.v <- n.a.v
    }
    tot_sqr <- 0
    for (a in Age.pred) {
      sqr <- (pred.N.total[a] - total_n[a])**2
      tot_sqr <- sum(tot_sqr, sqr)
      MSE <- sqrt(tot_sqr)
    }
    cat(paste0("c = ", c, ", MSE = ", MSE, "\n"))
    return(MSE)
  }
}
# migration estimating function (for population size) ## + b0 ## type = 1
{
  obs.N.total <- data.frame(Age.pred,N.list)
  names(obs.N.total) <- c("Age", "N")
  distN1 <- function(parm) {
    dist.N <- c()
    n.0 <- which(s1.db$Age == 0)
    n.0.h <- log(s1.db[n.0,]$h2)
    n.0.v <- hist(n.0.h, breaks = e.pred, plot = FALSE)$counts
    b0 <- parm
    c <- b0 # :(
    init.n.a.v <- n.0.v+c*F5
    pred.N.total <- c()
    for (a in 1:Age.mature) {
      n.a.v <- k.i.j.a[a,,]%*%init.n.a.v
      n.a.v <- n.a.v+c*F5
      pred.N.total <- c(pred.N.total, sum(n.a.v))
      init.n.a.v <- n.a.v
    }
    for (p in 1:nplot) {
      dist.N.p <- 0
      N.p <- obs.N.plot[[p]]
      max.a <- max(N.p$Age)
      for (j in 1:nrow(N.p))
      dist.N.p <- dist.N.p + (N.p$N[j]-pred.N.total[max.a-nrow(N.p)+j])^2
      dist.N <- c(dist.N, dist.N.p)
    }
    dist.N <- sqrt(sum(dist.N))
    cat(paste0("b0 = ", b0, ", dist = ", dist.N, "\n"))
    return(dist.N)
  }
}
# migration estimating function (for population size) ## + b0 + b1*Age ## type = 2
{
  distN2 <- function(parm) {
    dist.N <- c()
    n.0 <- which(s1.db$Age == 0)
    n.0.h <- log(s1.db[n.0,]$h2)
    n.0.v <- hist(n.0.h, breaks = e.pred, plot = FALSE)$counts
    b0 <- parm[1]
    b1 <- parm[2]
    c <- b0
    init.n.a.v <- n.0.v+c*F5
    pred.N.total <- c()
    for (a in 1:Age.mature) {
      n.a.v <- k.i.j.a[a,,]%*%init.n.a.v
      c <- b0 + b1*Age.pred[a]
      n.a.v <- n.a.v+c*F5
      pred.N.total <- c(pred.N.total, sum(n.a.v))
      init.n.a.v <- n.a.v
    }
    for (p in 1:nplot) {
      dist.N.p <- 0
      N.p <- obs.N.plot[[p]]
      max.a <- max(N.p$Age)
      for (j in 1:nrow(N.p))
      dist.N.p <- dist.N.p + (N.p$N[j]-pred.N.total[max.a-nrow(N.p)+j])^2
      dist.N <- c(dist.N, dist.N.p)
    }
    dist.N <- sqrt(sum(dist.N))
    cat(paste0("b0 = ", b0, ", b1 = ", b1, ", dist = ", dist.N, "\n"))
    return(dist.N)
  }
}
# migration estimating function (for population size) ## + b0 + b1*Age + b2*Age**2 ## type = 3
{
  distN3 <- function(parm) {
    dist.N <- c()
    n.0 <- which(s1.db$Age == 0)	#
    n.0.h <- log(s1.db[n.0,]$h2)	#
    n.0.v <- hist(n.0.h, breaks = e.pred, plot = FALSE)$counts
    b0 <- parm[1]
    b1 <- parm[2]
    b2 <- parm[3]
    c <- b0
    init.n.a.v <- n.0.v+c*F5
    pred.N.total <- c() # size structure vector by year, WM = with migration
    for (a in 1:Age.mature) {
      n.a.v <- k.i.j.a[a,,]%*%init.n.a.v
      c <- b0 + b1 * Age.pred[a] + b2 * Age.pred[a]**2
      n.a.v <- n.a.v+c*F5
      pred.N.total <- c(pred.N.total, sum(n.a.v))
      init.n.a.v <- n.a.v			# add size structure
    }
    for (p in 1:nplot) {
      dist.N.p <- 0
      N.p <- obs.N.plot[[p]]
      max.a <- max(N.p$Age)
      for (j in 1:nrow(N.p))
      dist.N.p <- dist.N.p + (N.p$N[j]-pred.N.total[max.a-nrow(N.p)+j])^2
      dist.N <- c(dist.N, dist.N.p)
    }
    dist.N <- sqrt(sum(dist.N))
    cat(paste0("b0 = ", b0, " , b1 = ", b1, " , b2 = ", b2,", dist = ", dist.N, "\n"))
    return(dist.N)
  }
}
# migration estimating function (for population size) ##  exp (+ b0 + b1*Age - b2*Age**2) ## type = 4
{
  distN4 <- function(parm) {
    dist.N <- c()
    n.0 <- which(s1.db$Age == 0)	#
    n.0.h <- log(s1.db[n.0,]$h2)	#
    n.0.v <- hist(n.0.h, breaks = e.pred, plot = FALSE)$counts
    b0 <- parm[1]
    b1 <- parm[2]
    b2 <- parm[3]
    c <- exp(b0)
    init.n.a.v <- n.0.v+c*F5
    pred.N.total <- c() # size structure vector by year, WM = with migration
    for (a in 1:Age.mature) {
      n.a.v <- k.i.j.a[a,,]%*%init.n.a.v
      c <- exp(b0 + b1 * Age.pred[a] - b2 * Age.pred[a]**2)
      n.a.v <- n.a.v+c*F5
      pred.N.total <- c(pred.N.total, sum(n.a.v))
      init.n.a.v <- n.a.v			# add size structure
    }
    for (p in 1:nplot) {
      dist.N.p <- 0
      N.p <- obs.N.plot[[p]]
      max.a <- max(N.p$Age)
      for (j in 1:nrow(N.p))
      dist.N.p <- dist.N.p + (N.p$N[j]-pred.N.total[max.a-nrow(N.p)+j])^2
      dist.N <- c(dist.N, dist.N.p)
    }
    dist.N <- sqrt(sum(dist.N))
    cat(paste0("b0 = ", b0, " , b1 = ", b1, " , b2 = ", b2,", dist = ", dist.N, "\n"))
    return(dist.N)
  }
}
# migration estimating function (for population size) ##  exp (+ b0 - b1*Age) ## type = 5
{
  distN5 <- function(parm) {
    dist.N <- c()
    n.0 <- which(s1.db$Age == 0)	#
    n.0.h <- log(s1.db[n.0,]$h2)	#
    n.0.v <- hist(n.0.h, breaks = e.pred, plot = FALSE)$counts
    b0 <- parm[1]
    b1 <- parm[2]
    c <- exp(b0)
    init.n.a.v <- n.0.v+c*F5
    pred.N.total <- c() # size structure vector by year, WM = with migration
    for (a in 1:Age.mature) {
      n.a.v <- k.i.j.a[a,,]%*%init.n.a.v
      c <- exp(b0 - b1 * Age.pred[a])
      n.a.v <- n.a.v+c*F5
      pred.N.total <- c(pred.N.total, sum(n.a.v))
      init.n.a.v <- n.a.v			# add size structure
    }
    for (p in 1:nplot) {
      dist.N.p <- 0
      N.p <- obs.N.plot[[p]]
      max.a <- max(N.p$Age)
      for (j in 1:nrow(N.p))
      dist.N.p <- dist.N.p + (N.p$N[j]-pred.N.total[max.a-nrow(N.p)+j])^2
      dist.N <- c(dist.N, dist.N.p)
    }
    dist.N <- sqrt(sum(dist.N))
    cat(paste0("b0 = ", " , b1 = ", b1,", dist = ", dist.N, "\n"))
    return(dist.N)
  }
}
# N plot function
{
  plot_N <- function(par, type) {
    lam.list <- c()			# create lambda vector
    n.0 <- which(s1.db$Age == 0)	# get row numbers with age = 0
    n.0.h <- log(s1.db[n.0,]$h2)	#
    n.0.v <- hist(n.0.h, breaks = e.pred, plot = FALSE)$counts
    n.0.v <- n.0.v*area_scale
    if (type == 4 | type == 5) {
      c <- exp(par[1])
    } else {
      c <- par[1]
    }
    init.n.a.v <- n.0.v+c*F5
    size.v.a.WM <- list(init.n.a.v) # size structure vector by year, WM = with migration
    n.list <- c(sum(init.n.a.v))
    for (a in 1:Age.mature) {
      n.a.v <- k.i.j.a[a,,]%*%init.n.a.v
      if (type == 0) {
        c <- 0
      } else if (type == 1) {
        c <- par[1]
      } else if (type == 2) {
        c <- par[1] + Age.pred[a]*par[2]
      } else if (type == 3) {
        c <- par[1] + Age.pred[a]*par[2] + par[3]*Age.pred[a]**2
      } else if (type == 4) {
        c <- exp(par[1] + Age.pred[a]*par[2] + par[3]*Age.pred[a]**2)
      } else if (type == 5) {
        c <- exp(par[1] + Age.pred[a]*par[2])
      } else {
        assert("Wrong type")
      }
      n.a.v <- n.a.v+c*F5
      lam.a <- sum(n.a.v)/sum(init.n.a.v)
      size.v.a.WM[[a+1]] <- n.a.v
      n.list <- c(n.list, sum(n.a.v))
      lam.list <- c(lam.list, lam.a)
      init.n.a.v <- n.a.v
    }
    #plot
    n.list <- n.list[1:100]
    N_change <- plot(0:100, total_n, type = 'l', col = 'red', xlab = 'successional age', ylab = 'pop. growth rate')
    lines(Age.pred, n.list, col = 'blue')
    return(N_change)
  }
}
#
# + c (1st year) ## type = 0
opt.lbf0 <- optim(0, MSE_0, method = 'L-BFGS-B', lower = 0, upper = 1000)
c_0 <- opt.lbf0$par
#
# + b0 ## type = 1
opt.lbf1 <- optim(0, distN1, method = "L-BFGS-B", lower = 0, upper = 500)
c_1 <- opt.lbf1$par
#
# + b0 + b1*Age ## type = 2
opt.lbf2 <- optim(c(0, -2), distN2, method = "L-BFGS-B", lower = c(0, -15), upper = c(500, 5))
c_2 <- opt.lbf2$par
#
# + b0 + b1*Age - b2*Age**2 ## type = 3
opt.lbf3 <- optim(c(0, -2, 1), distN3, method = "L-BFGS-B", lower = c(0, -15, -15), upper = c(500, 5, 5))
c_3 <- opt.lbf3$par
#
#  exp (+ b0 + b1*Age - b2*Age**2) ## type = 4
opt.lbf4 <- optim(c(0, 1, 1), distN4, method = "L-BFGS-B", lower = c(0, -15, -5), upper = c(500, 5, 5))
c_4 <- opt.lbf4$par
#
# exp (+ b0 - b2*Age**2) ## type = 5
opt.lbf5 <- optim(c(0, 1, 1), distN5, method = "L-BFGS-B", lower = c(0, -15), upper = c(500, 5))
c_5 <- opt.lbf5$par
