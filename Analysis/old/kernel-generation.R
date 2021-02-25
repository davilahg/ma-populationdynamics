nplot <- nlevels(as.factor(as.character(s1.db$plot)))
k.p.list <- list()
mat.p <- 5 # mature plot
s1.db <- subset(s1.db, !is.na(ln.h2))
s1.db$plot <- as.factor(as.character(s1.db$plot))
s1.db$Census <- as.factor(as.character(s1.db$Census))                  
for (p in 1:nplot) {
    Mim.p <- subset(s1.db, p == plot)
    min.a <- min(Mim.p$Age)
    max.a <- max(Mim.p$Age)
    if(p != 5) {
          # Growth
          GGG <- matrix(NA, ncol = max.a-min.a, nrow = M)
          for (a in (min.a+1):max.a)
            GGG[, a-min.a] <- predict(lmm.g, newdata = data.frame(ln.h1 = X.pred, Age = rep(a-1, M), plot = rep(p, M), Census = rep(Mim.p$Census[which(Mim.p$Age == a)][1], M)), type = "response", re.form = NULL)
          MUG <- matrix(NA, ncol = max.a-min.a, nrow = m)
          for (a in (min.a+1):max.a) {
            init.g <- 1
            for (i in 1:m) {
              MUG[i, a-min.a] <- mean(GGG[init.g:(init.g+sM-1), a-min.a])
              init.g <- init.g + sM
            }
          }
          sdg <- sigma(lmm.g)
          G <- array(NA, dim = c(max.a-min.a, m, m))
          for (a in (min.a+1):max.a) {
            for (i in 1:m) {
            G[(a-min.a), 1, i] <- pnorm(e.pred[2], MUG[i, (a-min.a)], sdg)
            for (j in 2:(m-1))
              G[(a-min.a), j, i] <- pnorm(e.pred[j+1], MUG[i, (a-min.a)], sdg) - pnorm(e.pred[j], MUG[i, (a-min.a)], sdg)
            G[(a-min.a), m, i] <- 1-pnorm(e.pred[m], MUG[i, (a-min.a)], sdg)
            }
          } 
          # Survival
          SSS <- matrix(NA, ncol = max.a-min.a, nrow = M)
          for (a in (min.a+1):max.a)
            SSS[ ,a-min.a] <- predict(glmm.s, type = "response", newdata = data.frame(ln.h1 = X.pred, Age = rep(a-1, M), plot = rep(p, M), Census = rep(Mim.p$Census[which(Mim.p$Age == a)][1], M)), re.form = NULL)
          MUS <- matrix(NA, ncol = max.a-min.a, nrow = m)
          for (a in (min.a+1):max.a) {
            init.s <- 1
            for (i in 1:m) {
              MUS[i, a-min.a] <- mean(SSS[init.s:(init.s+sM-1), a-min.a])
              init.s <- init.s + sM
            }
          }
          S <- MUS
          # Reproduction probability
          FF1 <- matrix(NA, ncol = max.a-min.a, nrow = M)
          for (a in (min.a+1):max.a) {
            if (all(p != as.integer(f1.db$plot))) {
              FF1[ ,a-min.a] <- predict(glmm.f1, type = "response", newdata = data.frame(ln.h1 = X.pred, Age = rep(a-1, M)), re.form = NA)
              } else {
                FF1[ ,a-min.a] <- predict(glmm.f1, type = "response", newdata = data.frame(ln.h1 = X.pred, Age = rep(a-1, M), plot = rep(as.factor(p), M)), re.form = NULL)
              }  
          }
          MUF1 <- matrix(NA, ncol = max.a-min.a, nrow = m)
          for (a in (min.a+1):max.a) {
            init.f1 <- 1
            for (i in 1:m) {
              MUF1[i, a-min.a] <- mean(FF1[init.f1:(init.f1+sM-1), a-min.a])
              init.f1 <- init.f1 + sM
            }
          }
          F1 <- MUF1
          # Fruitset by individual
          FF2 <- matrix(NA, ncol = max.a-min.a, nrow = M)
          for (a in (min.a+1):max.a) {
            if (all(p != as.integer(f2.db$plot))) {
              FF2[ ,a-min.a] <- predict(gam.f2$gam, type = "response", newdata = data.frame(ln.h1 = X.pred, Age = rep(a-1, M)), re.form = NA)
              } else {
                FF2[ ,a-min.a] <- predict(gam.f2$gam, type = "response", newdata = data.frame(ln.h1 = X.pred, Age = rep(a-1, M), plot = rep(as.factor(p), M)), re.form = NULL)
              }  
          }
          MUF2 <- matrix(NA, ncol = max.a-min.a, nrow = m)
          for (a in (min.a+1):max.a) {
            init.f2 <- 1
            for (i in 1:m) {
              MUF2[i, a-min.a] <- mean(FF2[init.f2:(init.f2+sM-1), a-min.a])
              init.f2 <- init.f2 + sM
            }
          }
          F2 <- MUF2
          # Seed number by fruit
          FF3 <- matrix(NA, ncol = max.a-min.a, nrow = M)
          for (a in (min.a+1):max.a) {
            if (all(p != f3.db$Plot)) {
              FF3[ ,a-min.a] <- predict(glmm.f3, type = "response", newdata = data.frame(ln.h1 = X.pred, Age = rep(a-1, M)), re.form = NA)
              } else {
                FF3[ ,a-min.a] <- predict(glmm.f3, type = "response", newdata = data.frame(ln.h1 = X.pred, Age = rep(a-1, M), Plot = rep(as.factor(p), M)), re.form = NULL)
              }  
          }
          MUF3 <- matrix(NA, ncol = max.a-min.a, nrow = m)
          for (a in (min.a+1):max.a) {
            init.f3 <- 1
            for (i in 1:m) {
              MUF3[i, a-min.a] <- mean(FF3[init.f3:(init.f3+sM-1), a-min.a])
              init.f3 <- init.f3 + sM
            }
          }
          F3 <- MUF3
          # establishment probability
          F4 <- predict(gam.f4, newdata = data.frame(x = f4.pred), type = "response")
          # recruit size
          F5 <- den.f5$y
          # kernel
          k.i.j.a <- s.i.j.a <- g.i.j.a <- p.i.j.a <- f.i.j.a <- F5.i.j.a <- array(NA, dim = c(max.a-min.a, m, m)) 
          s.i.a <- MUS
          f.i.a <- matrix(NA, nrow = m, ncol = max.a-min.a)
          for (a in (min.a+1):max.a) {
            for (i in 1:m) {
             f.i.a[i, a-min.a] <- F1[i, a-min.a]*F2[i, a-min.a]*F3[i, a-min.a]*F4[i]
             for (j in 1:m) {
              p.i.j.a[a-min.a, j, i] <- s.i.a[i, a-min.a]*G[a-min.a, j, i]
              f.i.j.a[a-min.a, j, i] <- f.i.a[i, a-min.a]*F5[j]
              k.i.j.a[a-min.a, j, i] <- p.i.j.a[a-min.a, j, i]+f.i.j.a[a-min.a, j, i]
             }
            }
           }
           k.p.list[[p]] <- k.i.j.a
          }
         }
save(k.p.list, file = "kernel-plot-list.RData")
         
      
      
      
      
