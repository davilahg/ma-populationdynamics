library(lme4)
library(AICcmodavg)

data <- read.csv("complete_data.csv")
lhmean <- mean(c(data$ln.h1, data$ln.h2), na.rm = TRUE)
lhsd <- sd(c(data$ln.h1, data$ln.h2), na.rm = TRUE)
data_scaled <- transform(data, ln.h1_scl = (ln.h1-lhmean)/lhsd, ln.h2_scl = (ln.h2-lhmean)/lhsd, Age_scl = scale(Age))
survival <- subset(data_scaled, !is.na(h1))

fix <- "sup~ln.h1_scl*Age_scl"
rand <- c("1", "Age_scl", "ln.h1_scl", "Age_scl:ln.h1_scl")
v.c <- c(1, 1, 1, 1)
v.p <- c(1, 1, 1, 1)
v <- c(v.c, v.p)
mods <- list()
cont <- 1
formula <- paste0(fix, "+")
if (any(v.c != 0)) {
	formula <- paste0(formula, "(")
 	for (c in 1:length(v.c)) {
		formula <- paste0(formula, rand[c])
	 	if (length(which(v.c != 0)) > 1 & c < length(v.c))
			formula <- paste0(formula, "+")
	}
	formula <- paste0(formula, "|", "Census)")
}
formula <- paste0(formula, "+")
formula <- paste0(formula, "(")
for (p in 1:length(v.p)) {
	formula <- paste0(formula, rand[p])
 	if (length(which(v.p != 0)) > 1 & p < length(v.p))
		formula <- paste0(formula, "+")
}
formula <- paste0(formula, "|", "PlotCode)")
mod <- glmer(as.formula(formula), survival, binomial)		 
if (length(mod@optinfo$conv$lme4) == 0) {
	mods[[cont]] <- mod
	cont <- cont + 1
}
par <- c(1, 2, 3, 4)
for (i in par) {
	choo <- choose(length(v.c) + length(v.p), i) # Calculates the number of sets with n elements than can be chosen from a set with k elements
	combs <- combn(length(v.c) + length(v.p), i) # Generate all combinations of the elements of x taken m at a time. If x is a positive integer, returns all combinations of the elements of seq(x) taken m at a time.
	for (j in 1:choo) {
		v.comb <- v
  		v.comb[combs[, j]] <- 0
		v.c.comb <- v.comb[1:length(v.c)]
		v.p.comb <- v.comb[(1 + length(v.c)):length(v)]
		formula <- paste0(fix, "+")
		if (any(v.c.comb != 0)) {
			formula <- paste0(formula, "(")
		 	for (c in 1:length(v.c.comb)) {
			  	if (v.c.comb[c] == 1) {
					formula <- paste0(formula, rand[c])
			  	} else {
			  		if (rand[c] == "1")
			  			formula <- paste0(formula, "0+")
			  	}
				if (length(which(v.c.comb != 0)) > 1 & any(na.omit(v.c.comb[c+1:length(v.c.comb)])[1] == 1) & c < length(v.c.comb))
					formula <- paste0(formula, "+")
			}
			formula <- paste0(formula, "|", "Census)")
		}
		if (any(v.c.comb != 0) & any(v.p.comb != 0))
			formula <- paste0(formula, "+")
		if (any(v.p.comb != 0)) {
			formula <- paste0(formula, "(")
		 	for (p in 1:length(v.p.comb)) {
			  	if (v.p.comb[p] == 1) {
					formula <- paste0(formula, rand[p])
			  	} else {
			  		if (rand[p] == "1")
			  			formula <- paste0(formula, "0+")
			  	}
			 	if (length(which(v.p.comb != 0)) > 1 & any(na.omit(v.p.comb[p+1:length(v.p.comb)])[1] == 1) & p < length(v.p.comb))
					formula <- paste0(formula, "+")
			}
			formula <- paste0(formula, "|", "PlotCode)")
		}
		cat(paste0("i = ", i, ", j = ", j, "\n"))
		mod <- glmer(as.formula(formula), survival, binomial)
		cat("\n")
		if (length(mod@optinfo$conv$lme4) == 0) { ##
			cat(paste0("CONV: ", formula, "\n"))
			mods[[cont]] <- mod
			cont <- cont + 1
		}
	}
	if (cont > 1)
		break
}


CONV: sup~ln.h1_scl*Age_scl+(1+ln.h1_scl+Age_scl:ln.h1_scl|Census)+(1+Age_scl+Age_scl:ln.h1_scl|PlotCode)
CONV: sup~ln.h1_scl*Age_scl+(1+Age_scl+ln.h1_scl|Census)+(1+ln.h1_scl+Age_scl:ln.h1_scl|PlotCode)
CONV: sup~ln.h1_scl*Age_scl+(1+Age_scl+ln.h1_scl|Census)+(1+Age_scl+ln.h1_scl|PlotCode) # !!!

> aictab(mods)

Model selection based on AICc:

      K    AICc Delta_AICc AICcWt Cum.Wt       LL
Mod3 16 3277.63       0.00      1      1 -1622.78
Mod1 16 3309.61      31.98      0      1 -1638.77
Mod2 16 3312.57      34.93      0      1 -1640.25

Warning message:
In aictab.AICglmerMod(mods) : 
Model names have been supplied automatically in the table


glmm.s <- glmer(sup~ln.h1_scl*Age_scl+(1+Age_scl+ln.h1_scl|Census)+(1+Age_scl+ln.h1_scl|PlotCode), survival, binomial)

