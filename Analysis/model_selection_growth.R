library(lme4)
library(AICcmodavg)

data <- read.csv("Data/cleaned_data/complete_data.csv")
lhmean <- mean(c(data$ln.h1, data$ln.h2), na.rm = TRUE)
lhsd <- sd(c(data$ln.h1, data$ln.h2), na.rm = TRUE)
data_scaled <- transform(data, ln.h1_scl = (ln.h1-lhmean)/lhsd, ln.h2_scl = (ln.h2-lhmean)/lhsd, Age_scl = scale(Age))
growth <- subset(data_scaled, !is.na(h1) & !is.na(h2) & sup == 1)

fix <- "ln.h2_scl~ln.h1_scl*Age_scl"
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
mod <- lmer(as.formula(formula), growth)		 
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
		mod <- lmer(as.formula(formula), growth)
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

CONV: ln.h2_scl~ln.h1_scl*Age_scl+(0++Age_scl+ln.h1_scl+Age_scl:ln.h1_scl|Census)+(1+Age_scl+Age_scl:ln.h1_scl|PlotCode) 	13*3+16*3 = 87
CONV: ln.h2_scl~ln.h1_scl*Age_scl+(1+ln.h1_scl+Age_scl:ln.h1_scl|Census)+(1+Age_scl+Age_scl:ln.h1_scl|PlotCode) 			13*3+16*3 = 87
CONV: ln.h2_scl~ln.h1_scl*Age_scl+(1+Age_scl+Age_scl:ln.h1_scl|Census)+(0++Age_scl+ln.h1_scl+Age_scl:ln.h1_scl|PlotCode)	13*3+16*3 = 87
CONV: ln.h2_scl~ln.h1_scl*Age_scl+(1+Age_scl+Age_scl:ln.h1_scl|Census)+(1+ln.h1_scl+Age_scl:ln.h1_scl|PlotCode)				13*3+16*3 = 87
CONV: ln.h2_scl~ln.h1_scl*Age_scl+(1+Age_scl+Age_scl:ln.h1_scl|Census)+(1+Age_scl+Age_scl:ln.h1_scl|PlotCode)				13*3+16*3 = 87
CONV: ln.h2_scl~ln.h1_scl*Age_scl+(1+Age_scl+ln.h1_scl+Age_scl:ln.h1_scl|Census)+(0++Age_scl+ln.h1_scl|PlotCode)			13*4+16*2 = 84

> aictab(aaa)

Model selection based on AICc:

      K    AICc Delta_AICc AICcWt Cum.Wt   Res.LL
Mod2 17 5278.00       0.00      1      1 -2621.95
Mod4 17 5360.84      82.84      0      1 -2663.37
Mod5 17 5377.24      99.24      0      1 -2671.58
Mod3 17 5383.73     105.73      0      1 -2674.82
Mod1 17 5483.25     205.25      0      1 -2724.58

Warning messages:
1: In aictab.AIClmerMod(aaa) : 
Model names have been supplied automatically in the table

2: In aictab.AIClmerMod(aaa) : 
Model selection for fixed effects is only appropriate with ML estimation:
REML (default) should only be used to select random effects for a constant set of fixed effects



lmm.g <- lmer(ln.h2_scl~ln.h1_scl*Age_scl+(1+ln.h1_scl+Age_scl:ln.h1_scl|Census)+(1+Age_scl+Age_scl:ln.h1_scl|PlotCode), growth)

