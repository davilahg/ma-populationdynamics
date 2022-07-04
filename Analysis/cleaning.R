
## data cleaning
#
{
    Age.mature <- 100 # mature forest age
    M <- 20000 # big-matrix dimension
    raw_canopy_data <- read.csv("./Data/raw_data/raw_canopy.csv")
    raw_canopy_data$UID <- as.factor(as.character(raw_canopy_data$UID))
    raw_canopy_data <- raw_canopy_data[which(raw_canopy_data$Species == "Mimosa acantholoba"), ]
    isa.rows <- which(raw_canopy_data$Plot != "ISA") # removing data from ISA and ISX, age is unknown
    isx.rows <- which(raw_canopy_data$Plot != "ISX")
    raw_canopy_data <- raw_canopy_data[isa.rows,]
    raw_canopy_data <- raw_canopy_data[isx.rows,]
    raw_canopy_data <- droplevels(raw_canopy_data)
    raw_canopy_data <- raw_canopy_data[!is.na(raw_canopy_data$UID),]
    raw_canopy_data <- transform(raw_canopy_data, Age_init = 0, Fst_cens = 2003)
    for (i in 1:nrow(raw_canopy_data)) {
        if (raw_canopy_data$Plot[i] == "DIV") {raw_canopy_data$Fst_cens[i] <- 2005}
    }
    for (j in 1:nrow(raw_canopy_data)) { # setting initial ages
        if (raw_canopy_data$Plot[j] == "RIC") {raw_canopy_data$Age_init[j] <- 30} #1 Ricardo
        else if (raw_canopy_data$Plot[j] == "ISP") {raw_canopy_data$Age_init[j] <- 23} #2 Isaìas Plano
        else if (raw_canopy_data$Plot[j] == "SEP") {raw_canopy_data$Age_init[j] <- 40} #3 Severino Plano = Severino I
        else if (raw_canopy_data$Plot[j] == "ABE") {raw_canopy_data$Age_init[j] <- 10} #4 Abel
        else if (raw_canopy_data$Plot[j] == "SEL") {raw_canopy_data$Age_init[j] <- Age.mature} #5 Severino Loma = Severino II
        else if (raw_canopy_data$Plot[j] == "ISC") {raw_canopy_data$Age_init[j] <- 18} #6 Isaìas Carretera
        else if (raw_canopy_data$Plot[j] == "BES") {raw_canopy_data$Age_init[j] <- 11} #7 Bernardo
        else if (raw_canopy_data$Plot[j] == "TOA") {raw_canopy_data$Age_init[j] <- 0} #8 Toño Arriba = Toño de Leòn 0
        else if (raw_canopy_data$Plot[j] == "MAL") {raw_canopy_data$Age_init[j] <- 36} #9 Malaquìas
        else if (raw_canopy_data$Plot[j] == "HIL") {raw_canopy_data$Age_init[j] <- 5} #10 Hilario
        else if (raw_canopy_data$Plot[j] == "MAR") {raw_canopy_data$Age_init[j] <- 3} #11 Marciano
        else if (raw_canopy_data$Plot[j] == "FID") {raw_canopy_data$Age_init[j] <- 9} #12 Wilfrido
        else if (raw_canopy_data$Plot[j] == "ESS") {raw_canopy_data$Age_init[j] <- 16} #13 Escolar = Escolar Secundario
        else if (raw_canopy_data$Plot[j] == "TOB") {raw_canopy_data$Age_init[j] <- 1} #14 Toño Abajo = Toño de Leòn I
        else if (raw_canopy_data$Plot[j] == "DIA") {raw_canopy_data$Age_init[j] <- 7} #15 Dino agua caliente
        else if (raw_canopy_data$Plot[j] == "DIV") {raw_canopy_data$Age_init[j] <- 56} #17 Dino Vàlvula
    }
 raw_canopy_data <- transform(raw_canopy_data, Age = 0) # setting ages per census year
 for (k in 1:nrow(raw_canopy_data)) {
	 raw_canopy_data$Age[k] <- raw_canopy_data$Census[k] - raw_canopy_data$Fst_cens[k] + raw_canopy_data$Age_init[k]
	 if (raw_canopy_data$Age[k] > Age.mature) {
	 	 raw_canopy_data$Age[k] <- Age.mature
	 }
 }
 #
 raw_canopy_data <- transform(raw_canopy_data, sup = rep(NA, nrow(raw_canopy_data)), h1 = rep(NA, nrow(raw_canopy_data)), h2 = rep(NA, nrow(raw_canopy_data)))
 raw_canopy_data$Sup[is.na(raw_canopy_data$Sup)] <- 0
 new_canopy_data <- as.data.frame(matrix(NA, nrow = 0, ncol = ncol(raw_canopy_data)))
 names(new_canopy_data) <- names(raw_canopy_data)
 for (i in levels(raw_canopy_data$UID)) {
 	new_canopy_data.i <- subset(raw_canopy_data, UID == i)
 	if (nrow(new_canopy_data.i) > 1) {
 		new_canopy_data.i <- new_canopy_data.i[order(new_canopy_data.i$Census), ]
 		flag <- 1             
 		for (j in 1:nrow(new_canopy_data.i)) {
 			if (new_canopy_data.i$Sup[j] == 0) { # Remove data following a registered 0 in survival. Sometimes there are trees wich are still being registered after death, resulting in more 0's than real and overestimating death probability
                #if (flag == 1)
                    #flag <- 0
 				new_canopy_data.i <- new_canopy_data.i[1:j, ]
                break
 			}
 			#if (new_canopy_data.i$Sup[j] == 1)
                #if (flag == 0) {
                    #cat(paste0(new_canopy_data.i$UID[1], "\n"))   ##
                    #break
                #}
 		}
    }
    if (nrow(new_canopy_data.i) > 1) {
  		new_canopy_data.i$sup[1:nrow(new_canopy_data.i)] <- new_canopy_data.i$Sup[1:nrow(new_canopy_data.i)]
  		new_canopy_data.i$h1[2:nrow(new_canopy_data.i)] <- new_canopy_data.i$Height[1:(nrow(new_canopy_data.i))-1]
        new_canopy_data.i$h2[1:nrow(new_canopy_data.i)] <- new_canopy_data.i$Height[1:nrow(new_canopy_data.i)]
        new_canopy_data <- rbind(new_canopy_data, new_canopy_data.i)
   }
  }
 new_canopy_data <- transform(new_canopy_data, ln.h1 = log(h1), ln.h2 = log(h2), Age = Age)
 used_vars <- names(new_canopy_data[c("Age", "ln.h1", "ln.h2", "h1", "h2", "PlotCode", "sup", "UID", "Census", "maxDBH")])
 canopy <- new_canopy_data[ ,used_vars]
 #write.table(canopy, file = "./Data/cleaned_data/canopy.csv", sep = ",", row.names = FALSE)
}

# additional mature forest data
{
    mature_data <- read.csv('./Data/raw_data/mature_data.csv')
    mature_data$UID <- as.factor(as.character(mature_data$UID))
    mature_data$Height <- as.numeric(mature_data$Height)
    mature_data <- transform(mature_data, Age = 100, sup = rep(NA, nrow(mature_data)), h1 = rep(NA, nrow(mature_data)), h2 = rep(NA, nrow(mature_data)))
    mature_data$Sup[is.na(mature_data$Sup)] <- 0
    new_mature_data <- as.data.frame(matrix(NA, nrow = 0, ncol = ncol(mature_data)))
    names(new_mature_data) <- names(mature_data)
    for (i in levels(mature_data$UID)) {
        new_mature_data.i <- subset(mature_data, UID == i)
        if (nrow(new_mature_data.i) > 1) {
            new_mature_data.i <- new_mature_data.i[order(new_mature_data.i$Census), ]
            for (j in 1:nrow(new_mature_data.i)) {
                if (new_mature_data.i$Sup[j] == 0) {
                    new_mature_data.i <- new_mature_data.i[1:j, ]
                    break
                }
            }
            new_mature_data.i$sup[1:nrow(new_mature_data.i)] <- new_mature_data.i$Sup[1:nrow(new_mature_data.i)]
            new_mature_data.i$h1[2:nrow(new_mature_data.i)] <- new_mature_data.i$Height[1:(nrow(new_mature_data.i))-1]
            new_mature_data.i$h2[1:nrow(new_mature_data.i)] <- new_mature_data.i$Height[1:nrow(new_mature_data.i)]
            new_mature_data <- rbind(new_mature_data, new_mature_data.i)
        }
    }
    new_mature_data <- transform(new_mature_data, ln.h1 = log(h1), ln.h2 = log(h2), Age = Age)
    used_vars <- names(new_mature_data[c("Age", "ln.h1", "ln.h2", "h1", "h2", "PlotCode", "sup", "UID", "Census", "maxDBH")])
    mature <- new_mature_data[ , used_vars]
    #write.table(mature, file = './Data/cleaned_data/mature.csv', sep = ",", row.names = FALSE)
}

# understory
{
 understory <- read.csv("./Data/raw_data/raw_understory.csv")
 understory$Plot <- as.factor(as.character(understory$Plot))
 understory <- transform(understory, h1 = rep(NA, nrow(understory)), sup = rep(NA, nrow(understory)), h2 = rep(NA, nrow(understory)), Age = Ecenso, PlotCode = Plot, Census = Año)
 #understory <- understory[ ,used_vars]
 understory <- subset(understory, Supv != "Exc" & Supv != "NE")
 new_understory_data <- as.data.frame(matrix(NA, nrow = 0, ncol = ncol(understory)))
 names(new_understory_data) <- names(understory)
 for (i in 1:nrow(understory))
    if (understory$Supv[i] == "Dosel") {understory$Supv[i] <- 1}
 understory$UID <- as.factor(as.character(understory$UID))
 understory <- droplevels(understory)
 for (i in levels(understory$UID)) {
    understory.i <- subset(understory, UID == i)
    if (nrow(understory.i) > 1) {
        for (j in 1:nrow(understory.i)) {
            if (understory.i$Supv[j] == 0) {
                understory.i <- understory.i[1:j, ]
                break
            }
        }
        if (nrow(understory.i) > 1) {
            understory.i$sup[1:nrow(understory.i)] <- understory.i$Supv[1:nrow(understory.i)]
            understory.i$h1[2:nrow(understory.i)] <- understory.i$Altura[1:(nrow(understory.i))-1]
            understory.i$h2[1:nrow(understory.i)] <- understory.i$Altura[1:nrow(understory.i)]
            new_understory_data <- rbind(new_understory_data, understory.i)
        }
    }
 }
 new_understory_data <- transform(new_understory_data, ln.h1 = log(h1), ln.h2 = log(h2), maxDBH = DB.max)
 understory <- new_understory_data[ ,used_vars]
 #write.table(understory, file = "./Data/cleaned_data/understory.csv", sep = ",", row.names = FALSE)
 complete_data <- rbind(canopy, understory, mature)
 #write.table(complete_data, file = "./Data/cleaned_data/complete_data.csv", sep = ",", row.names = FALSE)
}









 
 und.00 <- as.data.frame(matrix(nrow = 0, ncol = ncol(understory)))
 und.00 <- transform(und.00, sup = rep(NA, nrow(und.00)), h1 = rep(NA, nrow(und.00)), h2 = rep(NA, nrow(und.00)))
 und.s <- as.data.frame(matrix(NA, nrow = 0, ncol = ncol(und.00)))
 
 for (i in levels(und.00$id)) {
 	und.id.s <- subset(und.00, id == i)
 	if (nrow(und.id.s) > 1) {
 		for (j in 1:nrow(und.id.s)) {
 			if (und.id.s$Supv[j] == 0) {
 				und.id.s <- und.id.s[1:j, ]
 	 break
 			}
 		}
 	if (nrow(und.id.s) > 1) {
 		und.id.s$sup[1:nrow(und.id.s)] <- und.id.s$Supv[1:nrow(und.id.s)] # !!! AQUÍ SE LE ESTÁ SUMANDO UN UNO POR ALGUNA RAZÓN, VOY A RESTÁRSELO ABAJO Y LO GUARDARÉ, PERO HAY QUE ARREGLARLO
 		und.id.s$h1[2:nrow(und.id.s)] <- und.id.s$Altura[1:nrow(und.id.s)-1]
 		und.id.s$h2[1:nrow(und.id.s)] <- und.id.s$Altura[1:nrow(und.id.s)]
 		und.s <- rbind(und.s, und.id.s)
 	  }
   }
  }
  for (i in 1:nrow(und.s)) # AQUÍ LE QUITO EL UNO QUE SE LE SUMA ABAJO
    if (und.s$sup[i] == 1) {und.s$sup[i] <- 0} else if (und.s$sup[i] == 2) {
      und.s$sup[i] <- 1
    }
  und.s <- transform(und.s, ln.h1 = log(h1), ln.h2 = log(h2), Age = Ecenso, Census = Año) #All variables
  und.S <- und.s[,c("ln.h1", "ln.h2", "sup", "Age", "PlotCode", "id", "Census", "h1", "h2")] #Only used variables
  #write.table(und.s, file = "und.s.csv", sep = ",", row.names = FALSE)
  #write.table(und.S, file = "und.S.csv", sep = ",", row.names = FALSE)
}
und.S <- read.csv("/home/gerardo/MEGA/Tesis/Datos/und.S.csv")
und.S <- read.csv("/Users/gdh/Google\ Drive/L\ 2017\ Gerardo\ Dávila/Datos/und1.S.csv")
# agregar a los individuos de edad = 0
oldund <- read.csv("Sotobosque.csv")
oldund <- subset(oldund, Age == 0)
oldund <- transform(oldund, id = as.factor(paste0("p", Plot_code, "i", Tag)), h2 = h1)
oldund <- transform(oldund, h1 = rep(NA, nrow(oldund)), ln.h1 = rep(NA, nrow(oldund)), sup = rep(1, nrow(oldund)), ln.h2 = log(h2), Plot_code = rep(8, nrow(oldund)))
oldund <- oldund[,c("ln.h1", "ln.h2", "sup", "Age", "Plot_code", "id", "Census", "h1", "h2")]
und.S <- rbind(und.S, oldund)

# FUNCIÓN PARA CALCULAR h ANTES DEL PRIMER REGISTRO
 #Se necesita cargar "Growth" para correr lmm.g y utilizarlo para calcular glmm.s, usar script "Mim0.R"
{
 ## ln.h1 = altura antes de ser registrado por primera vez, ln.h2 = altura en su primer registro
 geth1 <- function(ln.h2, Age, Plot, Census) { ## Census? Plot_code? id?
  b0 <- fx.coef[1] + pt.coef$int[which(pt.coef$Plot_code == Plot)] + ifelse(Census != 0, cs.coef$int[which(cs.coef$Census == Census)], 0)
  b1 <- fx.coef[2] + pt.coef$h1[which(pt.coef$Plot_code == Plot)] + ifelse(Census != 0, cs.coef$h1[which(cs.coef$Census == Census)], 0)
  b2 <- fx.coef[3]
  b3 <- fx.coef[4]
   ln.h1 <- (ln.h2-b0-b2*Age)/(b1+b3*Age) #calcular el modelo general
   return(ln.h1)
   }
}

# CALCULAR TAMAÑO ANTERIOR PARA PRIMER REGISTRO
{
 fx.coef <- fixef(lmm.g)
 pt.coef <- ranef(lmm.g)$Plot_code
 pt.coef <- transform(pt.coef, Plot_code = as.factor(row.names(pt.coef)))
 colnames(pt.coef) <- c("int", "h1", "Plot_code")
 cs.coef <- ranef(lmm.g)$Census
 cs.coef <- transform(cs.coef, Census = as.factor(row.names(cs.coef)))
 colnames(cs.coef) <- c("int", "h1", "Census")
 censuses <- as.factor(row.names(cs.coef))
 newSurvival <- Survival
 for (i in levels(Survival$id)) {
   S.i <- subset(Survival, id == i)
   ln.h2 <- S.i$ln.h2[1]
   Age <- S.i$Age[1]
   Plot <- S.i$Plot_code[1]
   Census <- ifelse(any(censuses == S.i$Census[1]-1) == TRUE, S.i$Census[1]-1, 0)
   if(Census != 0) {
   newSurvival[which(newSurvival$id == i)[1],]$ln.h1[1] <- geth1(ln.h2, Age, Plot, Census)
   newSurvival[which(newSurvival$id == i)[1],]$h1[1] <- exp(newSurvival[which(newSurvival$id == i)[1],]$ln.h1[1])
   }
   if(Census == 0) {
   ln.h1.new <- geth1(ln.h2, Age, Plot, Census)
   new.row <- as.data.frame(c(S.i$Age[1]-1, ln.h1.new, ln.h2, exp(ln.h1.new), exp(ln.h2), Plot, 1, 0, S.i$Census[1]-1))
   rownames(new.row) <- colnames(newSurvival)
   colnames(new.row) <- "1"
   new.row[which(rownames(new.row) == "id"), 1] <- i
   row <- as.data.frame(matrix(rep(NA, length(rownames(new.row))), nrow = 1, ncol = length(rownames(new.row))))
   for (j in 1:length(row)) {
       if (rownames(new.row)[j] != "id") {
       row[1, j] <- as.numeric(new.row[j, 1])
       } else {
       row[1, j] <- new.row[j, 1]
       }
    }
    colnames(row) <- rownames(new.row)
    newSurvival <- rbind(newSurvival, row)
   }
 }
 #write.table(newSurvival, file = "nSur.csv", sep = ",", row.names = FALSE)
}


# CONSTRUCCIÓN DE MATRICES
{
 for (p in i:nPlot_code) {

 }
 # Crecimiento
 GGG <- matrix(NA, ncol = Age.mature, nrow = M)
 for (a in 1:Age.mature)### parece que así es el bueno, con paréntesis
  GGG[, a] <- predict(lmm.g, newdata = data.frame(ln.h1 = X.pred, Age = rep(a-1, M)), type = "response", re.form = NA)
 MUG <- matrix(NA, ncol = Age.mature, nrow = m)
 for (a in 1:Age.mature) {#### tambièn aquì
  init.g <- 1
  for (i in 1:m) {
     MUG[i, a] <- mean(GGG[init.g:(init.g+sM-1), a])
   init.g <- init.g + sM
  }
 }
 sdg <- sigma(lmm.g)
 G <- array(NA, dim = c(Age.mature, m, m))
 for (a in 1:Age.mature) {
  for (i in 1:m) {
   G[a, 1, i] <- pnorm(e.pred[2], MUG[i, a], sdg)
   for (j in 2:(m-1))
    G[a, j, i] <- pnorm(e.pred[j+1], MUG[i, a], sdg) - pnorm(e.pred[j], MUG[i, a], sdg)
   G[a, m, i] <- 1-pnorm(e.pred[m], MUG[i, a], sdg)
  }
 }
 # Supervivencia
 SSS <- matrix(NA, ncol = Age.mature, nrow = M)
 for (a in 1:Age.mature)
  SSS[ ,a] <- predict(glmm.s, type = "response", newdata = data.frame(ln.h1 = X.pred, Age = rep(a-1, M)), re.form = NA)
 MUS <- matrix(NA, ncol = Age.mature, nrow = m)
 for (a in 1:(Age.mature)) {
  init.s <- 1
  for (i in 1:m) {
     MUS[i, a] <- mean(SSS[init.s:(init.s+sM-1), a])
     init.s <- init.s + sM
  }
 }
 # Probabilidad de reproducirse
 FF1 <- matrix (NA, nrow = M, ncol = Age.mature)
 for (a in 1:Age.mature)
  FF1[ ,a] <- predict(glmm.f1, type = "response", newdata = data.frame(ln.h1 = X.pred, Age = rep(a-1, M)), re.form = NA)
 MUF1 <- matrix(NA, ncol = Age.mature, nrow = m)
 for (a in 1:Age.mature) {
    init.f1 <- 1
    for (i in 1:m) {
     MUF1[i, a] <- mean(FF1[init.f1:(init.f1+sM-1), a])
     init.f1 <- init.f1 + sM
  }
 }
 # Número de frutos que produce un individuo
 FF2 <- matrix (NA, nrow = M, ncol = Age.mature)
 for (a in 1:Age.mature)
  FF2[ ,a] <- predict(gam.f2$gam, type = "response", newdata = data.frame(ln.h1 = X.pred, Age = rep(a-1, M)), re.form = NA)
 MUF2 <- matrix(NA, ncol = Age.mature, nrow = m)
 for (a in 1:Age.mature) {
    init.f2 <- 1
    for (i in 1:m) {
     MUF2[i, a] <- mean(FF2[init.f2:(init.f2+sM-1), a])
     init.f2 <- init.f2 + sM
  }
 }
 # Número de semillas por fruto que produce un individuo
 FF3 <- matrix (NA, nrow = M, ncol = Age.mature)
 for (a in 1:Age.mature)
  FF3[ ,a] <- predict(glmm.f3, type = "response", newdata = data.frame(ln.h1 = X.pred, Age = rep(a-1, M)), re.form = NA)
 MUF3 <- matrix(NA, ncol = Age.mature, nrow = m)
 for (a in 1:Age.mature) {
    init.f3 <- 1
    for (i in 1:m) {
     MUF3[i, a] <- mean(FF3[init.f3:(init.f3+sM-1), a])
     init.f3 <- init.f3 + sM
  }
 }
 # Tamaño de los nuevos individuos
 ff5 <- den.f5$y
 FF5 <- matrix(NA, nrow = m, ncol = m)
 for (i in 1:m) {
  FF5[i,] <- ff5
 }
 Mlist <- list(G = G, MUS = MUS, MUF1 = MUF1, MUF2 = MUF2, MUF3 = MUF3)
 save(Mlist, file="Mlist.RData")
}
