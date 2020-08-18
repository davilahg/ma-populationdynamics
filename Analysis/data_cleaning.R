
## AJUSTE BASE DE DATOS CON ESTIMACIÓN INVERSA DEL TAMAÑO
# Este código sirve para ajustar los datos y las tasas vitales tomando en cuenta el
# valor de supervivencia del primer registro de los individuos en la base de datos.
# Primero se necesita calcular la función g(), posteriormente usar esa función para
# calcular el tamaño de los individuos una unidad de tiempo antes de su primer registro
# y poder asignar ese valor de altura al de la supervivencia de la misma fila.
#

## LINUX
setwd("/home/gerardo/MEGA/Tesis/Datos") # datos
#setwd("/home/gerardo/MEGA/Tesis/Gráficas") # graficas
## WINDOWS
setwd("/Users/gdh/Documents/MEGA/Tesis/Datos")
#setwd("/Users/gdh/Documents/MEGA/Tesis/Gráficas")
## MAC
setwd("/Users/gdh/Google\ Drive/L\ 2017\ Gerardo\ Dávila/Datos")
#setwd("/Users/gdh/Google\ Drive/L\ 2017\ Gerardo\ Dávila/Gráficas")
## EDGAR
setwd("/Users/Edgar/Google\ Drive/Edgar/Trabajo/FC/Dirección\ de\ tesis/L\ 2017\ Gerardo\ Dávila/Datos")
#setwd("/Users/Edgar/Google\ Drive/Edgar/Trabajo/FC/Dirección\ de\ tesis/L\ 2017\ Gerardo\ Dávila/Gráficas")


# Ajustes a la base de datos
{
 Age.mature <- 100
 M <- 20000 # big-matrix dimension
 MimDB <- read.csv("MimDB.csv")
 MimDB <- MimDB[which(MimDB$Species == "Mimosa acantholoba"), ]
 isa.rows <- which(MimDB$Plot != "ISA") # remover los datos con parcelas ISA y ISX, no se sabe la edad
 isx.rows <- which(MimDB$Plot != "ISX")
 MimDB <- MimDB[isa.rows,]
 MimDB <- MimDB[isx.rows,]
 MimDB <- droplevels(MimDB)

 MimDB <- MimDB[!is.na(MimDB$Plot),]
 MimDB$Height <- as.numeric(as.character(MimDB$Height))
 MimDB <- transform(MimDB, id = as.factor(paste0("p", Plot, "i", Tag, "q", Qua)), Age.init = 0, Fst.cens = 2003, plot = NA)
 for (i in 1:nrow(MimDB)) {
	 if (MimDB$Plot[i] == "DIV") {MimDB$Fst.cens[i] <- 2005}
 }
 for (j in 1:nrow(MimDB)) {
	 if (MimDB$Plot[j] == "RIC") {MimDB$Age.init[j] <- 30} #1 Ricardo
	 else if (MimDB$Plot[j] == "ISP") {MimDB$Age.init[j] <- 23} #2 Isaìas Plano
   else if (MimDB$Plot[j] == "SEP") {MimDB$Age.init[j] <- 40} #3 Severino Plano = Severino I
	 else if (MimDB$Plot[j] == "ABE") {MimDB$Age.init[j] <- 10} #4 Abel
	 else if (MimDB$Plot[j] == "SEL") {MimDB$Age.init[j] <- Age.mature} #5 Severino Loma = Severino II
	 else if (MimDB$Plot[j] == "ISC") {MimDB$Age.init[j] <- 18} #6 Isaìas Carretera
	 else if (MimDB$Plot[j] == "BES") {MimDB$Age.init[j] <- 11} #7 Bernardo
	 else if (MimDB$Plot[j] == "TOA") {MimDB$Age.init[j] <- 0} #8 Toño Arriba = Toño de Leòn 0. ¡¡FUEGO!!
	 else if (MimDB$Plot[j] == "MAL") {MimDB$Age.init[j] <- 36} #9 Malaquìas
 	 else if (MimDB$Plot[j] == "HIL") {MimDB$Age.init[j] <- 5} #10 Hilario
 	 else if (MimDB$Plot[j] == "MAR") {MimDB$Age.init[j] <- 3} #11 Marciano
 	 else if (MimDB$Plot[j] == "FID") {MimDB$Age.init[j] <- 9} #12 Wilfrido
 	 else if (MimDB$Plot[j] == "ESS") {MimDB$Age.init[j] <- 16} #13 Escolar = Escolar Secundario
 	 else if (MimDB$Plot[j] == "TOB") {MimDB$Age.init[j] <- 1} #14 Toño Abajo = Toño de Leòn I
 	 else if (MimDB$Plot[j] == "DIA") {MimDB$Age.init[j] <- 7} #15 Dino agua caliente
 	 #else if (MimDB$Plot[j] == "JUL") {MimDB$Age.init[j] <- 60} #16 Juliàn    ¡¡ no hay Julián para Mimosa, hay para Lysiloma !!
   else if (MimDB$Plot[j] == "DIV") {MimDB$Age.init[j] <- 56} #17 Dino Vàlvula
 }
 for (j in 1:nrow(MimDB)) {
	 if (MimDB$Plot[j] == "RIC") {MimDB$plot[j] <- 1}
	 else if (MimDB$Plot[j] == "ISP") {MimDB$plot[j] <- 2}
   else if (MimDB$Plot[j] == "SEP") {MimDB$plot[j] <- 3}
	 else if (MimDB$Plot[j] == "ABE") {MimDB$plot[j] <- 4}
	 else if (MimDB$Plot[j] == "SEL") {MimDB$plot[j] <- 5}
	 else if (MimDB$Plot[j] == "ISC") {MimDB$plot[j] <- 6}
	 else if (MimDB$Plot[j] == "BES") {MimDB$plot[j] <- 7}
	 else if (MimDB$Plot[j] == "TOA") {MimDB$plot[j] <- 8}
	 else if (MimDB$Plot[j] == "MAL") {MimDB$plot[j] <- 9}
 	 else if (MimDB$Plot[j] == "HIL") {MimDB$plot[j] <- 10}
 	 else if (MimDB$Plot[j] == "MAR") {MimDB$plot[j] <- 11}
 	 else if (MimDB$Plot[j] == "FID") {MimDB$plot[j] <- 12}
 	 else if (MimDB$Plot[j] == "ESS") {MimDB$plot[j] <- 13}
 	 else if (MimDB$Plot[j] == "TOB") {MimDB$plot[j] <- 14}
 	 else if (MimDB$Plot[j] == "DIA") {MimDB$plot[j] <- 15}
 	 #else if (MimDB$Plot[j] == "JUL") {MimDB$plot[j] <- 16}
   else if (MimDB$Plot[j] == "DIV") {MimDB$plot[i] <- 17}
 }
 MimDB <- transform(MimDB, Age.census = 0)
 for (k in 1:nrow(MimDB)) {
	 MimDB$Age.census[k] <- MimDB$Census[k] - MimDB$Fst.cens[k] + MimDB$Age.init[k]
	 if (MimDB$Age.census[k] > Age.mature) {
	 	 MimDB$Age.census[k] <- Age.mature
	 }
 }
 MimDB <- transform(MimDB, TotFrut = 0, TotSeeds = 0)
 for (n in 1:nrow(MimDB)) {
	 if (MimDB$Rep[n] == 1) {
		 MimDB$TotFrut[n] <- MimDB$FrutPerBranch[n] * MimDB$NumbBranch[n]
		 MimDB$TotSeeds[n] <- round(MimDB$AvrgSeeds[n] * MimDB$TotFrut[n])
	 }
 }
}

# Crear Mim.s y Mim.g
{
 MimDB <- transform(MimDB, sup = rep(NA, nrow(MimDB)), h1 = rep(NA, nrow(MimDB)), h2 = rep(NA, nrow(MimDB)))
 MimDB$Sup[is.na(MimDB$Sup)] <- 0
 Mim.s <- as.data.frame(matrix(NA, nrow = 0, ncol = ncol(MimDB)))
 for (i in levels(MimDB$id)) {
 	Mim.id.s <- subset(MimDB, id == i)
 	if (nrow(Mim.id.s) > 1) {
 		Mim.id.s <- Mim.id.s[order(Mim.id.s$Census), ]
 		for (j in 1:nrow(Mim.id.s)) {
 			if (Mim.id.s$Sup[j] == 0) { # Después de que se registra un 0, descartar los otros datos. A veces hay individuos que se siguen registrando, por lo que se tienen muchos 0 por individuo y se infla la mortalidad.
 				Mim.id.s <- Mim.id.s[1:j, ]
     break
 			}
 		}
    if (nrow(Mim.id.s) > 1) {
  		Mim.id.s$sup[1:nrow(Mim.id.s)] <- Mim.id.s$Sup[1:nrow(Mim.id.s)]
  		Mim.id.s$h1[2:nrow(Mim.id.s)] <- Mim.id.s$Height[1:(nrow(Mim.id.s))-1]
      Mim.id.s$h2[1:nrow(Mim.id.s)] <- Mim.id.s$Height[1:nrow(Mim.id.s)] ## La última fila no debe tener h1 ni h2. No sirve para S ni para G. La penúltima sólo tiene h1. Sirve sólo para S. Para G se necesita hacer un subconjunto de S tal que s == 1.
 	 	  Mim.s <- rbind(Mim.s, Mim.id.s)
   }
  }
 }
 Mim.s <- transform(Mim.s, ln.h1 = log(h1), ln.h2 = log(h2), Age = Age.census)
 used.vars <- names(Mim.s[c("Age", "ln.h1", "ln.h2", "h1", "h2", "plot", "sup", "id", "Census")])
 can.s <- Mim.s[ ,used.vars]
 #write.table(Mim.s, file = "Mim.s.csv", sep = ",", row.names = FALSE)
 #write.table(can.s, file = "can.s.csv", sep = ",", row.names = FALSE)
}

# Ajuste datos sotobosque
{
 understory <- read.csv("Sotobosque1.csv")
 understory$Plot <- as.factor(as.character(understory$Plot))
 undtoa <- subset(understory, Plot == 8 | Plot == 8.1)
 undtoa <- subset(undtoa, Año < 2010)
 undtoa <- droplevels(undtoa)
 und <- subset(understory, Plot != 8 & Plot != 8.1)  #descartar los datos de la parcela TOA después del fuego (a partir de 2010)
 und <- droplevels(und)
 understory <- rbind(und, undtoa)
 understory <- transform(understory, id = as.factor(paste0("p", Plot, "i", Etiqueta)), ln.h1 = log(Altura), sup = rep(NA, nrow(understory)), h1 = rep(NA, nrow(understory)), h2 = rep(NA, nrow(understory)))
 understory <- subset(understory, Supv != "Exc" & Supv != "NE" )
 for (i in 1:nrow(understory))
    if (understory$Supv[i] == "Dosel") {understory$Supv[i] <- 1}
 understory <- droplevels(understory)
 und.00 <- as.data.frame(matrix(nrow = 0, ncol = ncol(understory)))
 for (i in levels(understory$id)) { # agregar un 0 a la supervivencia de los individuos que tienen puros 1s, agregarlo al final
  	und.0 <- subset(understory, id == i)
  	if (und.0$Supv[length(und.0$id)] == 1) {
      newrow <- as.data.frame(matrix(nrow = 1, ncol = ncol(und.0))) # agregar una fila
      names(newrow) <- names(und.0)
      und.01 <- rbind(und.0, newrow)
  		und.01[1] #agregar 0 al final de la columna de Supv, ahora hay que ver que si se murieron y revivieron todos los anteriores sean 1, intentar con elseif
  		und.01$Año[nrow(und.01)] <- und.01$Año[nrow(und.01)]+1
      und.01$Plot[nrow(und.01)] <- und.01$Plot[nrow(und.01)-1]
      und.01$id[nrow(und.01)] <- und.01$id[nrow(und.01)-1]
      und.01$Ecenso[nrow(und.01)] <- und.01$Ecenso[nrow(und.01)-1]+1
      und.01$Supv[nrow(und.01)] <- 0
  		print("1")
  	    }
  	else if (und.0$Supv[length(und.0$id)] == 0) {
  		und.01 <- und.0
  		print("0")
  	    }
   	und.00 <- rbind(und.00, und.01)
  }
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
  und.s <- transform(und.s, ln.h1 = log(h1), ln.h2 = log(h2), Age = Ecenso, Census = Año, plot = Plot) #All variables
  und.S <- und.s[,c("ln.h1", "ln.h2", "sup", "Age", "plot", "id", "Census", "h1", "h2")] #Only used variables
  #write.table(und.s, file = "und.s.csv", sep = ",", row.names = FALSE)
  #write.table(und.S, file = "und.S.csv", sep = ",", row.names = FALSE)
}
und.S <- read.csv("/home/gerardo/MEGA/Tesis/Datos/und.S.csv")
und.S <- read.csv("/Users/gdh/Google\ Drive/L\ 2017\ Gerardo\ Dávila/Datos/und1.S.csv")
# agregar a los individuos de edad = 0
oldund <- read.csv("Sotobosque.csv")
oldund <- subset(oldund, Age == 0)
oldund <- transform(oldund, id = as.factor(paste0("p", plot, "i", Tag)), h2 = h1)
oldund <- transform(oldund, h1 = rep(NA, nrow(oldund)), ln.h1 = rep(NA, nrow(oldund)), sup = rep(1, nrow(oldund)), ln.h2 = log(h2), plot = rep(8, nrow(oldund)))
oldund <- oldund[,c("ln.h1", "ln.h2", "sup", "Age", "plot", "id", "Census", "h1", "h2")]
und.S <- rbind(und.S, oldund)

# FUNCIÓN PARA CALCULAR h ANTES DEL PRIMER REGISTRO
 #Se necesita cargar "Growth" para correr lmm.g y utilizarlo para calcular glmm.s, usar script "Mim0.R"
{
 ## ln.h1 = altura antes de ser registrado por primera vez, ln.h2 = altura en su primer registro
 geth1 <- function(ln.h2, Age, Plot, Census) { ## Census? plot? id?
  b0 <- fx.coef[1] + pt.coef$int[which(pt.coef$plot == Plot)] + ifelse(Census != 0, cs.coef$int[which(cs.coef$Census == Census)], 0)
  b1 <- fx.coef[2] + pt.coef$h1[which(pt.coef$plot == Plot)] + ifelse(Census != 0, cs.coef$h1[which(cs.coef$Census == Census)], 0)
  b2 <- fx.coef[3]
  b3 <- fx.coef[4]
   ln.h1 <- (ln.h2-b0-b2*Age)/(b1+b3*Age) #calcular el modelo general
   return(ln.h1)
   }
}

# CALCULAR TAMAÑO ANTERIOR PARA PRIMER REGISTRO
{
 fx.coef <- fixef(lmm.g)
 pt.coef <- ranef(lmm.g)$plot
 pt.coef <- transform(pt.coef, plot = as.factor(row.names(pt.coef)))
 colnames(pt.coef) <- c("int", "h1", "plot")
 cs.coef <- ranef(lmm.g)$Census
 cs.coef <- transform(cs.coef, Census = as.factor(row.names(cs.coef)))
 colnames(cs.coef) <- c("int", "h1", "Census")
 censuses <- as.factor(row.names(cs.coef))
 newSurvival <- Survival
 for (i in levels(Survival$id)) {
   S.i <- subset(Survival, id == i)
   ln.h2 <- S.i$ln.h2[1]
   Age <- S.i$Age[1]
   Plot <- S.i$plot[1]
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
 for (p in i:nplot) {

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
