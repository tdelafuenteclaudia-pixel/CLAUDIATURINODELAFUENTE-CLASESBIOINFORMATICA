###############################################################################
# TRABAJO R
###############################################################################

# 1. Cargar y explorar los datos
# Coloca archivo-entrada.txt en el directorio de trabajo de R.
datos <- read.table("archivo-entrada.txt",
                    header = TRUE, sep = "\t",
                    check.names = FALSE, stringsAsFactors = FALSE)

head(datos)
tail(datos)
dim(datos)

# Comprobacion: deben aparecer 619 filas (genes) y 7 columnas
# (nombre del gen + seis muestras).

control <- c("C2-TPM", "C1-TPM", "C3-TPM")
condicion <- c("1W5-TPM", "1W6-TPM", "1W5-COV")

###############################################################################
# 2. Calcular medias y desviaciones estandar por gen
###############################################################################

datos$Media_Control <- apply(datos[, control], 1, mean, na.rm = TRUE)
datos$DE_Control <- apply(datos[, control], 1, sd, na.rm = TRUE)
datos$Media_Condicion <- apply(datos[, condicion], 1, mean, na.rm = TRUE)
datos$DE_Condicion <- apply(datos[, condicion], 1, sd, na.rm = TRUE)

head(datos)

###############################################################################
# 3. Histogramas de las medias
###############################################################################

par(mfrow = c(1, 2), mar = c(5, 4, 4, 1))
hist(datos$Media_Control, breaks = 30, col = "steelblue",
     main = "Media de expresion - Control",
     xlab = "Media TPM", ylab = "Frecuencia")
hist(datos$Media_Condicion, breaks = 30, col = "darkorange",
     main = "Media de expresion - Condicion",
     xlab = "Media TPM", ylab = "Frecuencia")
par(mfrow = c(1, 1))

###############################################################################
# 4. Boxplot de las medias
###############################################################################

boxplot(datos$Media_Condicion, datos$Media_Control,
        names = c("Condicion", "Control"),
        col = c("darkorange", "steelblue"),
        main = "Distribucion de la expresion genica",
        ylab = "Media TPM")
legend("topleft", legend = c("Condicion", "Control"),
       fill = c("darkorange", "steelblue"), bty = "n")

###############################################################################
# 5. Normalidad y comparacion estadistica gen a gen
# Se aplica Shapiro-Wilk por gen y grupo. Si ambos grupos cumplen p > 0.05,
# se utiliza t de Welch (no asume varianzas iguales). Si alguno no cumple
# normalidad, se utiliza Wilcoxon-Mann-Whitney.
###############################################################################

shapiro_p <- function(x) {
  # Shapiro no puede evaluarse si todos los valores son iguales.
  if (length(unique(x[!is.na(x)])) < 2) return(NA_real_)
  tryCatch(shapiro.test(x)$p.value, error = function(e) NA_real_)
}

p_control <- numeric(nrow(datos))
p_condicion <- numeric(nrow(datos))
p_comparacion <- numeric(nrow(datos))
test_utilizado <- character(nrow(datos))

for (i in seq_len(nrow(datos))) {
  x <- as.numeric(datos[i, control])
  y <- as.numeric(datos[i, condicion])
  
  p_control[i] <- shapiro_p(x)
  p_condicion[i] <- shapiro_p(y)
  
  if (!is.na(p_control[i]) && !is.na(p_condicion[i]) &&
      p_control[i] > 0.05 && p_condicion[i] > 0.05) {
    prueba <- t.test(y, x, var.equal = FALSE)
    test_utilizado[i] <- "Welch t-test"
  } else {
    prueba <- wilcox.test(y, x, exact = FALSE)
    test_utilizado[i] <- "Wilcoxon"
  }
  p_comparacion[i] <- prueba$p.value
}

datos$Shapiro_Control <- p_control
datos$Shapiro_Condicion <- p_condicion
datos$Test_utilizado <- test_utilizado
datos$pvalue <- p_comparacion

###############################################################################
# 6. Genes diferencialmente expresados (p-value < 0.05)
###############################################################################

DEG <- datos[!is.na(datos$pvalue) & datos$pvalue < 0.05, ]

DEG$Expresion <- ifelse(DEG$Media_Condicion > DEG$Media_Control,
                        "Sobreexpresado", "Reprimido")

cat("Numero de genes diferencialmente expresados:", nrow(DEG), "\n")
cat("Genes sobreexpresados:",
    sum(DEG$Expresion == "Sobreexpresado"), "\n")
cat("Genes reprimidos:",
    sum(DEG$Expresion == "Reprimido"), "\n")

# Mostrar el listado de genes y su direccion de cambio.
print(DEG[, c("Gene Name", "Media_Control", "Media_Condicion",
              "pvalue", "Expresion")])

###############################################################################
# 7. Heatmap de los genes diferencialmente expresados
# Se transforma la expresion con log2(TPM + 1) y se estandariza cada gen
# (z-score por filas) para poder comparar el patron entre muestras.
###############################################################################

if (nrow(DEG) > 0) {
  muestras <- c(control, condicion)
  matriz <- as.matrix(DEG[, muestras])
  rownames(matriz) <- DEG$`Gene Name`
  
  matriz_log <- log2(matriz + 1)
  matriz_z <- t(scale(t(matriz_log)))
  matriz_z[is.na(matriz_z)] <- 0
  
  colores_muestras <- c(rep("steelblue", length(control)),
                        rep("darkorange", length(condicion)))
  
  heatmap(matriz_z,
          scale = "none",
          ColSideColors = colores_muestras,
          col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
          margins = c(6, 6),
          main = "Heatmap de genes diferencialmente expresados",
          xlab = "Muestras", ylab = "Genes",
          cexRow = ifelse(nrow(DEG) > 60, 0.35, 0.65),
          cexCol = 0.8)
  legend("topright", legend = c("Control", "Condicion"),
         fill = c("steelblue", "darkorange"), bty = "n", xpd = TRUE)
} else {
  cat("No se puede generar el heatmap porque no hay genes con p-value < 0.05.\n")
}
