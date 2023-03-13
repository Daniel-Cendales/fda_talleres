# Cargue de librerías
library(fda)

ruta <- './Datos/datos.txt'
datos1 <- read.table(colClasses = rep(c('NULL', 'numeric'), 
                                      c(6, 1)*571), file = ruta)
# Gráficos
t_rango <- seq(275, 560, by = 0.5)
matplot(x = t_rango, y = t(datos1), type = 'p', 
        pch = 19, cex = 0.05, frame = FALSE, 
        main = 'Datos observados', col = 'gray',
        xlab = 'Emisión espectral', ylab = 'Fluoresencia')

# Base
base_bsp <- create.bspline.basis(breaks = c(seq(275, 389, length = 28.5), 
                                            470, 520, 560), norder = 4)

# Ajuste
ajuste <- Data2fd(y = as.matrix(t(datos1)), basisobj = base_bsp, 
                  argvals = t_rango)
lines(ajuste, frame = FALSE, main = 'Curvas Suavizadas', 
      ylab = 'Fluorescencia', xlab = 'Emisión espectral')


# Matriz de ajustes individuales
y <- t(unname(datos1))
Y_HAT <- smooth.basis(argvals = t_rango, y = t(datos1), 
                      fdParobj = base_bsp)

# Obtener los valores de la función estimada :)
eval.fd(evalarg = t_rango, Y_HAT$fd) -> f_hat

# 1.1 Función de media
y_media <- apply(X = f_hat, MARGIN = 1, mean)
lines(x = t_rango, y = y_media, col = 'black', lwd = 3)
legend(x = 'topright', legend = c('Media', 'Bandas conf.'), lwd = 3 , 
       col = c('black', 'red'), lty = c(1, 3))

# 1.2 Función de varianza
apply(X = f_hat, MARGIN = 1, sd) -> sd_fhat
matplot(x = t_rango, y = cbind(y_media - 1.96*sd_fhat, 
                               y_media + 1.96*sd_fhat),
        add = TRUE, type = 'l', col = 'red', lty = 3, lwd = 3)

fbplot(fit = f_hat, x = t_rango, ylim = c(-10, 500), frame = FALSE, 
       main = 'Boxplot funcional', xlim = range(t_rango)) -> bp_funcional

# Observaciones ordenadas de 
order(bp_funcional$depth) -> centrales
matlines(y = f_hat[, centrales[1:26]], x = t_rango, lwd = 2, lty = 1)

## Vamos a ver si medio puedo calcular la profundidad modificada
combn(x = seq(2:ncol(y), m = 2, FUN = function()))


