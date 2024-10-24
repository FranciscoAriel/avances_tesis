# Lectura de Datos de influenza aviar

infl=read.csv("influenza.csv")
yt=ts(infl[,"total"],start = c(2020,1),frequency = 12)
n = length(yt)
p=1

## Valores iniciales: Semilla de la señal log(alpha_t)

library(statespacer)
yy = as.matrix(infl[,"total"])
fit <- statespacer(y = yy, local_level_ind = TRUE)
Qt=fit$system_matrices$Q$level

alfa_t = ts(fit$smoothed$level,frequency = 12,start = c(2020,1))
theta_ini = ts(log(fit$smoothed$level),frequency = 12,start = c(2020,1))

windows()
ts.plot(yt,theta_ini,alfa_t,col = c(1,2,3),lty = c(1,2,3),main = "Serie y señal filtrada")
legend("topleft",legend = c("Serie","log nivel","nivel"),col = c(1,2,3),lty = c(1,2,3))
abline(h = 0)
# Mi propuesta
# Distribución POISSON INFLADA CON CEROS
# Asumir probabilidades de que sea yt provenga de cero estructural como FIJAS
theta_0=theta_ini
zt <- as.numeric(yy == 0)

tmp=glm(zt~theta_0, family =binomial(link = "probit"))
pit = ts(predict.glm(tmp, type = "response"), frequency = 12, start = c(2020, 1))

## Función para obtener vector score
f1 <- function(theta) {
    unos = rep(1,n)
    den = -zt*(unos - pit)*exp(-exp(theta)+theta)
    nume = pit + (unos-pit)*exp(-exp(theta))
    parte1 = den/nume
    parte2 = (unos-zt)*(yy-exp(theta))
    sco = parte1 + parte2
    return(sco)
}
## función para obtener matriz hessiana
f2 = function(theta) {
    unos = rep(1, n)
    den = -zt * (unos - pit)*exp(theta-exp(theta))*((1-exp(theta))*(pit+(1-pit)*exp(-exp(theta)))+(1-pit)*exp(theta-exp(theta)))
    nume = (pit-(1-pit)*exp(-exp(theta)))^2
    parte1 = den / nume
    parte2 =(unos-zt)*exp(theta)
    res = parte1 - parte2
    hes = diag(n)
    diag(hes) <- res
    return(hes)
}
## Especificación de modelo según las diapositivas de la exposición
Zt = matrix(1,n,1)
Z = diag(n);
T = diag(n);
r = row(T); 
c = col(T);
idx = (r>=c);
T[idx] = 1;
R = diag(n);

Q = as.numeric(Qt)*diag(n);
## Valores iniciales
a1 <- alfa_t[1]
P1 <- 1000
P1_ = matrix(0,n,n);
P1_[1,1] = P1;
a1_ = matrix(0,n,1);
a1_[1,1] = a1;

## Varianza de la señal
Psi = Z %*% T %*% (P1_ + R %*% Q %*% t(R)) %*% t(T) %*% t(Z);
## MEdia de la señal
mu = Z %*% T %*% a1_;
Psi_inv = solve(Psi);
#Criterios de paro del algoritmo
dif = 1
tol = 0.01
nmax = 30
cont = 0
while(dif > tol & nmax > cont){
    score = f1(theta_0) - Psi_inv %*% (theta_0 - mu)
    A = f2(theta_0) - Psi_inv
    A_inv = solve(A) 
    theta_1 = theta_0 - A_inv%*%score
    cont = cont + 1;
    dif = max(abs(theta_1-theta_0))
    cat("i= ",cont,"Dif= ",dif,"\n")
    theta_0 = theta_1
}
## Varianza de la señal dados los datos
    A2 = f2(theta_1) - Psi_inv
    G = solve(-A2)
    var_g = diag(G)
li_theta = ts(theta_1 - 1.96 * sqrt(var_g), frequency = 12, start = c(2020, 1))
ls_theta = ts(theta_1+1.96*sqrt(var_g),frequency = 12,start = c(2020,1))
lambda_t = ts(exp(theta_1),frequency = 12,start = c(2020,1))
lambda_li = ts(exp(li_theta),frequency = 12,start = c(2020,1))
lambda_ls = ts(exp(ls_theta),frequency = 12,start = c(2020,1))
var_final = ts(var_g, frequency = 12, start = c(2020, 1))
logy = log(yt)
log0 = ts(rep(-Inf, n), frequency = 12, start = c(2020, 1))
log0[zt == TRUE] <- log(exp(-1))
#Gráficos
windows()
ts.plot(theta_1,li_theta,ls_theta,lty = c(1,2,2),main = "Señal suavizada")
points(log0, pch = 3)
points(logy)
legend("bottomright", legend = c("Señal", "Intervalo de confianza 95%", " o log(yt); + log(e^-1)"), lty = c(1, 2, NA))

## Valor esperado de Y dada la señal
Ey <- ts(lambda_t*(1-pit), frequency = 12, start = c(2020, 1))
Ey_li <- ts(lambda_li * (1 - pit), frequency = 12, start = c(2020, 1))
Ey_ls <- ts(lambda_ls * (1 - pit), frequency = 12, start = c(2020, 1))
windows()
ts.plot(Ey,Ey_li,Ey_ls, lty = c(1,2,2), main = "Serie con nivel estimado")
legend("topleft", legend = c("Nivel estimado", "Intervalo de confianza 95%"), lty = c(1, 2))
points(yt, col = "red")

