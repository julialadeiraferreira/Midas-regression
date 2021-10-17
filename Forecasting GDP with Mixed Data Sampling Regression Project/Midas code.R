#upload dictionaries


library(zoo)
library(stargazer)
library(xts)
library(gets)
library(aTSA)
library(forecast)
library(stats)
library(midasr)
library(tseries)
library(ggplot2)
library(dplyr)
library(scales)
library(quantmod)


# y data type: samples quarter 
y <- ts(data=teste_t$tri, frequency = 4, start=c(2003,2), end=c(2019,3))

# adf test: ho = has a unit root
adf_gdp <- adf.test(as.vector(y),k=0)
acf(y)
pacf(y)
nsdiffs(y) 

#AR(1) for real GDP (log dif, seasonal adj)
storage_arx1 <- ts(frequency = 4, start=c(2014,3), end=c(2019,3))
for(t in 1:21){ 
  model_arx1 <- arx(y[1:45+t], mc = TRUE, ar=1:1) 
  p_arx1 <- predict(model_arx1, newdata = y, n.ahead = 1, se.fit = TRUE)
  storage_arx1[t] <- p_arx1
  print(t)
  print(model_arx1)
}  

e_arx1 = ((y[46:66] - storage_arx1)^(2))
RMSE_arx1 <- mean(e_arx1)
print(RMSE_arx1)



#Cumulative MSE
e_movel_arx1 <- ts(frequency = 4, start=c(2014,3), end=c(2019,3))
for (e in 1:21){
  e_movel_arx1[e] <- mean(e_arx1[1:e])
}
print(e_movel_arx1)


#Information Criteria matrix
comparacao <- matrix(data=0, nrow=28, ncol=5, dimnames = list(c("ARX1", dimnames(teste_m)[[2]][2:28])))
colnames(comparacao) = c("AIC","BIC","RMSE","DM statistics ARX1","DM pval/ARX1")
comparacao[1,1] <- AIC(model_arx1)
comparacao[1,2] <- BIC(model_arx1)
comparacao[1,3] <- mean(e_arx1)


# x data type: monthly 
x_tri <- ts(data=teste_t[,3:29], frequency = 4, start=c(2003,2), end=c(2019,3))

# adf test (loop for each variable)
storage_adf <- matrix(data=0, nrow=28, ncol=2, dimnames = list(c(dimnames(teste_m)[[2]][2:28], "GDP")))
for (i in 1:27){
  adf <-adf.test(as.vector(x_tri[,i]),k=0)
  print(adf)
  storage_adf[i,1] <-adf$statistic
  storage_adf[i,2] <-adf$p.value
}
print(storage_adf)  

# Recursive MIDAS (LOOPS for variables and recursiveness) <- forecast matrix(0:0,nrow=22,ncol=27)
# (log dif, seasonal adj)
m <- matrix(0:0,nrow=21,ncol=27)
storage_m <- ts(m, frequency = 4, start=c(2014,3), end=c(2019,3)) # matrix
dimnames(storage_m)[[2]] = dimnames(teste_m)[[2]][2:28] #header 
mean_storage_m <-matrix(0:0,nrow=26, ncol=1) # matrix
for (w in 1:27){
  x1 <- ts(data=teste_m[,1+w], frequency = 12, start=c(2003,4), end=c(2019,9))
  dimnames(x1)
  for (i in 1:21){
    x1t <- ts(x1[1:(132+3*i)],frequency = 12, start=c(2003,4))
    x1p <- ts(x1[(133+3*i):(198)],frequency = 12, start=c(2014,4+3*i))
    yt <- ts(y[1:(44+i)],frequency = 4, start=c(2003,2))
    model_m <- midas_r(yt ~  fmls(x1t, 12,3, nealmon), start = list(x1t = c(0.5,1,-0.5)))
    p_midas <- forecast(model_m, list(x1t = x1p), method = "static")
    p_midas_final <- as.ts(p_midas$mean[1])
    storage_m[i,w] <- p_midas_final
    comparacao[w+1,1] <- AIC(model_m)
    comparacao[w+1,2] <- BIC(model_m)
  } 
  e_midas = ((y[46:66] - storage_m[,w])^(2))
  mean_storage_m[w] <- mean(e_midas)
  comparacao[w+1,3] <- mean(e_midas)
}
tabela_uni <- as.matrix(mean_storage_m)
dimnames(tabela_uni) <- list(list=dimnames(teste_m)[[2]][2:28])  
print(tabela_uni)
print(comparacao)

##############pooling
p1 <- matrix(0:0,nrow=21,ncol=18)
storage_pool1 <- ts(p1, frequency = 4, start=c(2014,3), end=c(2019,3)) 
dimnames(storage_pool1)[[2]] = c("PIM g","PIM e","PIM t","PIM bk","PIM int","PIM bc","PIM bcd","PIM bcnd","PMC","PMC a","PMC a comb","PMC a alim","PMC a veic","IBC-BR","Cap Inst","BM","M1","M2") #header da matriz de previsão
storage_pool1[1:21,1:15] <- storage_m[1:21,1:15]
storage_pool1[1:21,16:18] <- storage_m[1:21,25:27] 

#pool 1/n
p2 <- matrix(0:0,nrow=21,ncol=1)
storage_pool2 <- ts(p2, frequency = 4, start=c(2014,3), end=c(2019,3)) 
for (p in 1:21){
  storage_pool2[p] <-mean(storage_pool1[p,1:18])
}

e_pool <-((y[46:66] - storage_pool2)^(2))
RMSE_pool <- mean(e_pool)
print(RMSE_pool)

mariano_pool <- dm.test(e_pool, e_arx1, alternative = "less", h=1, power =1)

#pool 1/RMSE
peso <- matrix(0:0,nrow=1,ncol=18)
dimnames(peso)[[2]] = c("PIM g","PIM e","PIM t","PIM bk","PIM int","PIM bc","PIM bcd","PIM bcnd","PMC","PMC a","PMC a comb","PMC a alim","PMC a veic","IBC-BR","Cap Inst","BM","M1","M2") #header da matriz de previsão
peso[,1:15] <- 1/tabela_uni[1:15,]
peso[,16:18] <- 1/tabela_uni[25:27,]
soma_peso <- sum(peso)
w <- peso/soma_peso

p3 <- matrix(0:0,nrow=21,ncol=18)
storage_pool3 <- ts(p3, frequency = 4, start=c(2014,3), end=c(2019,3)) 
dimnames(storage_pool3)[[2]] = c("PIM g","PIM e","PIM t","PIM bk","PIM int","PIM bc","PIM bcd","PIM bcnd","PMC","PMC a","PMC a comb","PMC a alim","PMC a veic","IBC-BR","Cap Inst","BM","M1","M2") #header da matriz de previsão
for (p in 1:21){
  storage_pool3[p,] <- w*storage_pool1[p,]
}

p4 <- matrix(0:0,nrow=21,ncol=1)
storage_pool4 <- ts(p2, frequency = 4, start=c(2014,3), end=c(2019,3)) 
for (p in 1:21){
  storage_pool4[p] <-sum(storage_pool3[p,1:18])
}

e_pool_w <-((y[46:66] - storage_pool4)^(2))
RMSE_pool_w <- mean(e_pool_w)
print(RMSE_pool_w)

mariano_pool_w <- dm.test(e_pool_w, e_arx1, alternative = "less", h=1, power =1)




# weights MIDAS

m <- matrix(0:0,nrow=21,ncol=27)
storage_m <- ts(m, frequency = 4, start=c(2014,3), end=c(2019,3)) # matriz de previsão
dimnames(storage_m)[[2]] = dimnames(teste_m)[[2]][2:28] #header da matriz de previsão
mean_storage_m <-matrix(0:0,nrow=27, ncol=1) # matriz RMSE
weights <- matrix(0:0,nrow=13, ncol=27)
colnames(weights) <- dimnames(teste_m)[[2]][2:28]
dimnames(weights) = list(c(0:12))
for (w in 1:27){
  x1 <- ts(data=teste_m[,1+w], frequency = 12, start=c(2003,4), end=c(2019,9))
  dimnames(x1)
  for (i in 1:21){
    x1t <- ts(x1[1:(132+3*i)],frequency = 12, start=c(2003,4))
    x1p <- ts(x1[(133+3*i):(198)],frequency = 12, start=c(2014,4+3*i))
    yt <- ts(y[1:(44+i)],frequency = 4, start=c(2003,2))
    model_m <- midas_r(yt ~  fmls(x1t, 12,3, nealmon), start = list(x1t = c(0.5,1,-0.5)))
    p_midas <- forecast(model_m, list(x1t = x1p), method = "static")
    p_midas_final <- as.ts(p_midas$mean[1])
    storage_m[i,w] <- p_midas_final
    comparacao[w+1,1] <- AIC(model_m)
    comparacao[w+1,2] <- BIC(model_m)
  } 
  e_midas = ((y[46:66] - storage_m[,w])^(2))
  mean_storage_m[w] <- mean(e_midas)
  comparacao[w+1,3] <- mean(e_midas)
  weights[,w] <- matrix(model_m$midas_coefficients[2:14]/model_m$coefficients[2])
}
tabela_uni <- as.matrix(mean_storage_m)
dimnames(tabela_uni) <- list(list=dimnames(teste_m)[[2]][2:28])  
print(tabela_uni)
print(comparacao)

weights[,9:13]
weights[,c(14:15,25:27)]

#####end Weights


# Cumulative MSE for each variable
y_matriz <- matrix(data=y[46:66], nrow=21, ncol=27)
erro_matriz <- ((y_matriz-storage_m)^2)

e <- matrix(0:0,nrow=21,ncol=27)
erro_movel_matriz <- ts(e, frequency = 4, start=c(2014,3), end=c(2019,3))
for(k in 1:27){
  for (j in 1:21){
    erro_movel_matriz[j,k] <- mean(erro_matriz[1:j,k]) 
  }
}

colnames(erro_movel_matriz) <- dimnames(teste_m)[[2]][2:28]
print(erro_movel_matriz)

# Cumulative MSE ratio
e_movel_arx1_matriz <- matrix(data=e_movel_arx1, nrow=21, ncol=27)
um <- ts(matrix(1, nrow=21, ncol=1), frequency = 4, start=c(2003,2), end=c(2019,3))


# Diebold and Mariano test (out-of-sample)
for (w in 1:27){
mariano_arx1 <- dm.test(erro_matriz[,w], e_arx1, alternative = "less", h=1, power =1)
comparacao[w+1,4] <- mariano_arx1$statistic
comparacao[w+1,5] <- mariano_arx1$p.value
}


#Graphs
razao_arx1 <- erro_movel_matriz/e_movel_arx1_matriz


graf1 <- plot(razao_arx1[,1:8], lty = 2, main = "MSE cumulative ratio - Industrial Production", ylab= "ratio", plot.type = "single", col = c("blue","pink", "red", "orange", "brown", "green", "purple", "springgreen4"))
lines(um, lty = 1, lwd=2, col = "black");
legend("topright", ncol =2, c(list=dimnames(teste_m)[[2]][2:9],"1"), 
       fill=c("blue","pink", "red", "orange", "brown", "green", "purple", "springgreen4", "black"))

graf2 <- plot(razao_arx1[,9:13], lty = 2, main = "MSE cumulative ratio - Retail Sales", ylab= "ratio", plot.type = "single", col = c("blue","pink", "red", "orange", "brown"))
lines(um, lty = 1, lwd=2, col = "black");
legend("topright", ncol =2, c(list=dimnames(teste_m)[[2]][10:14],"1"), 
       fill=c("blue","pink", "red", "orange", "brown", "black"))

graf3 <- plot(razao_arx1[,c(14:15, 25, 26, 27)], lty = 2, main = "MSE cumulative ratio - IBC-BR, Ind. Cap. and Mon. Agregates", ylab= "ratio", ylim = c(0,2.5), plot.type = "single", col = c("blue","pink", "red", "orange", "brown", "green", "purple"))
lines(um, lty = 1, lwd=2, col = "black");
legend("topright", ncol =2, c(list=dimnames(teste_m)[[2]][15:16],list=dimnames(teste_m)[[2]][26:28],"1"), 
       fill=c("blue","pink", "red", "orange", "brown", "black"))






