# Simulation_variable


### exercicio 1

### gerar valores da normal pelo método da inversa ###


n = 100

u = runif(1)


x = 0


inv_normal = function(u){
  
  #u = runif(1)
  
  x = qnorm(u)
  
  return(x)
  
}

inv_normal(u)



for(i in 1:n){
  
  
  u = runif(1)
  
  x[i] = qnorm(u)
  
  
}


hist(x, xlab="Amostra", ylab="Frequência", border="black", main="")



### --- exercicio 2 ---


li = 1

ls = 2


fli = pnorm(li) # pnorm: Calcula a função de distribuição acumulada 

fls = pnorm(ls)



x = numeric(n)

inv_normal2 = function(u){
  
  x = qnorm(u, fli, fls)
  
  return(x)
  
}

#testando a função

inv_normal2(u=runif(1))


#iterações


for (i in 1:n){
  
  u = runif(1, fli, fls)  #uniforme(F(1), F(2)) para garantir o intervalo certo
  
  x[i] = qnorm(u)
  
}

hist(x, main="", xlab="Amostra", ylab="Frequência")



# Função densidade da Normal truncada
d_normal_trunc = function(x, li, ls) {
  dnorm(x) / (fls - fli)
}

hist(x, breaks = 30, probability = TRUE, ylab="Densidade",
     main = "",
     xlab = "Amostra", col = "gray", border = "black")

# Adicionando a curva de densidade teórica
curve(d_normal_trunc(x, li, ls), 
      col = "red", lwd = 2, add = TRUE)
legend("topright", legend = "Densidade Teórica", col="red", lwd=2, bty="0")




