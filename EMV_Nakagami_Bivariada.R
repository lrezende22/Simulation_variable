
library(ggplot2)


gerar_dados_nakagami_bivariada <- function(m1, m2, omega, n, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  beta <- omega / (m1 + m2)
  V1 <- rgamma(n, shape = m1, scale = beta)
  V2 <- rgamma(n, shape = m2, scale = beta)
  Z1 <- V1
  Z2 <- V1 + V2
  Y1 <- sqrt(Z1)
  Y2 <- sqrt(Z2)
  return(list(Y1 = Y1, Y2 = Y2))
}

# Estimador de Máxima Verossimilhança para Nakagami-m Bivariada
estimate_nakagami_biv <- function(Y1, Y2) {
  z1 <- Y1^2
  z2 <- Y2^2
  
  # Função de verossimilhança para Gama bivariada
  log_likelihood_gamma <- function(params) {
    alpha1 <- params[1]
    alpha2 <- params[2]
    beta <- params[3]
    
    if (any(params <= 0)) return(Inf)
    
    n <- length(z1)
    term1 <- -n * (alpha1 + alpha2) * log(beta) - n * lgamma(alpha1) - n * lgamma(alpha2)
    term2 <- (alpha1 - 1) * sum(log(z1))
    term3 <- (alpha2 - 1) * sum(log(z2 - z1))
    term4 <- -sum(z2) / beta
    
    return(-(term1 + term2 + term3 + term4))
  }
  
  
  initial_alpha1 <- mean(z1)^2 / var(z1)
  initial_alpha2 <- (mean(z2)^2 / var(z2)) - initial_alpha1
  initial_beta <- mean(z2) / (initial_alpha1 + initial_alpha2)
  
  initial_params <- c(
    pmax(initial_alpha1, 0.1),
    pmax(initial_alpha2, 0.1),
    pmax(initial_beta, 0.1)
  )
  
  # Otimização
  result <- optim(
    par = initial_params,
    fn = log_likelihood_gamma,
    method = "L-BFGS-B",
    lower = rep(1e-6, 3),
    control = list(maxit = 2000)
  )
  
  
  m1_est <- result$par[1]
  m2_est <- result$par[2]
  Omega_est <- result$par[3] * (m1_est + m2_est)
  log_lik <- -result$value
  
  return(list(m1 = m1_est, m2 = m2_est, Omega = Omega_est, log_lik = log_lik))
}


calcular_rmse <- function(true_params, n_list, reps = 300, global_seed = 1241) {
  set.seed(global_seed)
  m1_t <- true_params[1]
  m2_t <- true_params[2]
  omega_t <- true_params[3]
  
  rmse_dict <- list(n = c())
  for (param in c("m1", "m2", "omega")) {
    rmse_dict[[paste0("EMV_", param)]] <- c()
  }
  
  for (n in n_list) {
    erros <- list()
    
    for (i in 1:reps) {
      dados <- gerar_dados_nakagami_bivariada(m1_t, m2_t, omega_t, n)
      estimates <- estimate_nakagami_biv(dados$Y1, dados$Y2)
      
      # Calculando o erro para cada parâmetro
      erros$EMV_m1 <- c(erros$EMV_m1, estimates$m1 - m1_t)
      erros$EMV_m2 <- c(erros$EMV_m2, estimates$m2 - m2_t)
      erros$EMV_omega <- c(erros$EMV_omega, estimates$Omega - omega_t)
    }
    
    rmse_dict$n <- c(rmse_dict$n, n)
    
    # Calculando RMSE (REQM) para cada parâmetro
    for (param in c("EMV_m1", "EMV_m2", "EMV_omega")) {
      rmse_dict[[param]] <- c(rmse_dict[[param]], sqrt(mean(erros[[param]]^2)))
    }
  }
  
  return(rmse_dict)
}


true_params <- c(2.0, 1.5, 3.0)
n_list <- c(30, 40, 50, 100, 150, 200, 250, 300,400).

rmse_dict <- calcular_rmse(true_params, n_list, reps = 10000)

# Convertendo para DataFrame 
df_rmse <- data.frame(rmse_dict)





# ------------ Gerando gráficos separados para cada parâmetro ------------------



ggplot(df_rmse, aes(x = n, y = EMV_m1)) +
  geom_line(color = "#000080") +
  geom_point(color = "#000080") +
  labs(x = "Tamanho da amostra (n)", y = "RMSE para m1", title = "RMSE para m1 (EMV)") +
  theme_minimal()


ggplot(df_rmse, aes(x = n, y = EMV_m2)) +
  geom_line(color = "#000080") +
  geom_point(color = "#000080") +
  labs(x = "Tamanho da amostra (n)", y = "RMSE para m2", title = "RMSE para m2 (EMV)") +
  theme_minimal()


ggplot(df_rmse, aes(x = n, y = EMV_omega)) +
  geom_line(color = "#000080") +
  geom_point(color = "#000080") +
  labs(x = "Tamanho da amostra (n)", y = "RMSE para Omega", title = "RMSE para Omega (EMV)") +
  theme_minimal()

