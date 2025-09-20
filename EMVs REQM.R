library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(gridExtra)
library(grid)

# Função para gerar dados da distribuição Nakagami-m bivariada
gerar_dados_nakagami_bivariada <- function(m1, m2, omega, n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  V1 <- rgamma(n, shape = m1, rate = (m1 + m2)/omega)
  V2 <- rgamma(n, shape = m2, rate = (m1 + m2)/omega)
  Z1 <- V1
  Z2 <- V1 + V2
  Y1 <- sqrt(Z1)
  Y2 <- sqrt(Z2)
  return(list(Y1 = Y1, Y2 = Y2))
}

# Estimador EMV (numérico)
estimate_nakagami_biv <- function(Y1, Y2) {
  z1 <- Y1^2
  z2 <- Y2^2
  
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
  
  initial_params <- c(
    pmax(mean(z1)^2 / var(z1), 0.1),
    pmax((mean(z2)^2 / var(z2)) - mean(z1)^2 / var(z1), 0.1),
    pmax(mean(z2) / (mean(z1)^2 / var(z1) + (mean(z2)^2 / var(z2)) - mean(z1)^2 / var(z1)), 0.1)
  )
  
  result <- optim(par = initial_params, fn = log_likelihood_gamma, method = "L-BFGS-B", lower = rep(1e-6, 3), control = list(maxit = 2000))
  m1_est <- result$par[1]
  m2_est <- result$par[2]
  Omega_est <- result$par[3] * (m1_est + m2_est)
  
  return(list(m1 = m1_est, m2 = m2_est, Omega = Omega_est))
}

# C-Estimador (forma fechada)
estimadores_tipo_emv_nakagami_flex <- function(Y1, Y2) {
  n <- length(Y1)
  z1 <- Y1^2
  z2 <- Y2^2
  z2 <- pmax(z2, z1 + 1e-6)
  
  sum_log_z1 <- sum(log(z1))
  sum_log_z2 <- sum(log(z2))
  sum_z1_log_z1 <- sum(z1 * log(z1))
  sum_z2_log_z2 <- sum(z2 * log(z2))
  sum_z2 <- sum(z2)
  
  diff_z <- pmax(z2 - z1, 1e-10)
  frac_z2_log_z2 <- (z2 * log(z2)) / diff_z
  frac_z1_log_z1 <- (z1 * log(z1)) / diff_z
  
  C <- (1/n) * sum(frac_z2_log_z2) * sum_z2 * sum_log_z1 -
    sum(frac_z1_log_z1) * sum_z2_log_z2 -
    sum_log_z1 * sum_z2_log_z2
  
  D <- sum(frac_z1_log_z1) * sum_z2_log_z2 -
    (1/n) * sum(frac_z2_log_z2) * sum_z2 * sum_log_z1 +
    n * sum_z2_log_z2 +
    (1 + (1/n) * sum_log_z2) * sum_z2 * sum_log_z1
  
  m2_hat <- if (abs(C) < 1e-10 || !is.finite(C) || !is.finite(D)) {
    (mean(z2)^2) / var(z2) - (mean(z1)^2) / var(z1)
  } else {
    -D / C
  }
  
  m1_hat <- (1 / sum_log_z1) * ((m2_hat - 1) * sum(frac_z1_log_z1) - n)
  m1_hat <- max(min(m1_hat, 10), 0.5)
  m2_hat <- max(min(m2_hat, 10), 0.5)
  beta_hat <- sum_z2 / (n * (m1_hat + m2_hat))
  omega_hat <- beta_hat * (m1_hat + m2_hat)
  
  return(list(m1 = m1_hat, m2 = m2_hat, Omega = omega_hat))
}

# Método dos Momentos
metodo_momentos <- function(Y1, Y2) {
  z1 <- Y1^2
  z2 <- Y2^2
  m1_hat <- (mean(z1)^2) / var(z1)
  m_sum <- (mean(z2)^2) / var(z2)
  m2_hat <- max(m_sum - m1_hat, 0.5)
  m1_hat <- max(m1_hat, 0.5)
  omega_hat <- mean(z2)
  return(list(m1 = m1_hat, m2 = m2_hat, Omega = omega_hat))
}

# Função para simular RMSE
simular_rmse <- function(true_params, n_list, reps = 300, global_seed = 1241) {
  set.seed(global_seed)
  m1_t <- true_params[1]
  m2_t <- true_params[2]
  omega_t <- true_params[3]
  
  results <- map_dfr(n_list, function(n) {
    errors <- list(
      EMV_m1 = numeric(reps),
      EMV_m2 = numeric(reps),
      EMV_omega = numeric(reps),
      C_Estimador_m1 = numeric(reps),
      C_Estimador_m2 = numeric(reps),
      C_Estimador_omega = numeric(reps),
      Momentos_m1 = numeric(reps),
      Momentos_m2 = numeric(reps),
      Momentos_omega = numeric(reps)
    )
    
    for (i in 1:reps) {
      dados <- gerar_dados_nakagami_bivariada(m1_t, m2_t, omega_t, n)
      tryCatch({
        emv_numerico <- estimate_nakagami_biv(dados$Y1, dados$Y2)
        c_estimador <- estimadores_tipo_emv_nakagami_flex(dados$Y1, dados$Y2)
        momentos <- metodo_momentos(dados$Y1, dados$Y2)
        
        errors$EMV_m1[i] <- emv_numerico$m1 - m1_t
        errors$EMV_m2[i] <- emv_numerico$m2 - m2_t
        errors$EMV_omega[i] <- emv_numerico$Omega - omega_t
        
        errors$C_Estimador_m1[i] <- c_estimador$m1 - m1_t
        errors$C_Estimador_m2[i] <- c_estimador$m2 - m2_t
        errors$C_Estimador_omega[i] <- c_estimador$Omega - omega_t
        
        errors$Momentos_m1[i] <- momentos$m1 - m1_t
        errors$Momentos_m2[i] <- momentos$m2 - m2_t
        errors$Momentos_omega[i] <- momentos$Omega - omega_t
      }, error = function(e) {
        for (key in names(errors)) errors[[key]][i] <- NA
      })
    }
    
    rmse_values <- list(n = n)
    for (key in names(errors)) {
      rmse_values[[key]] <- sqrt(mean(errors[[key]]^2, na.rm = TRUE))
    }
    
    return(rmse_values)
  })
  
  return(results)
}

# Função para criar gráficos com tema cinza e grades brancas (geral para m1 e m2)
criar_grafico_tema_cinza <- function(df_rmse, parametro, eixo_y) {
  df_param <- df_rmse %>%
    select(n, contains(parametro)) %>%
    pivot_longer(cols = -n, names_to = "Metodo", values_to = "RMSE") %>%
    mutate(Metodo = gsub(paste0("_", parametro), "", Metodo),
           Metodo = case_when(
             Metodo == "EMV" ~ "EMV",
             Metodo == "C_Estimador" ~ "C-Estimador",
             Metodo == "Momentos" ~ "Método dos Momentos"
           ))
  
  ggplot(df_param, aes(x = n, y = RMSE, color = Metodo, linetype = Metodo, shape = Metodo)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3) +
    labs(x = "n", y = eixo_y) +
    scale_color_manual(values = c(
      "EMV" = "black",
      "C-Estimador" = "#800000",
      "Método dos Momentos" = "#7EC8E3"
    )) +
    scale_linetype_manual(values = c(
      "EMV" = "solid",
      "C-Estimador" = "solid",
      "Método dos Momentos" = "solid"
    )) +
    scale_shape_manual(values = c(
      "EMV" = 15,
      "C-Estimador" = 15,
      "Método dos Momentos" = 15
    )) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_blank(),
      legend.position = c(1, 1),  # Colocando a legenda no canto superior direito
      legend.justification = c(1, 1),  # Ajusta para o canto superior direito
      legend.title = element_blank(),
      legend.text = element_text(size = 11),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      legend.key = element_rect(fill = "white"),
      panel.background = element_rect(fill = "gray90", color = "gray90"),
      panel.grid.major = element_line(color = "white", linewidth = 0.5),
      panel.grid.minor = element_line(color = "white", linewidth = 0.3),
      axis.line = element_line(color = "white", size = 1)
    ) +
    scale_x_continuous(breaks = unique(df_rmse$n)) +
    ylim(0, NA)
}

# Função separada para o gráfico de Omega (diferenciar as curvas)
criar_grafico_omega <- function(df_rmse) {
  df_param_omega <- df_rmse %>%
    select(n, contains("omega")) %>%
    pivot_longer(cols = -n, names_to = "Metodo", values_to = "RMSE") %>%
    mutate(Metodo = gsub("_omega", "", Metodo),
           Metodo = case_when(
             Metodo == "EMV" ~ "EMV",
             Metodo == "C_Estimador" ~ "C-Estimador",
             Metodo == "Momentos" ~ "Método dos Momentos"
           ))
  
  ggplot(df_param_omega, aes(x = n, y = RMSE, color = Metodo, linetype = Metodo, shape = Metodo)) +
    geom_line(linewidth = 1.5) +
    geom_point(size = 4) +
    labs(x = "n", y = expression(REQM(hat(Omega)))) +
    scale_color_manual(values = c(
      "EMV" = "black",
      "C-Estimador" = "#800000",
      "Método dos Momentos" = "#7EC8E3"
    )) +
    scale_linetype_manual(values = c(
      "EMV" = "solid",
      "C-Estimador" = "dashed",  # Diferenciar as linhas de Omega
      "Método dos Momentos" = "dotted"
    )) +
    scale_shape_manual(values = c(
      "EMV" = 15,
      "C-Estimador" = 17,
      "Método dos Momentos" = 18
    )) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_blank(),
      legend.position = c(1, 1),  # Colocando a legenda no canto superior direito
      legend.justification = c(1, 1),  # Ajusta para o canto superior direito
      legend.title = element_blank(),
      legend.text = element_text(size = 11),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      legend.key = element_rect(fill = "white"),
      panel.background = element_rect(fill = "gray90", color = "gray90"),
      panel.grid.major = element_line(color = "white", linewidth = 0.5),
      panel.grid.minor = element_line(color = "white", linewidth = 0.3),
      axis.line = element_line(color = "white", size = 1)
    ) +
    scale_x_continuous(breaks = unique(df_rmse$n)) +
    ylim(0, NA)
}

# Parâmetros verdadeiros
true_params <- c(2.5, 3.0, 5.0)
n_list <- c(30, 60, 100, 150, 200, 250, 300)

# Simulação do RMSE
df_rmse <- simular_rmse(true_params, n_list, reps = 10000)

# Criar gráficos individuais
grafico_m1 <- criar_grafico_tema_cinza(df_rmse, "m1", expression(REQM(hat(m)[1])))
grafico_m2 <- criar_grafico_tema_cinza(df_rmse, "m2", expression(REQM(hat(m)[2])))
grafico_omega <- criar_grafico_omega(df_rmse)

# Exibir gráficos lado a lado com uma única legenda e título
grid.arrange(
  grafico_m1, grafico_m2, grafico_omega, 
  ncol = 3
)

# Salvar gráficos
ggsave("rmse_comparison_omega.pdf", 
       grid.arrange(
         grafico_m1, grafico_m2, grafico_omega, 
         ncol = 3
       ), 
       width = 15, height = 5, dpi = 300)
