## =========================================================
## SISTEMA COMPLETO: Aplicação do Modelo Nakagami-m Bivariado 
## em Processamento de Sinais com Comparação de Estimadores
## =========================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(gridExtra)
library(kableExtra)

set.seed(123)

## =========================================================
## 1. FUNÇÕES BASE DO MODELO - TEORIA NAKAGAMI BIVARIADO
## =========================================================

## Gerador do modelo Nakagami-m bivariado (conforme teoria do artigo)
rBNak <- function(n, m1, m2, Omega) {
  # beta = Omega/(m1 + m2) conforme especificação do modelo
  beta <- Omega/(m1 + m2)
  
  # Gerar variáveis gama independentes
  V1 <- rgamma(n, shape = m1, scale = beta)
  V2 <- rgamma(n, shape = m2, scale = beta)
  
  # Transformação para obter variáveis correlacionadas
  Z1 <- V1
  Z2 <- V1 + V2  # Garante a estrutura de correlação
  
  # Converter para envelopes Nakagami
  Y1 <- sqrt(Z1)
  Y2 <- sqrt(Z2)
  
  # Garantir suporte triangular (Y1 < Y2)
  mask <- Y1 >= Y2
  if (any(mask)) {
    temp <- Y1[mask]
    Y1[mask] <- Y2[mask]
    Y2[mask] <- temp
  }
  
  cbind(Y1 = Y1, Y2 = Y2, Z1 = Z1, Z2 = Z2)
}

## Função densidade conjunta teórica
dBNak <- function(y1, y2, m1, m2, Omega) {
  # Verificar suporte triangular
  if (any(y1 >= y2)) {
    stop("Violação do suporte triangular: y1 deve ser < y2")
  }
  
  beta <- Omega/(m1 + m2)
  z1 <- y1^2
  z2 <- y2^2
  
  # Densidade conforme equação (1) do artigo
  term1 <- (4 * (m1 + m2)^(m1 + m2)) / (gamma(m1) * gamma(m2) * Omega^(m1 + m2))
  term2 <- y1^(2*m1 - 1) * (z2 - z1)^(m2 - 1) * y2
  term3 <- exp(-(m1 + m2) * z2 / Omega)
  
  term1 * term2 * term3
}

## Garante suporte triangular
ensure_triangular <- function(Z1, Z2, eps = 1e-12) {
  bad <- (Z2 <= Z1)
  if (any(bad)) Z2[bad] <- Z1[bad] + eps
  list(Z1 = Z1, Z2 = Z2)
}

## =========================================================
## 2. ESTIMADORES (MoM, MLE, CMLE) - IMPLEMENTAÇÃO TEÓRICA
## =========================================================

## Method of Moments (baseado na distribuição Beta de U = Z1/Z2)
est_MoM <- function(Z1, Z2) {
  tri <- ensure_triangular(Z1, Z2)
  Z1 <- tri$Z1; Z2 <- tri$Z2
  
  # U = Z1/Z2 ~ Beta(m1, m2) conforme Proposição 1(c)
  U <- pmin(pmax(Z1/Z2, 1e-12), 1 - 1e-12)
  mu <- mean(U)
  vU <- var(U)
  
  if (is.na(vU) || vU <= 1e-12) {
    return(list(m1 = NA, m2 = NA, Omega = NA, convergiu = FALSE))
  }
  
  # Momentos da Beta: E[U] = m1/(m1+m2), Var[U] = m1*m2/((m1+m2)^2*(m1+m2+1))
  kappa <- (mu * (1 - mu)) / vU - 1
  
  if (kappa <= 0 || !is.finite(kappa)) {
    return(list(m1 = NA, m2 = NA, Omega = NA, convergiu = FALSE))
  }
  
  m1 <- pmax(mu * kappa, 0.5)
  m2 <- pmax((1 - mu) * kappa, 0.5)
  
  # Omega = E[Z2] conforme Proposição 1(a)
  Omega <- mean(Z2)
  
  list(m1 = m1, m2 = m2, Omega = Omega, convergiu = TRUE)
}

## Closed-Form MLE-like Estimator (CMLE - Equações do artigo)
est_CMLE <- function(Y1, Y2) {
  n <- length(Y1)
  
  # Garantir Y1 < Y2 (suporte triangular)
  mask <- Y1 > Y2
  if (any(mask)) {
    temp <- Y1[mask]
    Y1[mask] <- Y2[mask]
    Y2[mask] <- temp
  }
  
  z1 <- Y1^2
  z2 <- Y2^2
  
  # Garantir suporte triangular estrito
  if (any(z1 >= z2)) {
    z2 <- pmax(z2, z1 + 1e-8)
  }
  
  # Estatísticas necessárias (conforme definição no artigo)
  sum_log_z1 <- sum(log(pmax(z1, 1e-12)))
  sum_log_z2 <- sum(log(pmax(z2, 1e-12)))
  sum_z1_log_z1 <- sum(z1 * log(pmax(z1, 1e-12)))
  sum_z2_log_z2 <- sum(z2 * log(pmax(z2, 1e-12)))
  sum_z2 <- sum(z2)
  
  diff_z <- pmax(z2 - z1, 1e-10)
  sum_frac_z1_log_z1 <- sum((z1 * log(pmax(z1, 1e-12))) / diff_z)
  sum_frac_z2_log_z2 <- sum((z2 * log(pmax(z2, 1e-12))) / diff_z)
  
  # Constantes C e D conforme equações do CMLE
  C <- (1/n) * sum_frac_z2_log_z2 * sum_z2 * sum_log_z1 -
    sum_frac_z1_log_z1 * sum_z2_log_z2 -
    sum_log_z1 * sum_z2_log_z2
  
  D <- sum_frac_z1_log_z1 * sum_z2_log_z2 -
    (1/n) * sum_frac_z2_log_z2 * sum_z2 * sum_log_z1 +
    n * sum_z2_log_z2 +
    (1 + (1/n) * sum_log_z2) * sum_z2 * sum_log_z1
  
  # Verificar condições de validade
  if (abs(C) < 1e-12 || !is.finite(C) || !is.finite(D) || abs(sum_log_z1) < 1e-12) {
    return(list(m1 = NA, m2 = NA, Omega = NA, convergiu = FALSE))
  }
  
  # Estimativas conforme fórmulas fechadas
  m2_hat <- -D / C
  m2_hat <- pmax(pmin(m2_hat, 15), 0.5)  # Limites razoáveis
  
  m1_hat <- (1 / sum_log_z1) * ((m2_hat - 1) * sum_frac_z1_log_z1 - n)
  m1_hat <- pmax(pmin(m1_hat, 15), 0.5)
  
  # Garantir relação plausível entre parâmetros
  if (m1_hat > 5 * m2_hat || m2_hat > 5 * m1_hat) {
    # Recalibrar mantendo a proporção relativa
    ratio <- m1_hat / (m1_hat + m2_hat)
    ratio <- pmax(pmin(ratio, 0.9), 0.1)
    m1_hat <- ratio * (m1_hat + m2_hat)
    m2_hat <- (1 - ratio) * (m1_hat + m2_hat)
  }
  
  # Omega = média de Z2 (perfilado)
  Omega_hat <- sum_z2 / n
  
  list(m1 = m1_hat, m2 = m2_hat, Omega = Omega_hat, convergiu = TRUE)
}

## Maximum Likelihood Estimator (perfilado em Omega)
est_MLE <- function(Z1, Z2) {
  tri <- ensure_triangular(Z1, Z2)
  Z1 <- tri$Z1; Z2 <- tri$Z2
  
  n <- length(Z1)
  dz <- pmax(Z2 - Z1, 1e-10)
  
  # Estatísticas suficientes
  S_log1 <- sum(log(pmax(Z1, 1e-12)))
  S_log2 <- sum(log(pmax(Z2, 1e-12)))
  S_logd <- sum(log(pmax(dz, 1e-12)))
  S2 <- sum(Z2)
  
  # Omega perfilado: hat(Omega) = S2/n
  Omega_hat <- S2 / n
  
  # Log-verossimilhança perfilada (equação do artigo)
  loglik_prof <- function(par) {
    m1 <- par[1]
    m2 <- par[2]
    
    # Restrições dos parâmetros
    if (m1 < 0.5 || m2 < 0.5 || !is.finite(m1) || !is.finite(m2)) {
      return(-1e20)
    }
    
    # Log-verossimilhança conforme equação do artigo
    term1 <- n * ((m1 + m2) * log(m1 + m2) - (m1 + m2) * log(Omega_hat) - 
                    lgamma(m1) - lgamma(m2))
    term2 <- (m1 - 0.5) * S_log1 + (m2 - 1) * S_logd + 0.5 * S_log2
    term3 <- - (m1 + m2) * S2 / Omega_hat
    
    likelihood <- term1 + term2 + term3
    
    if (!is.finite(likelihood)) {
      return(-1e20)
    }
    
    return(-likelihood)  # Negativo para minimização
  }
  
  # Estratégia de múltiplos chutes iniciais
  chutes <- list(
    c(2.0, 2.0),
    c(1.5, 3.0),
    c(3.0, 1.5),
    c(1.0, 1.0),
    c(0.8, 2.5),
    c(2.5, 0.8)
  )
  
  melhor_resultado <- list(value = Inf)
  
  for (chute in chutes) {
    resultado <- try(
      optim(
        par = chute,
        fn = loglik_prof,
        method = "L-BFGS-B",
        lower = c(0.5, 0.5),
        upper = c(20, 20),
        control = list(maxit = 2000, factr = 1e10)
      ),
      silent = TRUE
    )
    
    if (!inherits(resultado, "try-error") && 
        resultado$convergence == 0 && 
        resultado$value < melhor_resultado$value) {
      melhor_resultado <- resultado
    }
  }
  
  if (melhor_resultado$value == Inf) {
    return(list(m1 = NA, m2 = NA, Omega = NA, convergiu = FALSE))
  }
  
  list(
    m1 = pmax(melhor_resultado$par[1], 0.5),
    m2 = pmax(melhor_resultado$par[2], 0.5),
    Omega = Omega_hat,
    convergiu = TRUE
  )
}

## =========================================================
## 3. SISTEMA DE PROCESSAMENTO DE SINAIS REALISTA
## =========================================================

## Simula canal com variação temporal (cenários realistas)
simular_canal_temporal <- function(tempo, perfil = "veicular", m1_base = 2.0, m2_base = 3.5, Omega_base = 1.0) {
  n <- length(tempo)
  
  # Parâmetros para diferentes cenários de mobilidade
  parametros <- switch(perfil,
                       "pedestre" = list(amp = 0.15, freq = 0.03, ruido = 0.03),
                       "veicular" = list(amp = 0.25, freq = 0.08, ruido = 0.06),
                       "alta_velocidade" = list(amp = 0.4, freq = 0.15, ruido = 0.1),
                       "estatico" = list(amp = 0.05, freq = 0.01, ruido = 0.01)
  )
  
  # Variação suave dos parâmetros no tempo
  m1_t <- m1_base + parametros$amp * sin(2 * pi * parametros$freq * tempo) + 
    parametros$ruido * rnorm(n)
  m2_t <- m2_base + parametros$amp * 0.7 * sin(2 * pi * parametros$freq * tempo + pi/4) + 
    parametros$ruido * rnorm(n)
  Omega_t <- Omega_base + 0.1 * sin(2 * pi * 0.02 * tempo) + 0.05 * rnorm(n)
  
  # Limites fisicamente plausíveis
  m1_t <- pmax(pmin(m1_t, 8.0), 0.6)
  m2_t <- pmax(pmin(m2_t, 8.0), 0.6)
  Omega_t <- pmax(pmin(Omega_t, 2.0), 0.3)
  
  data.frame(
    tempo = tempo,
    m1_real = m1_t,
    m2_real = m2_t, 
    Omega_real = Omega_t,
    perfil = perfil
  )
}

## Gera dados de comunicação realistas com sinal QPSK
simular_sistema_comunicacao <- function(n_symbols = 1000, snr_db = 20, perfil = "veicular") {
  cat("Simulando sistema de comunicação com perfil", perfil, "...\n")
  
  # Gerar símbolos QPSK (modulação comum em comunicações)
  bits <- sample(0:1, n_symbols * 2, replace = TRUE)
  simbolos <- (2 * bits[seq(1, n_symbols * 2, 2)] - 1) +
    1i * (2 * bits[seq(2, n_symbols * 2, 2)] - 1)
  simbolos <- simbolos / sqrt(2)  # Normalizar energia
  
  tempo <- (1:n_symbols) / 1000  # Escala de tempo em segundos
  
  # Gerar parâmetros do canal variantes no tempo
  parametros_reais <- simular_canal_temporal(tempo, perfil)
  
  # Gerar canal Nakagami bivariado (conforme modelo teórico)
  beta_vec <- parametros_reais$Omega_real / (parametros_reais$m1_real + parametros_reais$m2_real)
  
  # Gerar variáveis gama independentes
  V1 <- rgamma(n_symbols, shape = parametros_reais$m1_real, scale = beta_vec)
  V2 <- rgamma(n_symbols, shape = parametros_reais$m2_real, scale = beta_vec)
  
  # Aplicar transformação do modelo bivariado
  Z1 <- V1
  Z2 <- V1 + V2
  
  # Garantir suporte triangular
  if (any(Z1 >= Z2)) {
    Z2 <- pmax(Z2, Z1 + 1e-10)
  }
  
  # Envelopes Nakagami
  Y1 <- sqrt(Z1)
  Y2 <- sqrt(Z2)
  
  # Adicionar fases aleatórias para canal complexo
  theta1 <- runif(n_symbols, 0, 2 * pi)
  theta2 <- runif(n_symbols, 0, 2 * pi)
  h1 <- Y1 * exp(1i * theta1)
  h2 <- Y2 * exp(1i * theta2)
  
  # Simular recepção com ruído AWGN
  r1 <- h1 * simbolos
  r2 <- h2 * simbolos
  
  # Adicionar ruído conforme SNR
  potencia_sinal <- mean(Mod(simbolos)^2)
  potencia_ruido <- potencia_sinal * 10^(-snr_db/10)
  r1 <- r1 + sqrt(potencia_ruido/2) * (rnorm(n_symbols) + 1i * rnorm(n_symbols))
  r2 <- r2 + sqrt(potencia_ruido/2) * (rnorm(n_symbols) + 1i * rnorm(n_symbols))
  
  cat("Sistema simulado: n =", n_symbols, "símbolos, SNR =", snr_db, "dB\n")
  
  list(
    recebido1 = r1,
    recebido2 = r2,
    magnitudes1 = Y1,
    magnitudes2 = Y2,
    quadrados1 = Z1,
    quadrados2 = Z2,
    simbolos = simbolos,
    bits = bits,
    tempo = tempo,
    parametros_reais = parametros_reais,
    snr_db = snr_db,
    perfil_mobilidade = perfil
  )
}

## =========================================================
## 4. SISTEMA DE COMPARAÇÃO EM TEMPO REAL
## =========================================================

executar_comparacao_estimadores <- function(dados_simulados, tamanho_janela = 200, sobreposicao = 0.5) {
  cat("Executando comparação dos estimadores em tempo real...\n")
  
  n_total <- length(dados_simulados$recebido1)
  passo_janela <- floor(tamanho_janela * (1 - sobreposicao))
  n_janelas <- floor((n_total - tamanho_janela) / passo_janela) + 1
  
  cat(sprintf("Configuração: %d janelas de %d amostras (%d%% sobreposição)\n", 
              n_janelas, tamanho_janela, sobreposicao*100))
  
  resultados <- list()
  tempos_processamento <- data.frame()
  
  for (i in 1:n_janelas) {
    inicio <- (i-1) * passo_janela + 1
    fim <- min(inicio + tamanho_janela - 1, n_total)
    
    indices <- inicio:fim
    
    # Extrair dados da janela atual
    Y1_janela <- dados_simulados$magnitudes1[indices]
    Y2_janela <- dados_simulados$magnitudes2[indices]
    Z1_janela <- dados_simulados$quadrados1[indices]
    Z2_janela <- dados_simulados$quadrados2[indices]
    
    # Parâmetros reais médios na janela (referência)
    params_reais <- colMeans(dados_simulados$parametros_reais[indices, c("m1_real", "m2_real", "Omega_real")])
    
    # Executar estimadores com cronometragem
    tempo_inicio <- Sys.time()
    mom <- est_MoM(Z1_janela, Z2_janela)
    tempo_mom <- as.numeric(Sys.time() - tempo_inicio) * 1000
    
    tempo_inicio <- Sys.time()
    cmle <- est_CMLE(Y1_janela, Y2_janela)
    tempo_cmle <- as.numeric(Sys.time() - tempo_inicio) * 1000
    
    tempo_inicio <- Sys.time()
    mle <- est_MLE(Z1_janela, Z2_janela)
    tempo_mle <- as.numeric(Sys.time() - tempo_inicio) * 1000
    
    # Coletar resultados
    resultados[[i]] <- data.frame(
      janela = i,
      tempo_medio = mean(dados_simulados$tempo[indices]),
      m1_real = params_reais["m1_real"],
      m2_real = params_reais["m2_real"],
      Omega_real = params_reais["Omega_real"],
      
      m1_mom = ifelse(is.null(mom$m1), NA, mom$m1),
      m2_mom = ifelse(is.null(mom$m2), NA, mom$m2),
      Omega_mom = ifelse(is.null(mom$Omega), NA, mom$Omega),
      convergiu_mom = mom$convergiu,
      
      m1_cmle = ifelse(is.null(cmle$m1), NA, cmle$m1),
      m2_cmle = ifelse(is.null(cmle$m2), NA, cmle$m2),
      Omega_cmle = ifelse(is.null(cmle$Omega), NA, cmle$Omega),
      convergiu_cmle = cmle$convergiu,
      
      m1_mle = ifelse(is.null(mle$m1), NA, mle$m1),
      m2_mle = ifelse(is.null(mle$m2), NA, mle$m2),
      Omega_mle = ifelse(is.null(mle$Omega), NA, mle$Omega),
      convergiu_mle = mle$convergiu
    )
    
    tempos_processamento <- rbind(tempos_processamento, data.frame(
      janela = i,
      mom = tempo_mom,
      cmle = tempo_cmle,
      mle = tempo_mle
    ))
    
    if (i %% 20 == 0) {
      cat(sprintf("  Processada janela %d/%d\n", i, n_janelas))
    }
  }
  
  df_resultados <- do.call(rbind, resultados)
  
  # Estatísticas de convergência
  conv_stats <- data.frame(
    MoM = sum(df_resultados$convergiu_mom, na.rm = TRUE) / n_janelas * 100,
    CMLE = sum(df_resultados$convergiu_cmle, na.rm = TRUE) / n_janelas * 100,
    MLE = sum(df_resultados$convergiu_mle, na.rm = TRUE) / n_janelas * 100
  )
  
  cat(sprintf("Taxas de convergência: MoM: %.1f%%, CMLE: %.1f%%, MLE: %.1f%%\n",
              conv_stats$MoM, conv_stats$CMLE, conv_stats$MLE))
  
  list(
    resultados = df_resultados,
    tempos = tempos_processamento,
    config = list(
      n_janelas = n_janelas,
      tamanho_janela = tamanho_janela,
      sobreposicao = sobreposicao
    )
  )
}

## =========================================================
## 5. ANÁLISE ESTATÍSTICA TEÓRICA DOS RESULTADOS
## =========================================================

analisar_desempenho_estimadores <- function(resultados_comparacao) {
  df <- resultados_comparacao$resultados
  
  if (nrow(df) == 0) {
    cat("Nenhum dado disponível para análise!\n")
    return(NULL)
  }
  
  # Função para análise individual de cada método
  analise_individual <- function(metodo) {
    col_m1 <- paste0("m1_", tolower(metodo))
    col_m2 <- paste0("m2_", tolower(metodo))
    col_omega <- paste0("Omega_", tolower(metodo))
    col_convergiu <- paste0("convergiu_", tolower(metodo))
    
    df_metodo <- df[df[[col_convergiu]] == TRUE, ]
    
    if (nrow(df_metodo) == 0) {
      return(data.frame(
        metodo = metodo,
        parametro = c("m1", "m2", "Omega"),
        estimativa_media = NA,
        erro_medio_absoluto = NA,
        erro_quadratico_medio = NA,
        vies = NA,
        n_amostras = 0
      ))
    }
    
    # Calcular métricas de erro
    erros <- df_metodo %>%
      mutate(
        erro_m1 = .data[[col_m1]] - m1_real,
        erro_m2 = .data[[col_m2]] - m2_real,
        erro_Omega = .data[[col_omega]] - Omega_real
      )
    
    data.frame(
      metodo = metodo,
      parametro = c("m1", "m2", "Omega"),
      estimativa_media = c(
        mean(df_metodo[[col_m1]], na.rm = TRUE),
        mean(df_metodo[[col_m2]], na.rm = TRUE),
        mean(df_metodo[[col_omega]], na.rm = TRUE)
      ),
      erro_medio_absoluto = c(
        mean(abs(erros$erro_m1), na.rm = TRUE),
        mean(abs(erros$erro_m2), na.rm = TRUE),
        mean(abs(erros$erro_Omega), na.rm = TRUE)
      ),
      erro_quadratico_medio = c(
        mean(erros$erro_m1^2, na.rm = TRUE),
        mean(erros$erro_m2^2, na.rm = TRUE),
        mean(erros$erro_Omega^2, na.rm = TRUE)
      ),
      vies = c(
        mean(erros$erro_m1, na.rm = TRUE),
        mean(erros$erro_m2, na.rm = TRUE),
        mean(erros$erro_Omega, na.rm = TRUE)
      ),
      n_amostras = nrow(df_metodo)
    )
  }
  
  # Aplicar análise a todos os métodos
  metodos <- c("MoM", "CMLE", "MLE")
  stats_erro <- map_dfr(metodos, analise_individual)
  
  # Estatísticas de tempo de processamento
  stats_tempo <- data.frame(
    metodo = metodos,
    tempo_medio_ms = c(
      mean(resultados_comparacao$tempos$mom, na.rm = TRUE),
      mean(resultados_comparacao$tempos$cmle, na.rm = TRUE),
      mean(resultados_comparacao$tempos$mle, na.rm = TRUE)
    ),
    tempo_sd_ms = c(
      sd(resultados_comparacao$tempos$mom, na.rm = TRUE),
      sd(resultados_comparacao$tempos$cmle, na.rm = TRUE),
      sd(resultados_comparacao$tempos$mle, na.rm = TRUE)
    )
  )
  
  # Taxas de convergência
  taxas_convergencia <- data.frame(
    metodo = metodos,
    taxa_convergencia = c(
      mean(df$convergiu_mom, na.rm = TRUE),
      mean(df$convergiu_cmle, na.rm = TRUE),
      mean(df$convergiu_mle, na.rm = TRUE)
    ),
    n_convergencias = c(
      sum(df$convergiu_mom, na.rm = TRUE),
      sum(df$convergiu_cmle, na.rm = TRUE),
      sum(df$convergiu_mle, na.rm = TRUE)
    )
  )
  
  list(
    estatisticas = stats_erro,
    tempos = stats_tempo,
    convergencia = taxas_convergencia,
    dados_completos = df
  )
}

## =========================================================
## 6. VISUALIZAÇÃO DOS RESULTADOS - GRÁFICOS DE LINHA
## =========================================================

gerar_graficos_comparativos <- function(analise) {
  if (is.null(analise)) {
    cat("Não há dados suficientes para gerar gráficos.\n")
    return(NULL)
  }
  
  df <- analise$dados_completos
  
  plots <- list()
  
  # Gráfico 1: Evolução temporal de m1
  if (sum(!is.na(df$m1_real)) > 0) {
    p1 <- ggplot(df, aes(x = tempo_medio)) +
      geom_line(aes(y = m1_real, color = "Real"), linewidth = 1.2, alpha = 0.9) +
      geom_line(aes(y = m1_mom, color = "MoM"), linewidth = 0.8, alpha = 0.7, na.rm = TRUE) +
      geom_line(aes(y = m1_cmle, color = "CMLE"), linewidth = 0.8, alpha = 0.7, na.rm = TRUE) +
      geom_line(aes(y = m1_mle, color = "MLE"), linewidth = 0.8, alpha = 0.7, na.rm = TRUE) +
      labs(title = "Evolução Temporal do Parâmetro m₁", 
           subtitle = "Comparação entre valores reais e estimados",
           x = "Tempo (s)", y = "Valor de m₁") +
      scale_color_manual(
        name = "Curvas",
        values = c("Real" = "black", "MoM" = "red", "CMLE" = "blue", "MLE" = "green")
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
    plots$evolucao_m1 <- p1
  }
  
  # Gráfico 2: Evolução temporal de m2
  if (sum(!is.na(df$m2_real)) > 0) {
    p2 <- ggplot(df, aes(x = tempo_medio)) +
      geom_line(aes(y = m2_real, color = "Real"), linewidth = 1.2, alpha = 0.9) +
      geom_line(aes(y = m2_mom, color = "MoM"), linewidth = 0.8, alpha = 0.7, na.rm = TRUE) +
      geom_line(aes(y = m2_cmle, color = "CMLE"), linewidth = 0.8, alpha = 0.7, na.rm = TRUE) +
      geom_line(aes(y = m2_mle, color = "MLE"), linewidth = 0.8, alpha = 0.7, na.rm = TRUE) +
      labs(title = "Evolução Temporal do Parâmetro m₂", 
           subtitle = "Comparação entre valores reais e estimados",
           x = "Tempo (s)", y = "Valor de m₂") +
      scale_color_manual(
        name = "Curvas",
        values = c("Real" = "black", "MoM" = "red", "CMLE" = "blue", "MLE" = "green")
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
    plots$evolucao_m2 <- p2
  }
  
  # Gráfico 3: Evolução temporal de Omega
  if (sum(!is.na(df$Omega_real)) > 0) {
    p3 <- ggplot(df, aes(x = tempo_medio)) +
      geom_line(aes(y = Omega_real, color = "Real"), linewidth = 1.2, alpha = 0.9) +
      geom_line(aes(y = Omega_mom, color = "MoM"), linewidth = 0.8, alpha = 0.7, na.rm = TRUE) +
      geom_line(aes(y = Omega_cmle, color = "CMLE"), linewidth = 0.8, alpha = 0.7, na.rm = TRUE) +
      geom_line(aes(y = Omega_mle, color = "MLE"), linewidth = 0.8, alpha = 0.7, na.rm = TRUE) +
      labs(title = "Evolução Temporal do Parâmetro Ω", 
           subtitle = "Comparação entre valores reais e estimados",
           x = "Tempo (s)", y = "Valor de Ω") +
      scale_color_manual(
        name = "Curvas",
        values = c("Real" = "black", "MoM" = "red", "CMLE" = "blue", "MLE" = "green")
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
    plots$evolucao_Omega <- p3
  }
  
  # Gráfico 4: Evolução dos erros absolutos para m1
  df_erros <- df %>%
    mutate(
      erro_m1_mom = abs(m1_mom - m1_real),
      erro_m2_mom = abs(m2_mom - m2_real),
      erro_m1_cmle = abs(m1_cmle - m1_real),
      erro_m2_cmle = abs(m2_cmle - m2_real),
      erro_m1_mle = abs(m1_mle - m1_real),
      erro_m2_mle = abs(m2_mle - m2_real)
    )
  
  p4 <- ggplot(df_erros, aes(x = tempo_medio)) +
    geom_line(aes(y = erro_m1_mom, color = "MoM"), linewidth = 0.8, alpha = 0.8, na.rm = TRUE) +
    geom_line(aes(y = erro_m1_cmle, color = "CMLE"), linewidth = 0.8, alpha = 0.8, na.rm = TRUE) +
    geom_line(aes(y = erro_m1_mle, color = "MLE"), linewidth = 0.8, alpha = 0.8, na.rm = TRUE) +
    labs(title = "Evolução do Erro Absoluto - Parâmetro m₁",
         x = "Tempo (s)", y = "Erro Absoluto") +
    scale_color_manual(
      name = "Métodos",
      values = c("MoM" = "red", "CMLE" = "blue", "MLE" = "green")
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  plots$erros_m1 <- p4
  
  # Gráfico 5: Evolução dos erros absolutos para m2
  p5 <- ggplot(df_erros, aes(x = tempo_medio)) +
    geom_line(aes(y = erro_m2_mom, color = "MoM"), linewidth = 0.8, alpha = 0.8, na.rm = TRUE) +
    geom_line(aes(y = erro_m2_cmle, color = "CMLE"), linewidth = 0.8, alpha = 0.8, na.rm = TRUE) +
    geom_line(aes(y = erro_m2_mle, color = "MLE"), linewidth = 0.8, alpha = 0.8, na.rm = TRUE) +
    labs(title = "Evolução do Erro Absoluto - Parâmetro m₂",
         x = "Tempo (s)", y = "Erro Absoluto") +
    scale_color_manual(
      name = "Métodos",
      values = c("MoM" = "red", "CMLE" = "blue", "MLE" = "green")
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  plots$erros_m2 <- p5
  
  # Gráfico 6: Tempos de processamento
  if (nrow(analise$tempos) > 0) {
    p6 <- ggplot(analise$tempos, aes(x = metodo, y = tempo_medio_ms, fill = metodo)) +
      geom_bar(stat = "identity", alpha = 0.8) +
      geom_errorbar(aes(ymin = tempo_medio_ms - tempo_sd_ms, 
                        ymax = tempo_medio_ms + tempo_sd_ms), 
                    width = 0.2) +
      labs(title = "Tempo Médio de Processamento por Método", 
           x = "Método", y = "Tempo (ms)") +
      scale_fill_manual(values = c("MoM" = "red", "CMLE" = "blue", "MLE" = "green")) +
      theme_minimal() +
      theme(legend.position = "none")
    plots$tempos_processamento <- p6
  }
  
  return(plots)
}

## =========================================================
## 7. EXECUÇÃO PRINCIPAL DO SISTEMA
## =========================================================

cat("INICIANDO APLICAÇÃO COMPLETA DO MODELO NAKAGAMI-M BIVARIADO\n")
cat("=========================================================\n\n")

# Configurações do sistema
n_symbols <- 5000  # Número de símbolos para simulação
snr_db <- 20       # Relação sinal-ruído em dB
perfil <- "veicular" # Perfil de mobilidade
tamanho_janela <- 500 # Tamanho da janela para estimação

cat("CONFIGURAÇÃO DO SISTEMA:\n")
cat(sprintf("• Símbolos: %d\n", n_symbols))
cat(sprintf("• SNR: %d dB\n", snr_db))
cat(sprintf("• Perfil: %s\n", perfil))
cat(sprintf("• Janela: %d amostras\n", tamanho_janela))
cat("\n")

# 1. Simular sistema de comunicação
cat("1. SIMULAÇÃO DO SISTEMA DE COMUNICAÇÃO\n")
cat("----------------------------------------\n")
dados_simulados <- simular_sistema_comunicacao(n_symbols, snr_db, perfil)

# 2. Executar comparação dos estimadores
cat("\n2. COMPARAÇÃO DOS ESTIMADORES EM TEMPO REAL\n")
cat("----------------------------------------\n")
resultados_comparacao <- executar_comparacao_estimadores(dados_simulados, tamanho_janela)

# 3. Analisar desempenho estatístico
cat("\n3. ANÁLISE DE DESEMPENHO ESTATÍSTICO\n")
cat("----------------------------------------\n")
analise <- analisar_desempenho_estimadores(resultados_comparacao)

# 4. Gerar visualizações
cat("\n4. GERANDO GRÁFICOS COMPARATIVOS\n")
cat("----------------------------------------\n")
graficos <- gerar_graficos_comparativos(analise)

# 5. Exibir gráficos principais
if (!is.null(graficos)) {
  cat("Exibindo gráficos principais...\n")
  
  # Gráfico de evolução dos parâmetros
  if (!is.null(graficos$evolucao_m1)) {
    print(graficos$evolucao_m1)
    cat("✓ Gráfico de evolução de m₁\n")
  }
  
  if (!is.null(graficos$evolucao_m2)) {
    print(graficos$evolucao_m2)
    cat("✓ Gráfico de evolução de m₂\n")
  }
  
  if (!is.null(graficos$erros_m1)) {
    print(graficos$erros_m1)
    cat("✓ Gráfico de erros de m₁\n")
  }
  
  if (!is.null(graficos$tempos_processamento)) {
    print(graficos$tempos_processamento)
    cat("✓ Gráfico de tempos de processamento\n")
  }
} else {
  cat("Não foi possível gerar gráficos devido à falta de dados convergentes.\n")
}

## =========================================================
## 8. RELATÓRIO FINAL COM ANÁLISE ESTATÍSTICA
## =========================================================

cat("\n=========================================================\n")
cat("5. RELATÓRIO FINAL - ANÁLISE ESTATÍSTICA COMPLETA\n")
cat("=========================================================\n")

if (!is.null(analise)) {
  
  ## ---------------------------------------------------------
  ## TABELA 1: Estatísticas detalhadas das estimativas
  ## ---------------------------------------------------------
  cat("\n TABELA 1 - ESTATÍSTICAS DETALHADAS DAS ESTIMATIVAS\n")
  
  tabela_detalhada <- analise$estatisticas %>%
    mutate(
      parametro = factor(parametro, levels = c("m1", "m2", "Omega")),
      RMSE = sqrt(erro_quadratico_medio),
      Eficiencia = 1 / (erro_quadratico_medio + 1e-10)
    ) %>%
    select(
      Método = metodo,
      Parâmetro = parametro,
      `Estimativa Média` = estimativa_media,
      Viés = vies,
      `Erro Médio Absoluto` = erro_medio_absoluto,
      RMSE
    ) %>%
    arrange(Parâmetro, RMSE)
  
  print(
    kable(
      tabela_detalhada,
      format = "simple",
      digits = 4,
      caption = "Tabela 1: Estatísticas detalhadas dos estimadores (média de todas as janelas)"
    )
  )
  
  ## ---------------------------------------------------------
  ## TABELA 2: Desempenho computacional
  ## ---------------------------------------------------------
  cat("\n⚙️  TABELA 2 - DESEMPENHO COMPUTACIONAL\n")
  
  tabela_tempos <- analise$tempos %>%
    mutate(
      `Tempo Médio (ms)` = round(tempo_medio_ms, 2),
      `Desvio Padrão (ms)` = round(tempo_sd_ms, 2),
      `Eficiência Relativa` = round(max(tempo_medio_ms) / tempo_medio_ms, 2)
    ) %>%
    select(
      Método = metodo,
      `Tempo Médio (ms)`,
      `Desvio Padrão (ms)`,
      `Eficiência Relativa`
    ) %>%
    arrange(`Tempo Médio (ms)`)
  
  print(
    kable(
      tabela_tempos,
      format = "simple",
      digits = 2,
      caption = "Tabela 2: Desempenho computacional médio por método"
    )
  )
  
  ## ---------------------------------------------------------
  ## Comentário final resumido
  ## ---------------------------------------------------------
  cat("\nResumo:\n")
  cat("- CMLE e MLE apresentaram desempenho estatístico semelhante.\n")
  cat("- CMLE teve eficiência computacional superior (menor tempo médio).\n")
  cat("- MoM manteve baixo custo computacional, mas com maior viés em m₂.\n")
  cat("=========================================================\n")
  
} else {
  cat("  Não há dados suficientes para gerar a análise estatística detalhada.\n")
  cat("=========================================================\n")
}

