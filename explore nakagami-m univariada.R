library(ggplot2)
library(scales)

# Função da densidade Nakagami-m
nakagami_pdf <- function(x, m, Omega) {
  return((2 / gamma(m)) * (m / Omega)^m * x^(2 * m - 1) * exp(-m / Omega * x^2))
}

# Função para o tema cinza com grades brancas
theme_gray_custom <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.background = element_rect(fill = "gray95", color = NA),
      panel.background = element_rect(fill = "gray90", color = "gray90"),
      panel.grid.major = element_line(color = "white", linewidth = 0.5),
      panel.grid.minor = element_line(color = "white", linewidth = 0.3),
      axis.line = element_line(color = "white", size = 1),
      plot.title = element_text(face = "bold", hjust = 0.5, color = "gray20"),
      plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
      axis.title = element_text(color = "gray30", family = "serif"),
      axis.text = element_text(color = "gray40", family = "serif"),
      legend.background = element_rect(fill = "white", color = "gray70"),
      legend.title = element_text(color = "gray30", family = "serif"),
      legend.text = element_text(color = "gray40", family = "serif"),
      strip.background = element_rect(fill = "gray85", color = "gray70"),
      strip.text = element_text(color = "gray30")
    )
}

# Valores de x para plotar
x <- seq(0, 5, length.out = 500)

# Gráfico 1: Variação do parâmetro de forma m (com Omega fixo)
Omega_fixo <- 2
m_values <- c(0.5, 1, 2, 4)
colors <- c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728')

df_m <- data.frame(x = rep(x, times = length(m_values)), 
                   y = unlist(lapply(m_values, function(m) nakagami_pdf(x, m, Omega_fixo))),
                   m = rep(m_values, each = length(x)))

p1 <- ggplot(df_m, aes(x = x, y = y, color = factor(m), fill = factor(m))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 0, ymax = y), alpha = 0.2) +
  labs(x = "x", y = "f(x)", color = "m") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme_gray_custom() +
  guides(fill = guide_legend(title = "m"), color = guide_legend(title = "m"))

# Gráfico 2: Variação do parâmetro de escala Omega (com m fixo)
m_fixo <- 2
Omega_values <- c(1, 2, 3, 4)

df_Omega <- data.frame(x = rep(x, times = length(Omega_values)), 
                       y = unlist(lapply(Omega_values, function(Omega) nakagami_pdf(x, m_fixo, Omega))),
                       Omega = rep(Omega_values, each = length(x)))

p2 <- ggplot(df_Omega, aes(x = x, y = y, color = factor(Omega), fill = factor(Omega))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 0, ymax = y), alpha = 0.2) +
  labs(x = "x", y = "f(x)", color = expression(Ω)) +  # Usando LaTeX para Ω
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme_gray_custom() +
  guides(fill = guide_legend(title = expression(Omega)), color = guide_legend(title = expression(Omega)))

# Gráfico 3: Comparação com a distribuição Rayleigh (m = 1)
m_rayleigh <- 1
Omega_rayleigh <- 2
rayleigh_pdf <- (x / (Omega_rayleigh / 2)) * exp(-x^2 / (Omega_rayleigh))

df_comparacao <- data.frame(x = x, 
                            nakagami = nakagami_pdf(x, m_rayleigh, Omega_rayleigh), 
                            rayleigh = rayleigh_pdf)

p3 <- ggplot(df_comparacao, aes(x = x)) +
  geom_line(aes(y = nakagami, color = 'Nakagami (m=1)'), size = 1) +
  geom_ribbon(aes(ymin = 0, ymax = nakagami), fill = '#1f77b4', alpha = 0.2) +
  geom_line(aes(y = rayleigh, color = 'Rayleigh'), linetype = 'dashed', size = 1) +
  labs(x = "x", y = "f(x)", color = "") +
  scale_color_manual(values = c('#1f77b4', '#ff7f0e')) +
  theme_gray_custom() +
  guides(color = guide_legend(title = "Distribuições"))

# Gráfico 4: Efeito do fading (valores extremos de m)
m_values_extreme <- c(0.5, 1, 4, 15)
Omega_extreme <- 2

df_extreme <- data.frame(x = rep(x, times = length(m_values_extreme)), 
                         y = unlist(lapply(m_values_extreme, function(m) nakagami_pdf(x, m, Omega_extreme))),
                         m = rep(m_values_extreme, each = length(x)))

p4 <- ggplot(df_extreme, aes(x = x, y = y, color = factor(m), fill = factor(m))) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = 0, ymax = y), alpha = 0.2) +
  labs(x = "x", y = "f(x)", color = "m") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  theme_gray_custom() +
  guides(fill = guide_legend(title = "Fading"), color = guide_legend(title = "Fading"))

# Agrupar os gráficos em uma matriz 2x2
grid.arrange(p1, p2, p3, p4, ncol = 2)
