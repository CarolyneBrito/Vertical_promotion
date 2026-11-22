pacman::p_load(foreign, tidyverse, survival, cmprsk, patchwork, scales)
dados_org <- read_csv("Banco de dados/dados_ossos.csv")

# Gráfico FIA GERAL ------
setwd("~/Estudos/Doutorado/Materias/Sobrevivência/Projeto Vertical Promotion/Artigo/Figuras")
cores <- c("black", "#d00000","#0349A5", "green" ,"orange","brown","yellow","purple")

pdf("fia_inicial.pdf", width = 8, height = 5)
par(mar = c(4.5, 3.5, 1.5, 1), oma = c(0, 0.1, 0, 0), mgp = c(1.9, 0.5, 0))
fia <- cuminc(ftime   = dados_org$time, fstatus = dados_org$status, cencode = 0)

plot.cuminc(fia, lty = 1, lwd = 2, ylim = c(0,0.7),
            cex.axis =1.4, cex.lab=1.4, cex.main = 1.6, cex = 1.4,  
            col = cores, xlab = "Tempo (dias)", ylab = "FIA",
            curvlab = c("Morte por câncer", "Morte por outras causas")) #8 5
dev.off()

# Gráfico FIA Variáveis Explicativas----
# Gráfico FIA Sexo
pdf("fia_graf_cov_1.pdf", width = 3800/300, height = 3000/300)
layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE))

par(mar = c(5, 3, 4, 1), oma = c(0, 0.1, 0, 0), mgp = c(1.9, 0.5, 0))

fia_sexo <- cuminc(ftime   = dados_org$time, fstatus = dados_org$status,
                   group   = dados_org$sexo_cat, cencode = 0)

#Causa 1
plot.cuminc(fia_sexo[c(1,2)], lty = 1, lwd = 2, ylim = c(0,0.7),
            cex.axis =1.4, cex.lab=1.4, cex.main = 1.6, cex = 1.4, col = cores,
            xlab = "Tempo (dias)", ylab = "FIA", main = "Morte por câncer",
            curvlab = c("Feminino", "Masculino"))
text(600, 0.7, "Sexo", font = 2, cex=1.4)
#Causa 2
plot.cuminc(fia_sexo[c(3,4)], lty = 1, lwd = 2, ylim = c(0,0.7),
            cex.axis =1.4, cex.lab=1.4, cex.main = 1.6, cex = 1.4, col = cores, 
            xlab = "Tempo (dias)", ylab = "FIA",
            curvlab = c("Feminino", "Masculino"), main = "Morte por outras causas")
text(600, 0.7, "Sexo", font = 2, cex=1.4)

# Gráfico FIA CATEATEND
fia_cat <- cuminc(ftime   = dados_org$time, fstatus = dados_org$status,
                  group   = dados_org$catend_cat, cencode = 0)

#Causa 1
plot.cuminc(fia_cat[c(1:2)], lty = 1, lwd = 2, ylim = c(0,0.7),
            cex.axis =1.4, cex.lab=1.4, cex.main = 1.6, cex = 1.4, col = cores, 
            xlab = "Tempo (dias)", ylab = "FIA",
            curvlab = c("Convênio/Particular", "SUS"),
            main = "Morte por câncer")
text(2150, 0.7, "Caráter de atendimento", font = 2, cex=1.4)
#Causa 2
plot.cuminc(fia_cat[c(3,4)], lty = 1, lwd = 2, ylim = c(0,0.7),
            cex.axis =1.4, cex.lab=1.4, cex.main = 1.6, cex = 1.4, col = cores, 
            xlab = "Tempo (dias)", ylab = "FIA",
            curvlab = c("Convênio/Particular", "SUS"),
            main = "Morte por outras causas")
text(2150, 0.7, "Caráter de atendimento", font = 2, cex=1.4)

dev.off()

# Gráfico FIA RADIO 
pdf("fia_graf_cov_2.pdf", width = 3800/300, height = 4150/300)
layout(matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, byrow = TRUE))

par(mar = c(5, 3, 4, 1), oma = c(0, 0.1, 0, 0), mgp = c(1.9, 0.5, 0))

fia_radio <- cuminc(ftime   = dados_org$time, fstatus = dados_org$status,
                    group   = dados_org$radio_cat, cencode = 0)

#Causa 1
plot.cuminc(fia_radio[c(1:2)], lty = 1, lwd = 2, ylim = c(0,0.85),
            cex.axis =1.4, cex.lab=1.4, cex.main = 1.6, cex = 1.4, 
            wh = c(0, 0.85), col = cores, 
            xlab = "Tempo (dias)", ylab = "FIA",
            curvlab = c("Com radioterapia", "Sem radioterapia"),
            main = "Morte por câncer")
text(1000, 0.85, "Radioterapia", font = 2, cex=1.4)
#Causa 2
plot.cuminc(fia_radio[c(3,4)], lty = 1, lwd = 2, ylim = c(0,0.85),
            cex.axis =1.4, cex.lab=1.4, cex.main = 1.6, cex = 1.4, col = cores, 
            xlab = "Tempo (dias)", ylab = "FIA",
            curvlab = c("Com radioterapia", "Sem radioterapia"),
            main = "Morte por outras causas")
text(1000, 0.85, "Radioterapia", font = 2, cex=1.4)

# Gráfico FIA QUIMIO
fia_quimio <- cuminc(ftime   = dados_org$time, fstatus = dados_org$status,
                     group   = dados_org$quimio_cat, cencode = 0)

#Causa 1
plot.cuminc(fia_quimio[c(1:2)], lty = 1, lwd = 2, ylim = c(0,0.75),
            cex.axis =1.4, cex.lab=1.4, cex.main = 1.6, cex = 1.4,  
            wh = c(0, 0.75), col = cores, 
            xlab = "Tempo (dias)", ylab = "FIA",
            curvlab = c("Com quimioterapia", "Sem quimioterapia"),
            main = "Morte por câncer")
text(1100, 0.75, "Quimioterapia", font = 2, cex=1.4)
#Causa 2
plot.cuminc(fia_quimio[c(3,4)], lty = 1, lwd = 2, ylim = c(0,0.75),
            cex.axis =1.4, cex.lab=1.4, cex.main = 1.6, cex = 1.4, col = cores, 
            xlab = "Tempo (dias)", ylab = "FIA",
            curvlab = c("Com quimioterapia", "Sem quimioterapia"),
  main = "Morte por outras causas")
text(1100, 0.75, "Quimioterapia", font = 2, cex=1.4)

# Gráfico FIA CIRURGIA  
fia_cirurgia <- cuminc(
  ftime   = dados_org$time,
  fstatus = dados_org$status,
  group   = dados_org$cirurgia_cat,
  cencode = 0)

#Causa 1
plot.cuminc(fia_cirurgia[c(1:2)], lty = 1, lwd = 2, ylim = c(0,0.85),
            cex.axis =1.4, cex.lab=1.4, cex.main = 1.6, cex = 1.4,  
            wh = c(0, 0.85), col = cores, 
            xlab = "Tempo (dias)", ylab = "FIA",
            curvlab = c("Com cirurgia","Sem cirurgia"),
            main = "Morte por câncer")
text(690, 0.85, "Cirurgia", font = 2, cex=1.4)
#Causa 2
plot.cuminc(fia_cirurgia[c(3,4)], lty = 1, lwd = 2, ylim = c(0,0.85),
            cex.axis =1.4, cex.lab=1.4, cex.main = 1.6, cex = 1.4, col = cores, 
            xlab = "Tempo (dias)", ylab = "FIA",
            curvlab = c("Com cirurgia", "Sem cirurgia"),
            main = "Morte por outras causas")
text(690, 0.85, "Cirurgia", font = 2, cex=1.4)
dev.off()

# Gráfico FIA RECIDIVA 
pdf("fia_graf_cov_3.pdf", width = 3800/300, height = 4150/300)
layout(matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, byrow = TRUE))

par(mar = c(5, 3, 4, 1), oma = c(0, 0.1, 0, 0), mgp = c(1.9, 0.5, 0))

fia_recidiva <- cuminc(ftime   = dados_org$time, fstatus = dados_org$status,
                       group   = dados_org$rec_cat, cencode = 0)

#Causa 1
plot.cuminc(fia_recidiva[c(1:2)], lty = 1, lwd = 2, ylim = c(0,0.9),
            cex.axis =1.4, cex.lab=1.4, cex.main = 1.6, cex = 1.4,  
            wh = c(0, 0.9), col = cores, 
            xlab = "Tempo (dias)", ylab = "FIA",
            curvlab = c("Com recidiva", "Sem recidiva"),
            main = "Morte por câncer")
text(1500, 0.9, "Recidiva do câncer", font = 2, cex=1.4)
#Causa 2
plot.cuminc(fia_recidiva[c(3,4)], lty = 1, lwd = 2, ylim = c(0,0.9),
            cex.axis =1.4, cex.lab=1.4, cex.main = 1.6, cex = 1.4, col = cores, 
            xlab = "Tempo (dias)", ylab = "FIA", 
            curvlab = c("Com recidiva", "Sem recidiva"),
            main = "Morte por outras causas")
text(1500, 0.9, "Recidiva do câncer", font = 2, cex=1.4)

# Gráfico FIA ECGRUP
fia_ec <- cuminc(ftime   = dados_org$time, fstatus = dados_org$status,
                 group   = dados_org$ecgrup_cat, cencode = 0)

#Causa 1
plot.cuminc(fia_ec[c(1:3)], lty = 1, lwd = 2, ylim = c(0,1),
            cex.axis =1.4, cex.lab=1.4, cex.main = 1.6, cex = 1.4,  
            wh = c(0, 1), col = cores, 
            xlab = "Tempo (dias)", ylab = "FIA",
            curvlab = c("I", "II/III", "IV"),
            main = "Morte por câncer")
text(1200, 1, "Estágio tumoral", font = 2, cex=1.4)
#Causa 2
plot.cuminc(fia_ec[c(4:6)], lty = 1, lwd = 2, ylim = c(0,1),
            cex.axis =1.4, cex.lab=1.4, cex.main = 1.6, cex = 1.4, col = cores, 
            xlab = "Tempo (dias)", ylab = "FIA", 
            curvlab = c("I", "II/III", "IV"),
            main = "Morte por outras causas")
text(1200, 1, "Estágio tumoral", font = 2, cex=1.4)

# Gráfico FIA M 
fia_m <- cuminc(ftime   = dados_org$time, fstatus = dados_org$status,
                group   = dados_org$m_cat, cencode = 0)

#Causa 1
plot.cuminc(fia_m[c(1:2)], lty = 1,lwd = 2, ylim = c(0,0.95),
            cex.axis =1.4, cex.lab=1.4, cex.main = 1.6, cex = 1.4, 
            wh = c(0, 0.95), col = cores, 
            xlab = "Tempo (dias)", ylab = "FIA",
            curvlab = c("Ausência de metástase", "Presença de metástase"),
            main = "Morte por câncer")
text(800, 0.95, "Metástase", font = 2, cex=1.4)
#Causa 2
plot.cuminc(fia_m[c(3:4)], lty = 1, lwd = 2, ylim = c(0,0.95),
            cex.axis =1.4, cex.lab=1.4, cex.main = 1.6, cex = 1.4, col = cores, 
            xlab = "Tempo (dias)", ylab = "FIA",
            curvlab = c("Ausência de metástase", "Presença de metástase"),
            main = "Morte por outras causas")
text(800, 0.95, "Metástase", font = 2, cex=1.4)

dev.off()

# Descritiva Categórica ----
vars_cat <- c("catend_cat", "ecgrup_cat", "m_cat", "radio_cat",
              "quimio_cat", "cirurgia_cat", "sexo_cat", "rec_cat")

nomes_cat <- c("Caráter de atendimento", "Estágio tumoral", "Metástase",
               "Radioterapia", "Quimioterapia", "Cirurgia", "Sexo",
               "Recidiva do câncer")

# FUNÇÃO PARA OS GRÁFICOS CATEGÓRICOS
grafico_cat <- function(dados, variavel, nome_grafico, limite_y) {
  dados %>%
    select(all_of(variavel)) %>%
    filter(!is.na(.data[[variavel]])) %>%
    count(.data[[variavel]], name = "n") %>%
    mutate(percentual = n / sum(n) * 100,
           categoria = as.character(.data[[variavel]])) %>%
    #Barras
    ggplot(aes(x = reorder(categoria, n), y = n)) +
    geom_col(width = 0.65, fill = "grey70", color = "grey25", linewidth = 0.4) +
    #Frequência absoluta e percentual
    geom_text(aes(label = paste0(n," (",round(percentual, 1),"%)")),
              hjust = -0.1, color = "black",size = 4, fontface = "plain") +
    coord_flip() +
    labs(title = nome_grafico, x = NULL, y = "Frequência") +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", size = 12,
                                    hjust = 0.5, color = "black"),
          axis.title = element_text(size = 11, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.text.x = element_text(size = 10, color = "black"),
          axis.line = element_line(linewidth = 0.5, color = "black"),
          axis.ticks = element_line(linewidth = 0.4, color = "black"),
          panel.grid = element_blank(),
          plot.margin = margin(8, 8, 8, 8)) +
    scale_y_continuous(limits = c(0, limite_y),
                       expand = expansion(mult = c(0, 0.03)))
}

limites_y <- c(3400, 1300, 2900, 3100, 2600, 2800, 2100, 2900)
graficos_cat <- Map(
  function(var, nome, limite) {
    grafico_cat(dados = dados_org, variavel = var, 
                nome_grafico = nome, limite_y = limite) },
  vars_cat, nomes_cat, limites_y)

names(graficos_cat) <- nomes_cat

fig_cat <- wrap_plots(graficos_cat, ncol = 2)
fig_cat

ggsave("descritiva_categoricas.pdf", fig_cat,
       width = 10, height = 10, units = "in", dpi = 300, bg = "white")

# TABELA DESCRITIVA DAS CATEGÓRICAS
tabela_cat <- dados_org %>%
  select(all_of(vars_cat)) %>%
  pivot_longer(cols = everything(), names_to = "Variavel",values_to = "Categoria") %>%
  filter(!is.na(Categoria)) %>%
  count(Variavel, Categoria, name = "n") %>%
  group_by(Variavel) %>%
  mutate(percentual = 100 * n / sum(n)) %>%
  ungroup()
tabela_cat

# Descritiva Quantitativa ----
vars_quant <- c("IDADE", "time")
nomes_quant <- c("Idade", "Tempo de acompanhamento")

dados_quant <- dados_org %>%
  select(all_of(vars_quant)) %>%
  pivot_longer(cols = everything(), names_to = "Variavel", values_to = "Valor") %>%
  mutate(Variavel = factor(Variavel, levels = vars_quant, labels = nomes_quant))

estatisticas <- dados_quant %>%
  group_by(Variavel) %>%
  summarise(media = mean(Valor, na.rm = TRUE),
            mediana = median(Valor, na.rm = TRUE),
            .groups = "drop")

#Histograma
grafico_hist <- ggplot(dados_quant, aes(x = Valor)) +
  #Histograma
  geom_histogram(bins = 20, fill = "grey75", color = "grey25",
                 linewidth = 0.45, alpha = 0.9) +
  # Linhas de média e mediana
  geom_vline(data = estatisticas, aes(xintercept = media, linetype = "Média"),
             color = "#2471A3",linewidth = 0.8) +
  geom_vline(data = estatisticas, aes(xintercept = mediana, linetype = "Mediana"),
             color = "#C0392B",linewidth = 0.8) +
  facet_wrap(~Variavel, scales = "free", ncol = 2) +
  labs(x = NULL, y = "Densidade", linetype = NULL) +
  scale_linetype_manual(values = c("Média" = "solid", "Mediana" = "dashed")) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12, color = "black"),
        axis.title = element_text(size = 11, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
        axis.line = element_line(linewidth = 0.5, color = "black"),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.text = element_text(size = 12),
        legend.key.width = unit(1.2, "cm"),
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        plot.margin = margin(8, 8, 5, 8))
grafico_hist

#Boxplot
grafico_box <- ggplot(dados_quant, aes(x = "", y = Valor)) +
  geom_boxplot(width = 0.35, fill = "grey85", color = "#2C3E50",
               linewidth = 0.8, outlier.shape = NA) +
  # Mediana
  stat_summary(fun = median, geom = "point", shape = 21, size = 2.5,
               fill = "#C0392B", color = "white", stroke = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.5,
               fill = "#2471A3", color = "white", stroke = 0.6) +
  facet_wrap(~Variavel, scales = "free", ncol = 2) +
  labs(x = NULL, y = "Valor") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12, color = "black"),
        axis.title = element_text(size = 11, color = "black"),
        axis.text = element_text(size = 11, color = "black"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line = element_line(linewidth = 0.5, color = "black"),
        panel.grid = element_blank(),
        plot.margin = margin(5, 8, 8, 8))

grafico_quant_completo <- grafico_box/grafico_hist +
  plot_layout(heights = c(0.85, 1.2))
grafico_quant_completo

ggsave("descritiva_quantitativas.pdf", grafico_quant_completo,
       width = 10, height = 6, units = "in", dpi = 300, bg = "white")

# TABELA DESCRITIVA DAS QUANTITATIVAS
tabela_quant <- dados_org %>%
  select(all_of(vars_quant)) %>%
  pivot_longer(cols = everything(), names_to = "Variavel", values_to = "Valor") %>%
  mutate(Variavel = factor(Variavel, levels = vars_quant, labels = nomes_quant)) %>%
  group_by(Variavel) %>%
  summarise(
    n = sum(!is.na(Valor)),
    Média = mean(Valor, na.rm = TRUE),
    Mediana = median(Valor, na.rm = TRUE),
    `Desvio-padrão` = sd(Valor, na.rm = TRUE),
    Mínimo = min(Valor, na.rm = TRUE),
    `1º Quartil` = quantile(Valor, 0.25, na.rm = TRUE),
    `3º Quartil` = quantile(Valor, 0.75, na.rm = TRUE),
    Máximo = max(Valor, na.rm = TRUE))
tabela_quant

# Descritiva Anual ----
grafico_ano <- dados_org %>%
  filter(!is.na(ano)) %>%
  count(ano) %>%
  ggplot(aes(x = ano, y = n)) +
  geom_line(linewidth = 0.8, color = "black") +
  geom_point(size = 2.5, color = "black") +
  geom_text(aes(label = n), vjust = -0.7, size = 4) +
  labs(x = "Ano", y = "Número de pacientes") +
  theme_classic() +
  theme(axis.title = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        axis.line = element_line(linewidth = 0.5, color = "black"),
        panel.grid = element_blank(),
        plot.margin = margin(8, 8, 8, 8)) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 10)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08)))

grafico_ano

grafico_ano_status <- dados_org %>%
  filter(!is.na(ano), !is.na(status_cat)) %>%
  count(ano, status_cat) %>%
  ggplot(aes(x = ano, y = n, group = status_cat, color = status_cat)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.3) +
  labs(x = "Ano",y = "Número de pacientes", color = NULL) +
  theme_classic() +
  theme(axis.title = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        legend.position = c(0.22, 0.90),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.background = element_rect(fill = "white", linewidth = 0.3),
        legend.key.width = unit(1.2, "cm"),
        panel.grid = element_blank()) +
  scale_color_manual(values = c("Censura" = "black",
                                "Morte por câncer" = "#0349A5",
                                "Morte por outras causas" = "#d00000")) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 10)) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0)))

grafico_ano_status

grafico_ano_completo <- grafico_ano +  grafico_ano_status 
grafico_ano_completo
ggsave("descritiva_ano.pdf", grafico_ano_completo, width = 13, height = 4.5,
       units = "in", dpi = 300, bg = "white")

# Descritiva Status ----
grafico_tempo_status <- dados_org %>%
  filter(!is.na(time), !is.na(status_cat)) %>%
  ggplot(aes(x = status_cat, y = time)) +
  # Boxplot
  geom_boxplot(width = 0.35, fill = "grey85", color = "#2C3E50",
               linewidth = 0.8, outlier.shape = NA) +
  # Mediana
  stat_summary(fun = median, geom = "point", shape = 21, size = 2.5,
               fill = "#C0392B", color = "white", stroke = 0.6) +
  # Média
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.5,
               fill = "#2471A3", color = "white", stroke = 0.6) +
  labs(x = "Status", y = "Tempo de acompanhamento") +
  theme_classic(base_size = 11) +
  theme(axis.title = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        axis.line = element_line(linewidth = 0.5, color = "black"),
        axis.ticks = element_line(linewidth = 0.4, color = "black"),
        panel.grid = element_blank(),
        plot.margin = margin(5, 8, 8, 8))

grafico_tempo_status

grafico_cat_s <- function(dados, variavel, nome_grafico, limite_y) {
  dados %>% 
    select(all_of(variavel)) %>%
    filter(!is.na(.data[[variavel]])) %>%
    count(.data[[variavel]], name = "n") %>%
    mutate(percentual = n / sum(n) * 100,
           categoria = as.character(.data[[variavel]])) %>%
    #Barras
    ggplot(aes(x = reorder(categoria, n), y = n)) +
    geom_col(width = 0.65, fill = "grey70", color = "grey25", linewidth = 0.4) +
    #Frequência absoluta e percentual
    geom_text(aes(label = paste0(n," (",round(percentual, 1),"%)")),
              hjust = 0.5, vjust = -0.5, color = "black",size = 4.5, 
              fontface = "plain") +
    labs(title = nome_grafico, x = "Status", y = "Frequência") +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", size = 12, 
                                    hjust = 0.5, color = "black"),
          axis.title = element_text(size = 15, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.text.x = element_text(size = 12, color = "black"),
          axis.line = element_line(linewidth = 0.5, color = "black"),
          axis.ticks = element_line(linewidth = 0.4, color = "black"),
          panel.grid = element_blank(),
          plot.margin = margin(8, 8, 8, 8)) +
    scale_y_continuous(limits = c(0, 1350), expand = expansion(mult = c(0, 0.03)))
}
grafico_status <- Map(
  function(var, nome, limite) {
    grafico_cat_s(dados = dados_org, variavel = var, 
                  nome_grafico = nome, limite_y = limite) },
  "status_cat", "", 1600)

grafico_status_completo <- grafico_status[[1]] +  grafico_tempo_status 

grafico_status_completo

ggsave("descritiva_status.pdf", grafico_status_completo,
       width = 13, height = 4.5, units = "in", dpi = 300, bg = "white")

# Comparação de grupos ----
# Sexo
cuminc(ftime = dados_org$time, fstatus = dados_org$status,
       group = dados_org$sexo_cat, cencode = 0)

# Catend
cuminc(ftime = dados_org$time, fstatus = dados_org$status,
       group = dados_org$catend_cat, cencode = 0)

# Radio
cuminc(ftime = dados_org$time, fstatus = dados_org$status,
       group = dados_org$radio_cat, cencode = 0)

#Quimio
cuminc(ftime = dados_org$time, fstatus = dados_org$status,
       group = dados_org$quimio_cat, cencode = 0)

#Cirurgia
cuminc(ftime = dados_org$time, fstatus = dados_org$status,
       group = dados_org$cirurgia_cat, cencode = 0)

#Recnenhum
cuminc(ftime = dados_org$time, fstatus = dados_org$status,
       group = dados_org$rec_cat, cencode = 0)

#M
cuminc(ftime = dados_org$time, fstatus = dados_org$status,
       group = dados_org$m_cat, cencode = 0)

#ECGroup
cuminc(ftime = dados_org$time, fstatus = dados_org$status,
       group = dados_org$ecgrup_cat, cencode = 0)

gray_pairwise <- function(dados, grupo1, grupo2) {
  dados_pair <- dados[dados$ecgrup_cat %in% c(grupo1, grupo2), ]
  dados_pair$ecgrup_cat <- factor(dados_pair$ecgrup_cat, levels = c(grupo1, grupo2))
  resultado <- cuminc(ftime = dados_pair$time, fstatus = dados_pair$status,
                      group = dados_pair$ecgrup_cat, cencode = 0)
  data.frame( Comparacao = paste(grupo1, "×", grupo2),
              p_causa1 = resultado$Tests[1, "pv"],
              p_causa2 = resultado$Tests[2, "pv"])
}

# Comparações
resultado_pairwise <- rbind(gray_pairwise(dados_org, "I", "II/III"),
                            gray_pairwise(dados_org, "I", "IV"),
                            gray_pairwise(dados_org, "II/III", "IV"))

resultado_pairwise$p_causa1_holm <- p.adjust(resultado_pairwise$p_causa1,
                                             method = "holm")

resultado_pairwise$p_causa2_holm <- p.adjust(resultado_pairwise$p_causa2,
                                             method = "holm")

resultado_pairwise