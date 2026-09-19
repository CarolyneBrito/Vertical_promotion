pacman::p_load(foreign, tidyverse, survival, cmprsk, patchwork, scales)

# Salvamento do banco em RDS ------
dados <- lapply(2000:2022, function(ano) {
  cat("Lendo ano:", ano, "\n") 
  read.dbf(paste0("RHC_GERAL", "/RHC_SP_", ano, ".dbf"), as.is = TRUE) %>% 
    mutate(ano = ano)
}) %>% 
  bind_rows()

# Limpeza do banco em RDS -----
dados_org<- dados %>% mutate(time = as.numeric(DTULTINFO - dados$DTDIAG),
                             status = case_when(ULTINFO %in% c(1,2) ~ 0,
                                                ULTINFO == 3 ~ 1,
                                                ULTINFO == 4 ~ 2)) %>% 
  filter(time >= 0) %>% 
  mutate(time = time+1) %>% 
  filter(TOPOGRUP %in% c("C40","C41")) %>% 
  select(time, status, 
         IDADE, SEXO, TOPO, TOPOGRUP, CATEATEND, RADIO, QUIMIO, CIRURGIA,
         NENHUM, RECNENHUM, PERDASEG, ECGRUP, M, ano) %>% 
  filter(ECGRUP %in% c("I","II","III","IV"),
         M %in% c("0", "1", "1A", "1B"))   %>% 
  mutate(sexo_cat = case_when(SEXO == 1 ~ "Masculino",
                              SEXO == 2 ~ "Feminino",
                              TRUE ~ NA_character_),
         idade_cat = case_when(IDADE <= 60 ~ "Até 60 anos",
                               IDADE > 60 ~ "Mais de 60 anos",
                               TRUE ~ NA_character_),
         catend_cat = case_when(CATEATEND %in% c(1,3) ~ "Convênio/Particular",
                                CATEATEND == 2 ~ "SUS",
                                TRUE ~ NA_character_),
         ecgrup_cat = case_when(ECGRUP %in% c("II","III") ~ "II/III",
                                TRUE ~ ECGRUP),
         m_cat = case_when(M == "0" ~ "Ausência de metástase",
                           TRUE ~ "Presença de metástase"),
         radio_cat = case_when(RADIO == 1 ~ "Com radioterapia",
                               RADIO == 0 ~ "Sem radioterapia"),
         quimio_cat = case_when(QUIMIO == 1 ~ "Com quimioterapia",
                                QUIMIO == 0 ~ "Sem quimoterapia"),
         cirurgia_cat = case_when(CIRURGIA == 1 ~ "Com cirurgia",
                                  CIRURGIA == 0 ~ "Sem cirurgia"),
         rec_cat = case_when(RECNENHUM == 1 ~ "Com recidiva",
                             RECNENHUM == 0 ~ "Sem recidiva"),
         status_cat = case_when(status == 0 ~ "Censura",
                                status == 1 ~ "Morte por câncer",
                                status == 2 ~ "Morte por outras causas"),
         sexo_reg = case_when(SEXO == 1 ~ 0,
                              SEXO == 2 ~ 1,
                              TRUE ~ NA_integer_),
         idade_reg = case_when(IDADE <= 60 ~ 0,
                               IDADE > 60 ~ 1,
                               TRUE ~ NA_integer_),
         catend_reg = case_when(CATEATEND %in% c(1,3) ~ 0,
                                CATEATEND == 2 ~ 1,
                                TRUE ~ NA_integer_),
         ecgrup_2_3_reg = case_when(ECGRUP %in% c("II","III") ~ 1,
                                 TRUE ~ 0),
         ecgrup_4_reg = case_when(ECGRUP == "IV" ~ 1,
                               TRUE ~ 0),
         m_reg = case_when(M == "0" ~ 0,
                           TRUE ~ 1))
table(dados_org$status, dados_org$PERDASEG, useNA="always")
table(dados_org$sexo_cat, dados_org$sexo_reg, useNA="always")
table(dados_org$idade_cat, dados_org$idade_reg, useNA="always")
table(dados_org$catend_cat, dados_org$catend_reg, useNA="always")
table(dados_org$ecgrup_cat, dados_org$ecgrup_2_3_reg, useNA="always")
table(dados_org$ecgrup_cat, dados_org$ecgrup_4_reg, useNA="always")
table(dados_org$m_cat, dados_org$m_reg, useNA="always")

table(dados$IDMITOTIC, useNA="always")
write.csv(dados_org, "dados_ossos.csv", row.names = FALSE)

# Banco em CSV --------
dados_org <- read_csv("dados_ossos.csv")
# Gráfico FIA GERAL ------
setwd("~/Estudos/Doutorado/Matérias/Sobrevivência/Projeto Vertical Promotion/Gráficos")
cores <- c("black", "#d00000","#0349A5", "green" ,"orange","brown","yellow","purple")

fia <- cuminc(
  ftime   = dados_org$time,
  fstatus = dados_org$status,
  cencode = 0)

plot.cuminc(
  fia,
  lty = 1,
  lwd = 2,
  ylim = c(0,0.7),
  cex.axis =1.4,
  cex.lab=1.4,
  cex.main = 1.6,
  cex = 1.4,  
  col = cores, 
  xlab = "Tempo (dias)",
  ylab = "FIA",
  curvlab = c("Morte por câncer",
              "Morte por outras causas"),
  main="Câncer de ossos: CIDs C40 e C41"
)


# Gráfico FIA Variáveis Explicativas----
# Gráfico FIA Sexo
png("fia_graf_cov.png",
    width = 5800,
    height = 4400,
    res = 300)
layout(matrix(c(
  1, 2, 0, 3, 4,
  5, 6, 0, 7, 8,
  9, 10, 0, 11, 12,
  13, 14, 0, 15, 16), nrow = 4, byrow = TRUE),
  widths = c(1, 1, 0.15, 1, 1))

par(mar = c(5, 3, 4, 1),
    oma = c(0, 0.1, 0, 0),
    mgp = c(1.9, 0.5, 0))

fia_sexo <- cuminc(
  ftime   = dados_org$time,
  fstatus = dados_org$status,
  group   = dados_org$sexo_cat,
  cencode = 0)

#Causa 1
plot.cuminc(
  fia_sexo[c(1,2)],
  lty = 1,
  lwd = 2,
  ylim = c(0,0.7),
  cex.axis =1.4,
  cex.lab=1.4,
  cex.main = 1.6,
  cex = 1.4,  
  col = cores, 
  xlab = "Tempo (dias)",
  ylab = "FIA",
  main = "Morte por câncer",
  curvlab = c("Feminino",
              "Masculino"))
text(700, 0.7, "Sexo", font = 2, cex=1.4)
#Causa 2
plot.cuminc(
  fia_sexo[c(3,4)],
  lty = 1,
  lwd = 2,
  ylim = c(0,0.7),
  cex.axis =1.4,
  cex.lab=1.4,
  cex.main = 1.6,
  cex = 1.4,  
  col = cores, 
  xlab = "Tempo (dias)",
  ylab = "FIA",
  curvlab = c("Feminino",
              "Masculino"),
  main = "Morte por outras causas")
text(700, 0.7, "Sexo", font = 2, cex=1.4)

# Gráfico FIA CATEATEND
fia_cat <- cuminc(
  ftime   = dados_org$time,
  fstatus = dados_org$status,
  group   = dados_org$catend_cat,
  cencode = 0)

#Causa 1
plot.cuminc(
  fia_cat[c(1:2)],
  lty = 1,
  lwd = 2,
  ylim = c(0,0.7),
  cex.axis =1.4,
  cex.lab=1.4,
  cex.main = 1.6,
  cex = 1.4,  
  col = cores, 
  xlab = "Tempo (dias)",
  ylab = "FIA",
  curvlab = c("Convênio/Particular",
              "SUS"),
  main = "Morte por câncer")
text(2350, 0.7, "Caráter de atendimento", font = 2, cex=1.4)
#Causa 2
plot.cuminc(
  fia_cat[c(3,4)],
  lty = 1,
  lwd = 2,
  ylim = c(0,0.7),
  cex.axis =1.4,
  cex.lab=1.4,
  cex.main = 1.6,
  cex = 1.4,  
  col = cores, 
  xlab = "Tempo (dias)",
  ylab = "FIA",
  curvlab = c("Convênio/Particular",
              "SUS"),
  main = "Morte por outras causas")
text(2350, 0.7, "Caráter de atendimento", font = 2, cex=1.4)

# Gráfico FIA RADIO 
fia_radio <- cuminc(
  ftime   = dados_org$time,
  fstatus = dados_org$status,
  group   = dados_org$RADIO,
  cencode = 0)

#Causa 1
plot.cuminc(
  fia_radio[c(1:2)],
  lty = 1,
  lwd = 2,
  ylim = c(0,0.85),
  cex.axis =1.4,
  cex.lab=1.4,
  cex.main = 1.6,
  cex = 1.4,  
  wh = c(0, 0.85),
  col = cores, 
  xlab = "Tempo (dias)",
  ylab = "FIA",
  curvlab = c("Sem radioterapia",
              "Com radioterapia"),
  main = "Morte por câncer")
text(1450, 0.85, "Radioterapia", font = 2, cex=1.4)
#Causa 2
plot.cuminc(
  fia_radio[c(3,4)],
  lty = 1,
  lwd = 2,
  ylim = c(0,0.85),
  cex.axis =1.4,
  cex.lab=1.4,
  cex.main = 1.6,
  cex = 1.4,  
  col = cores, 
  xlab = "Tempo (dias)",
  ylab = "FIA",
  curvlab = c("Sem radioterapia",
              "Com radioterapia"),
  main = "Morte por outras causas")
text(1450, 0.85, "Radioterapia", font = 2, cex=1.4)


# Gráfico FIA QUIMIO
fia_quimio <- cuminc(
  ftime   = dados_org$time,
  fstatus = dados_org$status,
  group   = dados_org$QUIMIO,
  cencode = 0)

#Causa 1
plot.cuminc(
  fia_quimio[c(1:2)],
  lty = 1,
  lwd = 2,
  ylim = c(0,0.75),
  cex.axis =1.4,
  cex.lab=1.4,
  cex.main = 1.6,
  cex = 1.4,  
  wh = c(0, 0.75),
  col = cores, 
  xlab = "Tempo (dias)",
  ylab = "FIA",
  curvlab = c("Sem quimioterapia",
              "Com quimioterapia"),
  main = "Morte por câncer")
text(1550, 0.75, "Quimioterapia", font = 2, cex=1.4)
#Causa 2
plot.cuminc(
  fia_quimio[c(3,4)],
  lty = 1,
  lwd = 2,
  ylim = c(0,0.75),
  cex.axis =1.4,
  cex.lab=1.4,
  cex.main = 1.6,
  cex = 1.4,  
  col = cores, 
  xlab = "Tempo (dias)",
  ylab = "FIA",
  curvlab = c("Sem quimioterapia",
              "Com quimioterapia"),
  main = "Morte por outras causas")
text(1550, 0.75, "Quimioterapia", font = 2, cex=1.4)

# Gráfico FIA CIRURGIA  
fia_cirurgia <- cuminc(
  ftime   = dados_org$time,
  fstatus = dados_org$status,
  group   = dados_org$CIRURGIA,
  cencode = 0)

#Causa 1
plot.cuminc(
  fia_cirurgia[c(1:2)],
  lty = 1,
  lwd = 2,
  ylim = c(0,0.85),
  cex.axis =1.4,
  cex.lab=1.4,
  cex.main = 1.6,
  cex = 1.4,  
  wh = c(0, 0.85),
  col = cores, 
  xlab = "Tempo (dias)",
  ylab = "FIA",
  curvlab = c("Sem cirurgia",
              "Com cirurgia"),
  main = "Morte por câncer")
text(1070, 0.85, "Cirurgia", font = 2, cex=1.4)

#Causa 2
plot.cuminc(
  fia_cirurgia[c(3,4)],
  lty = 1,
  lwd = 2,
  ylim = c(0,0.85),
  cex.axis =1.4,
  cex.lab=1.4,
  cex.main = 1.6,
  cex = 1.4,  
  col = cores, 
  xlab = "Tempo (dias)",
  ylab = "FIA",
  curvlab = c("Sem cirurgia",
              "Com cirurgia"),
  main = "Morte por outras causas")
text(1070, 0.85, "Cirurgia", font = 2, cex=1.4)

# Gráfico FIA RECIDIVA  
fia_recidiva <- cuminc(
  ftime   = dados_org$time,
  fstatus = dados_org$status,
  group   = dados_org$RECNENHUM,
  cencode = 0)

#Causa 1
plot.cuminc(
  fia_recidiva[c(1:2)],
  lty = 1,
  lwd = 2,
  ylim = c(0,0.9),
  cex.axis =1.4,
  cex.lab=1.4,
  cex.main = 1.6,
  cex = 1.4,  
  wh = c(0, 0.9),
  col = cores, 
  xlab = "Tempo (dias)",
  ylab = "FIA",
  curvlab = c("Com recidiva",
              "Sem recidiva"),
  main = "Morte por câncer")
text(2150, 0.9, "Recidiva do câncer", font = 2, cex=1.4)
#Causa 2
plot.cuminc(
  fia_recidiva[c(3,4)],
  lty = 1,
  lwd = 2,
  ylim = c(0,0.9),
  cex.axis =1.4,
  cex.lab=1.4,
  cex.main = 1.6,
  cex = 1.4,  
  col = cores, 
  xlab = "Tempo (dias)",
  ylab = "FIA",
  curvlab = c("Com recidiva",
              "Sem recidiva"),
  main = "Morte por outras causas")
text(2150, 0.9, "Recidiva do câncer", font = 2, cex=1.4)

# Gráfico FIA ECGRUP
fia_ec <- cuminc(
  ftime   = dados_org$time,
  fstatus = dados_org$status,
  group   = dados_org$ecgrup_cat,
  cencode = 0)

#Causa 1
plot.cuminc(
  fia_ec[c(1:3)],
  lty = 1,
  lwd = 2,
  ylim = c(0,1),
  cex.axis =1.4,
  cex.lab=1.4,
  cex.main = 1.6,
  cex = 1.4,  
  wh = c(0, 1),
  col = cores, 
  xlab = "Tempo (dias)",
  ylab = "FIA",
  curvlab = c("I",
              "II/III",
              "IV"),
  main = "Morte por câncer")
text(1850, 1, "Estágio tumoral", font = 2, cex=1.4)
#Causa 2
plot.cuminc(
  fia_ec[c(4:6)],
  lty = 1,
  lwd = 2,
  ylim = c(0,1),
  cex.axis =1.4,
  cex.lab=1.4,
  cex.main = 1.6,
  cex = 1.4,  
  col = cores, 
  xlab = "Tempo (dias)",
  ylab = "FIA",
  curvlab = c("I",
              "II/III",
              "IV"),
  main = "Morte por outras causas")
text(1850, 1, "Estágio tumoral", font = 2, cex=1.4)

# Gráfico FIA M 
fia_m <- cuminc(
  ftime   = dados_org$time,
  fstatus = dados_org$status,
  group   = dados_org$m_cat,
  cencode = 0)

#Causa 1
plot.cuminc(
  fia_m[c(1:2)],
  lty = 1,
  lwd = 2,
  ylim = c(0,0.95),
  cex.axis =1.4,
  cex.lab=1.4,
  cex.main = 1.6,
  cex = 1.4, 
  wh = c(0, 0.95),
  col = cores, 
  xlab = "Tempo (dias)",
  ylab = "FIA",
  curvlab = c("Presença de metástase",
              "Ausência de metástase"),
  main = "Morte por câncer")
text(1320, 0.95, "Metástase", font = 2, cex=1.4)
#Causa 2
plot.cuminc(
  fia_m[c(3:4)],
  lty = 1,
  lwd = 2,
  ylim = c(0,0.95),
  cex.axis =1.4,
  cex.lab=1.4,
  cex.main = 1.6,
  cex = 1.4,  
  col = cores, 
  xlab = "Tempo (dias)",
  ylab = "FIA",
  curvlab = c("Presença de metástase",
              "Ausência de metástase"),
  main = "Morte por outras causas")
text(1320, 0.95, "Metástase", font = 2, cex=1.4)

dev.off()

# Descritiva Categórica ----
vars_cat <- c("catend_cat", "ecgrup_cat", "m_cat", "radio_cat",
              "quimio_cat", "cirurgia_cat", "sexo_cat", "rec_cat",
              "status_cat")

nomes_cat <- c("Caráter de atendimento", "Estágio tumoral", "Metástase",
               "Radioterapia", "Quimioterapia", "Cirurgia", "Sexo",
               "Recidiva do câncer", "Status")

# FUNÇÃO PARA OS GRÁFICOS CATEGÓRICOS
grafico_cat <- function(dados, variavel, nome_grafico) {
  dados %>%
    select(all_of(variavel)) %>%
    filter(!is.na(.data[[variavel]])) %>%
    count(.data[[variavel]], name = "n") %>%
    mutate(
      percentual = n / sum(n) * 100,
      categoria = as.character(.data[[variavel]])) %>%
    #Barras
    ggplot(aes(x = reorder(categoria, n), y = n)) +
    geom_col(width = 0.65, fill = "grey70",
             color = "grey25", linewidth = 0.4) +
    #Frequência absoluta e percentual
    geom_text(aes(label = paste0(n," (",round(percentual, 1),"%)")),
              hjust = 1, color = "black",size = 3.4, fontface = "plain") +
    coord_flip() +
    labs(title = nome_grafico, x = NULL, y = "Frequência") +
    theme_classic() +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 12,
        hjust = 0.5,
        color = "black"),
      axis.title = element_text(
        size = 10,
        color = "black"),
      axis.text.y = element_text(
        size = 9,
        color = "black"),
      axis.text.x = element_text(
        size = 8.5,
        color = "black"),
      axis.line = element_line(
        linewidth = 0.5,
        color = "black"),
      axis.ticks = element_line(
        linewidth = 0.4,
        color = "black"),
      panel.grid = element_blank(),
      plot.margin = margin(8, 8, 8, 8)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)))
}

graficos_cat <- Map(
  function(var, nome) {
    grafico_cat(
      dados = dados_org,
      variavel = var,
      nome_grafico = nome) },
  vars_cat,
  nomes_cat)

names(graficos_cat) <- nomes_cat

fig_cat <- wrap_plots(graficos_cat, ncol = 3)
fig_cat

ggsave("descritiva_categoricas.png",
       fig_cat, width = 22, height = 12, 
       units = "in", dpi = 300, bg = "white")

# TABELA DESCRITIVA DAS CATEGÓRICAS
tabela_cat <- dados_org %>%
  select(all_of(vars_cat)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Variavel",
    values_to = "Categoria") %>%
  filter(!is.na(Categoria)) %>%
  count(Variavel, Categoria, name = "n") %>%
  group_by(Variavel) %>%
  mutate(
    percentual = 100 * n / sum(n)) %>%
  ungroup()
tabela_cat

# Descritiva Quantitativa ----
vars_quant <- c("IDADE", "time")
nomes_quant <- c("Idade", "Tempo")

dados_quant <- dados_org %>%
  select(all_of(vars_quant)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Variavel",
    values_to = "Valor") %>%
  mutate(
    Variavel = factor(
      Variavel,
      levels = vars_quant,
      labels = nomes_quant))

# FUNÇÃO PARA OS GRÁFICOS QUANTITATIVOS - Histograma
grafico_hist <- ggplot(dados_quant, aes(x = Valor)) +
  geom_histogram(
    bins = 20,
    fill = "grey70",
    color = "black",
    linewidth = 0.5) +
  facet_wrap(
    ~ Variavel,
    scales = "free",
    ncol = 3) +
  labs(
    x = NULL,
    y = "Frequência") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 12,
      color = "black"),
    axis.title = element_text(
      size = 11,
      color = "black"),
    axis.text = element_text(
      size = 9,
      color = "black"),
    axis.line = element_line(
      linewidth = 0.5),
    panel.grid = element_blank(),
    plot.margin = margin(8, 8, 8, 8))

grafico_hist

# FUNÇÃO PARA OS GRÁFICOS QUANTITATIVOS - Boxplot
grafico_box <- ggplot(dados_quant, aes(x = "", y = Valor)) +
  geom_boxplot(
    width = 0.35,
    fill = "grey70",
    color = "black",
    linewidth = 0.7,
    outlier.shape = 21,
    outlier.size = 2,
    outlier.fill = "white",
    outlier.color = "black") +
  facet_wrap(
    ~ Variavel,
    scales = "free",
    ncol = 3) +
  labs(
    x = NULL,
    y = "Valor") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 12),
    axis.title = element_text(
      size = 11),
    axis.text = element_text(
      size = 9,
      color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank())

grafico_box

# FIGURA FINAL
grafico_quant_completo <-
  grafico_hist +
  grafico_box 

grafico_quant_completo

ggsave(
  "Figura_descritiva_quantitativas.png",
  grafico_quant_completo,
  width = 20,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)


# TABELA DESCRITIVA DAS QUANTITATIVAS
tabela_quant <- dados_org %>%
  select(all_of(vars_quant)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Variavel",
    values_to = "Valor") %>%
  mutate(
    Variavel = factor(
      Variavel,
      levels = vars_quant,
      labels = nomes_quant)) %>%
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
  ggplot(
    aes(x = ano, y = n)) +
  geom_line(
    linewidth = 0.8,
    color = "black") +
  geom_point(
    size = 2.5,
    color = "black") +
  geom_text(
    aes(label = n),
    vjust = -0.7,
    size = 3.2) +
  labs(x = "Ano", y = "Número de pacientes") +
  theme_classic() +
  theme(
    axis.title = element_text(
      size = 11,
      color = "black"),
    axis.text = element_text(
      size = 9,
      color = "black"),
    axis.line = element_line(
      linewidth = 0.5,
      color = "black"),
    panel.grid = element_blank(),
    plot.margin = margin(8, 8, 8, 8)) +
  scale_x_continuous(
    breaks = scales::breaks_pretty(n = 10)) +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.08)))

grafico_ano

grafico_ano_status <- dados_org %>%
  filter(!is.na(ano), !is.na(status_cat)) %>%
  count(ano, status_cat) %>%
  ggplot(aes(x = ano, y = n,
             group = status_cat, color = status_cat)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.3) +
  labs(
    x = "Ano",
    y = "Número de pacientes",
    color = "Status") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    legend.position = c(0.2, 0.90),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.background = element_rect(
      fill = "white",
      color = "grey70",
      linewidth = 0.3),
    legend.key.width = unit(1.2, "cm"),
    panel.grid = element_blank()) +
  scale_color_manual(
    values = c( 
      "Censura" = "black",
      "Morte por câncer" = "#0349A5",
      "Morte por outras causas" = "#d00000")) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 10)) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = expansion(
      mult = c(0, 0)))

grafico_ano_status


grafico_tempo_status <- dados_org %>%
  filter(!is.na(time), !is.na(status_cat)) %>%
  ggplot(aes(x = status_cat, y = time)) +
  geom_boxplot(
    width = 0.5,
    fill = "grey70",
    color = "black",
    linewidth = 0.7,
    outlier.shape = 21,
    outlier.size = 2) +
  labs(
    x = "Status",
    y = "Tempo de acompanhamento") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    panel.grid = element_blank())

grafico_tempo_status

grafico_ano_completo <-
  grafico_ano +  grafico_ano_status + grafico_tempo_status

grafico_ano_completo

ggsave(
  "descritiva_ano.png",
  grafico_ano_completo,
  width = 20,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)

