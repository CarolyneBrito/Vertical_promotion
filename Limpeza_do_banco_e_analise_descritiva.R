pacman::p_load(foreign, tidyverse, survival, cmprsk)

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


# Gráfico FIA SEXO ----
png("fia_graf_cov.png",
    width = 5800,
    height = 5500,
    res = 300)
layout(matrix(c(
  1, 2, 0, 3, 4,
  5, 6, 0, 7, 8,
  9, 10, 0, 11, 12,
  13, 14, 0, 15, 16,
  17, 18, 0, 19,20), nrow = 5, byrow = TRUE),
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


# Gráfico FIA IDADE ----
fia_idade <- cuminc(
  ftime   = dados_org$time,
  fstatus = dados_org$status,
  group   = dados_org$idade_cat,
  cencode = 0)

#Causa 1
plot.cuminc(
  fia_idade[c(1:2)],
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
  curvlab = c("Até 60 anos",
              "Mais de 60 anos"),
  main = "Morte por câncer")

#Causa 2
plot.cuminc(
  fia_idade[c(3,4)],
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
  curvlab = c("Até 60 anos",
              "Mais de 60 anos"),
  main = "Morte por outras causas")


# Gráfico FIA CATEATEND ----
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


# Gráfico FIA RADIO ----
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
  ylim = c(0,0.8),
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

#Causa 2
plot.cuminc(
  fia_radio[c(3,4)],
  lty = 1,
  lwd = 2,
  ylim = c(0,0.8),
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

#dev.off()

# Gráfico FIA QUIMIO ----
# png("fia_graf_2.png",
#     width = 5800,
#     height = 3300,
#     res = 300)
# layout(matrix(c(
#   1, 2, 0, 3, 4,
#   5, 6, 0, 7, 8,
#   9, 10, 0, 11, 12), nrow = 3, byrow = TRUE),
#   widths = c(1, 1, 0.15, 1, 1))
# 
# par(mar = c(5, 3, 4, 1),
#     oma = c(0, 0.1, 0, 0),
#     mgp = c(1.9, 0.5, 0))

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

# Gráfico FIA CIRURGIA  ----
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
  ylim = c(0,0.8),
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

#Causa 2
plot.cuminc(
  fia_cirurgia[c(3,4)],
  lty = 1,
  lwd = 2,
  ylim = c(0,0.8),
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


# Gráfico FIA RECIDIVA  ----
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
  ylim = c(0,0.85),
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

#Causa 2
plot.cuminc(
  fia_recidiva[c(3,4)],
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
  curvlab = c("Com recidiva",
              "Sem recidiva"),
  main = "Morte por outras causas")


# Gráfico FIA ECGRUP  ----
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
  ylim = c(0,0.9),
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

#Causa 2
plot.cuminc(
  fia_ec[c(4:6)],
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
  curvlab = c("I",
              "II/III",
              "IV"),
  main = "Morte por outras causas")

# Gráfico FIA M  ----
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
  ylim = c(0,0.9),
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

#Causa 2
plot.cuminc(
  fia_m[c(3:4)],
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
  curvlab = c("Presença de metástase",
              "Ausência de metástase"),
  main = "Morte por outras causas")


dev.off()




