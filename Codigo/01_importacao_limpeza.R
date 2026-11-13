pacman::p_load(foreign, tidyverse, survival, cmprsk, patchwork, scales)
setwd("~/Estudos/Doutorado/Materias/Sobrevivência/Projeto Vertical Promotion")

# Salvamento do banco em RDS ------
dados <- lapply(2000:2022, function(ano) {
  cat("Lendo ano:", ano, "\n") 
  read.dbf(paste0("Banco de dados", "/RHC_SP_", ano, ".dbf"), as.is = TRUE) %>% 
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
         RECNENHUM, PERDASEG, ECGRUP, M, ano) %>% 
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

# Critérios de exclusão -----
dados_ex<- dados %>% mutate(time = as.numeric(DTULTINFO - dados$DTDIAG),
                            status = case_when(ULTINFO %in% c(1,2) ~ 0,
                                               ULTINFO == 3 ~ 1,
                                               ULTINFO == 4 ~ 2)) %>% 
  filter(time >= 0) %>% 
  mutate(time = time+1) %>% 
  filter(TOPOGRUP %in% c("C40","C41")) %>% 
  select(time, status, 
         IDADE, SEXO, TOPOGRUP, CATEATEND, RADIO, QUIMIO, CIRURGIA,
         RECNENHUM, PERDASEG, ECGRUP, M, ano)

for (v in names(dados_ex)) {
  x <- dados_ex[[v]]
  if (is.character(x) ||
      is.factor(x) ||
      is.logical(x) ||
      length(unique(x[!is.na(x)])) <= 20) {
    
    cat("\n--------------------------------------------\n")
    cat("VARIÁVEL:", v, "\n")
    
    print(
      sort(
        table(x, useNA = "ifany"),
        decreasing = TRUE
      )
    )
  }
}

# Banco em CSV --------
write.csv(dados_org, "dados_ossos.csv", row.names = FALSE)
dados_org <- read_csv("Banco de dados/dados_ossos.csv")