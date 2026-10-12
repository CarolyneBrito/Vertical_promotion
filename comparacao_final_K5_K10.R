# =============================================================================
# ANALISE FINAL: K=5 (BIC) versus K=10 (AIC)
# Executar APOS definir fit_vertical_pt_em(), bootstrap_full(),
# predict_profile_quantities(), plot_nicolaie_figures() e dados_org.
# =============================================================================

# ------------------------------- Data analysis -------------------------------
mel <- dados_org %>%
  dplyr::select(time, status, sexo_reg, IDADE, ano, catend_reg,
                ecgrup_2_3_reg, ecgrup_4_reg, m_reg,
                RADIO, QUIMIO, CIRURGIA) %>%
  dplyr::mutate(status = dplyr::case_when(
    status == 0 ~ 2L,
    status == 1 ~ 1L,
    status == 2 ~ 3L
  ))

# Modelos piecewise-exponential: K=1,...,10
modelos <- lapply(1:10, function(k) {
  cat('\n-----------------------------------\n')
  cat('Iniciando modelo com K =', k, '\n')
  cat('Horario de inicio:', format(Sys.time(), '%H:%M:%S'), '\n')
  inicio <- Sys.time()
  fit <- fit_vertical_pt_em(mel, K = k, trace = FALSE, compute_se = FALSE)
  tempo <- as.numeric(difftime(Sys.time(), inicio, units = 'secs'))
  cat('Modelo K =', k, 'finalizado!\n')
  cat('Convergiu:', fit$converged, '\n')
  cat('Iteracoes do EM:', fit$iterations, '\n')
  cat('Tempo gasto:', round(tempo, 2), 'segundos\n')
  cat('-----------------------------------\n')
  fit$tempo <- tempo
  fit
})
names(modelos) <- paste0('K=', 1:10)

# Comparacao por log-verossimilhanca, AIC, AICc e BIC
comparacao <- do.call(rbind, lapply(1:10, function(k) {
  fit <- modelos[[k]]
  npar <- length(fit$gamma) + length(fit$lambda) + length(fit$beta)
  n <- nrow(fit$data)
  aic <- -2 * fit$logLik + 2 * npar
  data.frame(
    K = k, logLik = fit$logLik, npar = npar,
    AIC = aic,
    AICc = if (n > npar + 1) aic + 2*npar*(npar+1)/(n-npar-1) else Inf,
    BIC = -2 * fit$logLik + log(n)*npar,
    convergiu = fit$converged,
    tempo_seg = fit$tempo
  )
}))
print(comparacao)
write.csv(comparacao, 'comparacao_modelos_K1_K10.csv', row.names = FALSE)

# Diagnosticos descritivos dos intervalos
print(lapply(modelos, function(fit) fit$d_interval))
print(lapply(modelos, function(fit) round(fit$lambda, 4)))
print(lapply(modelos, function(fit) round(fit$cuts, 2)))
print(data.frame(
  K = 1:length(modelos),
  Cura_media = sapply(modelos, function(fit) mean(fit$cure)),
  Eventos_min = sapply(modelos, function(fit) min(fit$d_interval)),
  Eventos_max = sapply(modelos, function(fit) max(fit$d_interval))
))

# Selecao principal pelo BIC e sensibilidade ao modelo escolhido pelo AIC
melhor_k <- comparacao$K[which.min(ifelse(comparacao$convergiu,
                                          comparacao$BIC, Inf))]
cat('Melhor modelo pelo BIC: K =', melhor_k, '\n')
fit_k5 <- modelos[['K=5']]
fit_k10 <- modelos[['K=10']]
stopifnot(fit_k5$converged, fit_k10$converged)
fit_melanoma <- fit_k5  # preserva o nome anterior, agora modelo principal
print(fit_k5)
print(fit_k10)

# -------------------------- Perfis de pacientes ------------------------------
# Tres cenarios retirados de combinacoes REALMENTE observadas na amostra:
# Estagio I sem metastase; II/III sem metastase; IV com metastase.
# Para cada estrato, seleciona-se um paciente proximo das medianas de idade e ano.
# As demais covariaveis sao preservadas, nao inventadas.
mel_analise <- fit_k5$raw_data  # mesma amostra completa efetivamente usada no ajuste
vars_perfil <- c('sexo_reg', 'IDADE', 'ano', 'catend_reg',
                 'ecgrup_2_3_reg', 'ecgrup_4_reg', 'm_reg',
                 'RADIO', 'QUIMIO', 'CIRURGIA')

escolher_perfil <- function(ec23, ec4, metastase) {
  idx <- which(mel_analise$ecgrup_2_3_reg == ec23 &
               mel_analise$ecgrup_4_reg == ec4 &
               mel_analise$m_reg == metastase)
  if (!length(idx)) stop('Nao ha pacientes na amostra com este perfil clinico.')
  sub <- mel_analise[idx, , drop = FALSE]
  # Distancia padronizada das medianas no proprio estrato
  esc_idade <- max(stats::sd(mel_analise$IDADE), 1)
  esc_ano <- max(stats::sd(mel_analise$ano), 1)
  distancia <- abs(sub$IDADE - stats::median(sub$IDADE))/esc_idade +
    abs(sub$ano - stats::median(sub$ano))/esc_ano
  perfil <- sub[which.min(distancia), vars_perfil, drop = FALSE]
  rownames(perfil) <- NULL
  perfil
}

perfis <- list(
  'Estagio I - sem metastase' = escolher_perfil(0, 0, 0),
  'Estagio II/III - sem metastase' = escolher_perfil(1, 0, 0),
  'Estagio IV - com metastase' = escolher_perfil(0, 1, 1)
)

cat('\nNumero de pacientes nos estratos considerados:\n')
print(c(
  'Estagio I - sem metastase' = sum(mel_analise$ecgrup_2_3_reg == 0 &
                                    mel_analise$ecgrup_4_reg == 0 & mel_analise$m_reg == 0),
  'Estagio II/III - sem metastase' = sum(mel_analise$ecgrup_2_3_reg == 1 &
                                         mel_analise$ecgrup_4_reg == 0 & mel_analise$m_reg == 0),
  'Estagio IV - com metastase' = sum(mel_analise$ecgrup_2_3_reg == 0 &
                                     mel_analise$ecgrup_4_reg == 1 & mel_analise$m_reg == 1)
))
print(do.call(rbind, perfis))

# Mantem o antigo nome do perfil principal para facilitar o reaproveitamento
profile_nicolaie <- perfis[['Estagio II/III - sem metastase']]
perfil_referencia <- profile_nicolaie

# ----------------------- Predicoes por perfil --------------------------------
# Os tempos do modelo estao em ANOS, embora os dados brutos estejam em DIAS.
time_grid <- seq(0, max(mel_analise$time)/365.25, length.out = 301)

pred_perfis_k5 <- lapply(perfis, function(p) {
  predict_profile_quantities(fit_k5, p, time_grid)
})
pred_perfis_k10 <- lapply(perfis, function(p) {
  predict_profile_quantities(fit_k10, p, time_grid)
})

# Curvas de CIF por perfil: 3 perfis x 2 causas, K=5 versus K=10
plot_cif_perfis <- function(arquivo = NULL) {
  if (!is.null(arquivo)) grDevices::pdf(arquivo, width = 11, height = 11)
  op <- graphics::par(no.readonly = TRUE)
  graphics::par(mfrow = c(3, 2), mar = c(4, 4, 3, 1), las = 1)
  vars <- c('cif_melanoma_population', 'cif_other_population')
  titulos <- c('Morte por cancer osseo', 'Morte por outras causas')
  for (nm in names(perfis)) {
    a <- pred_perfis_k5[[nm]]
    b <- pred_perfis_k10[[nm]]
    for (j in seq_along(vars)) {
      v <- vars[j]
      graphics::plot(a$time, a[[v]], type = 'l', col = 'steelblue', lwd = 2,
                     ylim = c(0, 1), xlab = 'Tempo desde o diagnostico (anos)',
                     ylab = 'Incidencia acumulada', main = paste(nm, '-', titulos[j]))
      graphics::lines(b$time, b[[v]], col = 'darkred', lty = 2, lwd = 2)
      graphics::legend('topleft', c('K=5', 'K=10'),
                       col = c('steelblue', 'darkred'),
                       lty = c(1, 2), lwd = 2, bty = 'n')
    }
  }
  graphics::par(op)
  if (!is.null(arquivo)) grDevices::dev.off()
}
plot_cif_perfis('CIF_perfis_K5_K10.pdf')
plot_cif_perfis()

# Comparacao numerica: maxima diferenca absoluta de cada CIF por perfil
resumo_diferencas_perfis <- do.call(rbind, lapply(names(perfis), function(nm) {
  a <- pred_perfis_k5[[nm]]
  b <- pred_perfis_k10[[nm]]
  data.frame(
    perfil = nm,
    cura_k5 = a$cure_probability[1],
    cura_k10 = b$cure_probability[1],
    diferenca_cura = b$cure_probability[1] - a$cure_probability[1],
    max_dif_cif_cancer = max(abs(b$cif_melanoma_population - a$cif_melanoma_population)),
    max_dif_cif_outras = max(abs(b$cif_other_population - a$cif_other_population))
  )
}))
print(resumo_diferencas_perfis)
write.csv(resumo_diferencas_perfis, 'comparacao_perfis_K5_K10.csv', row.names = FALSE)

# ------------------------ CIF populacional geral ------------------------------
# Nao e o perfil medio: media das CIFs individuais em TODOS os pacientes.
# O predict() original integra individualmente e pode demorar varios minutos.
prever_geral <- function(fit, times) {
  cat('\nCalculando CIF geral para K =', length(fit$lambda), '...\n')
  inicio <- Sys.time()
  cif <- predict(fit, times, type = 'cif')
  S <- predict(fit, times, type = 'survival')
  ans <- data.frame(
    time = cif$times,
    survival_population = colMeans(S),
    cif_melanoma_population = colMeans(cif$cause1_melanoma),
    cif_other_population = colMeans(cif$cause3_other)
  )
  cat('Finalizado em', round(as.numeric(difftime(Sys.time(), inicio,
                                             units = 'secs')), 1), 'segundos\n')
  ans
}
pred_geral_k5 <- prever_geral(fit_k5, time_grid)
pred_geral_k10 <- prever_geral(fit_k10, time_grid)

plot_cif_geral <- function(arquivo = NULL) {
  if (!is.null(arquivo)) grDevices::pdf(arquivo, width = 12, height = 4)
  op <- graphics::par(no.readonly = TRUE)
  graphics::par(mfrow = c(1, 3), mar = c(4, 4, 3, 1), las = 1)
  vars <- c('cif_melanoma_population', 'cif_other_population',
            'survival_population')
  titulos <- c('CIF geral - cancer osseo', 'CIF geral - outras causas',
               'Sobrevivencia geral')
  for (j in seq_along(vars)) {
    v <- vars[j]
    graphics::plot(pred_geral_k5$time, pred_geral_k5[[v]], type = 'l',
                   col = 'steelblue', lwd = 2, ylim = c(0, 1),
                   xlab = 'Tempo desde o diagnostico (anos)',
                   ylab = 'Probabilidade', main = titulos[j])
    graphics::lines(pred_geral_k10$time, pred_geral_k10[[v]],
                    col = 'darkred', lwd = 2, lty = 2)
    graphics::legend('topright', c('K=5', 'K=10'),
                     col = c('steelblue', 'darkred'),
                     lty = c(1, 2), lwd = 2, bty = 'n')
  }
  graphics::par(op)
  if (!is.null(arquivo)) grDevices::dev.off()
}
plot_cif_geral('CIF_geral_K5_K10.pdf')
plot_cif_geral()

# Diagnostico de integracao numerica: deve ser proximo de zero.
# A integral e trapezoidal e pode ter pequeno erro.
print(c(
  erro_identidade_k5 = max(abs(with(pred_geral_k5,
     survival_population + cif_melanoma_population + cif_other_population - 1))),
  erro_identidade_k10 = max(abs(with(pred_geral_k10,
     survival_population + cif_melanoma_population + cif_other_population - 1)))
))
print(data.frame(
  K = c(5, 10),
  cura_media = c(mean(fit_k5$cure), mean(fit_k10$cure)),
  max_dif_cif_cancer = rep(max(abs(pred_geral_k5$cif_melanoma_population -
                                    pred_geral_k10$cif_melanoma_population)), 2),
  max_dif_cif_outras = rep(max(abs(pred_geral_k5$cif_other_population -
                                    pred_geral_k10$cif_other_population)), 2)
))

# ---------------------------- Bootstrap -------------------------------------
# Bootstrap completo apenas para o perfil de referencia em K=5 e K=10.
# Nao repetir o bootstrap para cada perfil: pode demorar MUITO.
# Para um piloto: Sys.setenv(B_BOOT = '50') antes de executar esta secao.
B_BOOT <- as.integer(Sys.getenv('B_BOOT', '500'))

cat('\nIniciando bootstrap K=5\n')
boot_k5 <- bootstrap_full(
  fit_k5, profile = perfil_referencia, times = time_grid,
  B = B_BOOT, seed = 2026, trace = TRUE
)
cat('\nIniciando bootstrap K=10\n')
boot_k10 <- bootstrap_full(
  fit_k10, profile = perfil_referencia, times = time_grid,
  B = B_BOOT, seed = 2026, trace = TRUE
)
boot_melanoma <- boot_k5  # preserva nome anterior

print(boot_k5$parameter_table, digits = 4, row.names = FALSE)
print(boot_k10$parameter_table, digits = 4, row.names = FALSE)
cat('Bootstrap K=5:', boot_k5$successful, 'de', B_BOOT, '\n')
cat('Bootstrap K=10:', boot_k10$successful, 'de', B_BOOT, '\n')

# Probabilidade de cura do perfil de referencia, com IC percentil em cada K
cure_k5 <- boot_k5$bands$cure_probability
cure_k10 <- boot_k10$bands$cure_probability
cure_summary <- data.frame(
  K = c(5, 10), perfil = rep('Estagio II/III - sem metastase', 2),
  estimativa = c(cure_k5$estimate[1], cure_k10$estimate[1]),
  lower_95 = c(cure_k5$lower[1], cure_k10$lower[1]),
  upper_95 = c(cure_k5$upper[1], cure_k10$upper[1])
)
print(cure_summary)
write.csv(cure_summary, 'cura_perfil_referencia_K5_K10.csv', row.names = FALSE)

# Graficos originais de quatro paineis para o mesmo perfil, separados por K.
# A legenda do codigo original deve ser revisada para a causa 'cancer osseo'.
plot_nicolaie_figures(fit_k5, boot_k5,
                      file = 'vertical_promotion_K5_referencia.pdf')
plot_nicolaie_figures(fit_k10, boot_k10,
                      file = 'vertical_promotion_K10_referencia.pdf')

# ----------------------------- Diagnosticos ---------------------------------
diagnosticos <- lapply(list('K=5' = fit_k5, 'K=10' = fit_k10), function(fit) {
  list(
    em_converged = fit$converged,
    loglik_monotone = all(diff(fit$logLik_history) >= -1e-7),
    events_per_interval = fit$d_interval,
    events_by_cause = table(factor(fit$data$status[fit$data$status != 2],
                                   levels = c(1, 3)))
  )
})
diagnosticos[['K=5']]$bootstrap_successful <- boot_k5$successful
diagnosticos[['K=5']]$bootstrap_failure_rate <- boot_k5$failure_rate
diagnosticos[['K=10']]$bootstrap_successful <- boot_k10$successful
diagnosticos[['K=10']]$bootstrap_failure_rate <- boot_k10$failure_rate
print(diagnosticos)

# Checagem: componente de causa e identico, pois somente K foi modificado.
cat('Maior diferenca dos beta:', max(abs(fit_k5$beta - fit_k10$beta)), '\n')

# Importante: bandas bootstrap sao para o PERFIL de referencia e
# condicionais aos cortes originais. As CIFs GERAIS acima sao pontuais,
# padronizadas pela distribuicao observada das covariaveis.
