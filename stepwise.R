# ============================================================
# SELECAO BACKWARD PELO BIC - MODELO VERTICAL COM K = 5
# ============================================================

# Variaveis mantidas obrigatoriamente em cada componente
obrigatorias_cura <- c(
  "IDADE", "ano", "ecgrup_2_3_reg", "ecgrup_4_reg", "m_reg"
)

obrigatorias_causa <- c(
  "log_time", "IDADE", "ano",
  "ecgrup_2_3_reg", "ecgrup_4_reg", "m_reg"
)

# Funcao para calcular o BIC
calcular_bic <- function(fit) {
  p <- length(fit$gamma) +
    length(fit$lambda) +
    length(fit$beta)
  
  -2 * fit$logLik + log(nrow(fit$data)) * p
}

# Funcao para retirar uma covariavel de uma formula
retirar_variavel <- function(formula, variavel) {
  termos <- attr(terms(formula), "term.labels")
  
  reformulate(setdiff(termos, variavel))
}

# Funcao de selecao backward
backward_bic <- function(fit_inicial) {
  
  atual <- fit_inicial
  historico <- data.frame()
  etapa <- 0
  
  repeat {
    
    etapa <- etapa + 1
    bic_atual <- calcular_bic(atual)
    
    # Variaveis candidatas em cada componente
    vars_cura <- setdiff(
      attr(terms(atual$count_formula), "term.labels"),
      obrigatorias_cura
    )
    
    vars_causa <- setdiff(
      attr(terms(atual$cause_formula), "term.labels"),
      obrigatorias_causa
    )
    
    candidatos <- list()
    
    # Testar retirada no componente de cura
    for (v in vars_cura) {
      
      f_cura <- retirar_variavel(atual$count_formula, v)
      
      ajuste <- try(
        fit_vertical_pt_em(
          mel,
          K = 5,
          cuts = fit_k5$cuts,
          count_formula = f_cura,
          cause_formula = atual$cause_formula,
          scaling = fit_k5$scale_info,
          compute_se = FALSE
        ),
        silent = TRUE
      )
      
      if (!inherits(ajuste, "try-error") && ajuste$converged) {
        candidatos[[paste0("Cura:", v)]] <- ajuste
      }
    }
    
    # Testar retirada no componente de causa
    for (v in vars_causa) {
      
      f_causa <- retirar_variavel(atual$cause_formula, v)
      
      ajuste <- try(
        fit_vertical_pt_em(
          mel,
          K = 5,
          cuts = fit_k5$cuts,
          count_formula = atual$count_formula,
          cause_formula = f_causa,
          scaling = fit_k5$scale_info,
          compute_se = FALSE
        ),
        silent = TRUE
      )
      
      if (!inherits(ajuste, "try-error") && ajuste$converged) {
        candidatos[[paste0("Causa:", v)]] <- ajuste
      }
    }
    
    # Se nao houver candidatos, encerrar
    if (length(candidatos) == 0) break
    
    # BIC de cada candidato
    bics <- sapply(candidatos, calcular_bic)
    
    # Melhor candidato da etapa
    melhor <- which.min(bics)
    bic_novo <- bics[melhor]
    
    # Registrar resultados
    historico <- rbind(
      historico,
      data.frame(
        Etapa = etapa,
        Variavel_retirada = names(bics)[melhor],
        BIC_anterior = bic_atual,
        BIC_novo = bic_novo,
        Diferenca_BIC = bic_novo - bic_atual,
        Aceita = bic_novo < bic_atual
      )
    )
    
    cat("\nEtapa:", etapa, "\n")
    cat("Candidata:", names(bics)[melhor], "\n")
    cat("BIC anterior:", round(bic_atual, 4), "\n")
    cat("BIC novo:", round(bic_novo, 4), "\n")
    
    # Parar se nenhuma exclusao melhorar o BIC
    if (bic_novo >= bic_atual - 1e-8) break
    
    # Atualizar modelo
    atual <- candidatos[[melhor]]
  }
  
  list(
    modelo = atual,
    historico = historico
  )
}

# Executar selecao
resultado_backward <- backward_bic(fit_k5)

# Modelo selecionado
fit_backward <- resultado_backward$modelo

# Historico das exclusoes
print(resultado_backward$historico, row.names = FALSE)

# Formulas finais
fit_backward$count_formula
fit_backward$cause_formula
