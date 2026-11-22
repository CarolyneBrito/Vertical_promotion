# Modelo vertical de tempo de promoção com riscos competitivos e fração de cura

Este repositório possui os dados tratados, os arquivos em R e os resultados da pesquisa **A Promotion-Time Cure Model for Competing Risks under Vertical Modeling: An Application to Bone Neoplasm Data**. 
O trabalho propõe a combinação de um modelo de cura por tempo de promoção com a modelagem vertical de riscos competitivos e apresenta uma aplicação a pacientes com neoplasias malignas dos ossos e das cartilagens articulares (CID-10 C40 e C41) registrados no estado de São Paulo, no período de 2000 a 2022.

## Dados e fonte oficial

Os dados originais são provenientes do Registro Hospitalar de Câncer do Estado de São Paulo (RHC/SP), disponibilizado pela Fundação Oncocentro de São Paulo (FOSP):
https://fosp.saude.sp.gov.br/fosp/diretoria-adjunta-de-informacao-e-epidemiologia/rhc-registro-hospitalar-de-cancer/download-de-arquivos/

As bases oficiais, originalmente em formato `.DBF`, reúnem informações de todas as categorias de câncer registradas nas instituições participantes do estado de São Paulo. Esses arquivos originais não foram incluídos no GitHub devido ao seu tamanho elevado. É possível obtê-los no endereço oficial acima. O script `01_importacao_limpeza.R` contém os códigos para importação, filtragem, limpeza e organização dessas bases.
A base final utilizada nas análises, com 2.613 pacientes, está disponível neste repositório em: `Banco de Dados/dados_ossos.csv`
O arquivo contém os dados filtrados e tratados para as análises. Informações complementares sobre sua origem, variáveis e organização estão no arquivo `README_dados.txt`, na mesma pasta.

## Organização dos arquivos
```text
Vertical_promotion/
├── README.md
├── Artigo/
│   ├── Artigos base/
│   ├── Figuras/
│   ├── artigo.pdf
│   ├── artigo.tex
│   └── referencias.bib
├── Banco de Dados/
|   ├── DICIONARIO_DADOS_RHC-SP_FOSP
│   ├── README_dados.txt
│   └── dados_ossos.csv
├── Codigo/
│   ├── 01_importacao_limpeza.R
│   ├── 02_analise_descritiva.R
│   └── 03_modelagem.R
└── Resultados_codigo/
    ├── tabelas de resultados (.csv)
    ├── objetos e resultados da aplicação (.rds)
    └── figuras da modelagem na aplicação (.pdf)
```
Os itens descritivos de `Resultados_codigo/` representam grupos de arquivos, não nomes literais de arquivos.

Na pasta `Resultados_codigo` estão organizados os produtos gerados pela aplicação: tabelas com resultados em formato `.csv`, arquivos `.rds` contendo objetos e resultados produzidos no R e as três figuras da aplicação.

## Ordem de execução dos programas

Os scripts foram escritos em R e devem ser executados na seguinte ordem:

| Ordem | Programa | Descrição |
|:--:|---|---|
| **1** | `01_importacao_limpeza.R` | Download dos arquivos oficiais da FOSP, importação dos dados `.DBF`, seleção dos registros de câncer ósseo (C40/C41), limpeza e organização da base; ao final, gera e salva `dados_ossos.csv`. |
| **2** | `02_analise_descritiva.R` | Análise descritiva e exploratória da base tratada, incluindo tabelas, figuras e análises não paramétricas de sobrevivência e riscos competitivos. |
| **3** | `03_modelagem.R` | Ajuste do modelo vertical de tempo de promoção com fração de cura, comparação de especificações, estimação dos parâmetros e erros-padrão, bootstrap, predições e resultados da aplicação. |


## Ambiente computacional e versões

As análises foram realizadas em R 4.6.1 (2026-06-24 ucrt), em Windows 11 x64, na plataforma `x86_64-w64-mingw32/x64`. 
As versões abaixo correspondem ao ambiente informado por `sessionInfo()` e `packageVersion()`:

| Pacote | Versão |
|---|---|
| `boot` | 1.3.32 |
| `cmprsk` | 2.2-12 |
| `dplyr` | 1.2.1 |
| `foreign` | 0.8-91 |
| `ggplot2` | 4.0.3 |
| `lubridate` | 1.9.5 |
| `patchwork` | 1.3.2 |
| `scales` | 1.4.0 |
| `survival` | 3.8-6 |
| `tidyverse` | 2.0.0 |
| `xtable` | 1.8-8 |
| `pacman` | 0.5.1 |


## Disponibilidade e acesso

Repositório: https://github.com/CarolyneBrito/Vertical_promotion
O repositório é privado e requer solicitação de autorização de acesso. Os arquivos `.DBF` originais permanecem disponíveis para download diretamente no portal público da FOSP.
