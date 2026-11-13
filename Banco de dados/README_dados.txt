DADOS DO PROJETO – CÂNCER ÓSSEO E RISCOS COMPETITIVOS
=================================================

1. Origem dos dados

Os dados utilizados neste estudo são provenientes do Registro Hospitalar de Câncer do Estado de São Paulo (RHC/SP), disponibilizado pela Fundação Oncocentro de São Paulo (FOSP), e abrangem o período de 2000 a 2022.

As bases oficiais são disponibilizadas em arquivos no formato .DBF, com registros de diferentes tipos de câncer (CID) em todo o estado de São Paulo. Devido ao tamanho elevado desses arquivos, as bases originais não foram incluídas neste repositório. Caso seja necessário obtê-las, o download pode ser realizado diretamente no site da FOSP:

https://fosp.saude.sp.gov.br/fosp/diretoria-adjunta-de-informacao-e-epidemiologia/rhc-registro-hospitalar-de-cancer/download-de-arquivos/

2. Base tratada disponível no repositório

O arquivo "Banco de Dados/dados_ossos.csv" contém a base filtrada, limpa e organizada utilizada nas análises. Foram selecionados pacientes com neoplasias malignas dos ossos e das cartilagens articulares, identificadas pelos códigos C40 e C41 da CID-10, e mantidas as observações com informações completas para as variáveis selecionadas no estudo. A base analítica final contém 2.613 pacientes.

O período de diagnóstico vai de 2000 a 2022. A idade ao diagnóstico varia de 0 a 90 anos (média de 29 anos e mediana de 22 anos). O tempo de acompanhamento, medido em dias desde o diagnóstico até o evento ou a última informação disponível, varia de 1 a 9.480 dias (média de 2.089 e mediana de 1.219 dias).

Na base tratada, a variável "status" identifica três situações:
- 0: censura;
- 1: morte por câncer ósseo;
- 2: morte por outras causas.

Na amostra, 1.232 pacientes (47,1%) foram censurados, 1.132 (43,3%) apresentaram óbito por câncer ósseo e 249 (9,5%) apresentaram óbito por outras causas. A base também contém informações demográficas, clínicas e assistenciais, como sexo, idade, ano do diagnóstico, localização tumoral, estágio clínico, metástase, tipo de atendimento, cirurgia, quimioterapia, radioterapia, recidiva e perda de seguimento, além de variáveis recodificadas para as análises.

3. Preparação e análise dos dados

A ordem de execução dos programas é:

1) 01_importacao_limpeza.R – realiza o download das bases oficiais da FOSP, a importação dos arquivos .DBF, a seleção dos registros de interesse, a limpeza e a organização das variáveis. Ao final, salva a base tratada "dados_ossos.csv", disponibilizada neste repositório.

2) 02_analise_descritiva.R – utiliza a base já tratada para produzir a análise descritiva, incluindo resumos, tabelas e representações gráficas das características dos pacientes e dos desfechos.

3) 03_modelagem.R – utiliza a base tratada para os ajustes e resultados da modelagem de sobrevivência com riscos competitivos e fração de cura.

Para reproduzir apenas as análises posteriores à preparação dos dados, pode-se utilizar diretamente o arquivo "dados_ossos.csv", sem baixar novamente todas as bases oficiais. Para reproduzir o processamento desde os dados brutos, devem-se obter os arquivos .DBF no site da FOSP e executar primeiro "01_importacao_limpeza.R". Os caminhos de leitura e gravação dos arquivos devem ser ajustados ao ambiente local.

4. Acesso e referência à fonte

O repositório do projeto é privado; é necessário solicitar autorização para acessar os arquivos nele disponibilizados:

https://github.com/CarolyneBrito/Vertical_promotion

Fonte dos dados: Fundação Oncocentro de São Paulo (FOSP), Registro Hospitalar de Câncer do Estado de São Paulo (RHC/SP), registros de 2000 a 2022. Para definições e códigos das variáveis originais, recomenda-se consultar o dicionário de dados disponibilizado pela instituição.
