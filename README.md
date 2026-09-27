# Simulação do Canal 3GPP TR 38.901

Implementação Orientada a Objetos em Python de um modelo estocástico de canal de propagação em frequências de 0.5 a 100 GHz, baseado nas diretrizes do relatório técnico 3GPP TR 38.901. O código abrange desde a parametrização de larga escala até a geração de desvanecimento rápido, dispersão Doppler e modelagem em banda básica do sinal recebido.

Este projeto foi desenvolvido para a disciplina de Comunicações Móveis (PPGEE2166) do Departamento de Engenharia Elétrica da Universidade de Brasília (UnB).

## Funcionalidades Simuladas (Cenário UMi - Urban Microcell)
* **Geometria 3D:** Cálculo de distâncias diretas e projetadas com probabilidade de Linha de Visada (LoS / NLoS).
* **Parâmetros de Larga Escala:** Sorteio log-normal de Espalhamento de Atraso (DS), Espalhamento Angular (AoA, ZoA, AoD, ZoD), Sombreamento e Fator de Rice.
* **Pequena Escala (Multi-percurso):** Geração do Perfil de Atraso de Potência (PDP) exponencial e espectros angulares espaciais.
* **Efeito Doppler:** Inserção de mobilidade vetorial 3D e cálculo do desvanecimento rápido no tempo.
* **Seletividade do Canal:** Convolução de pulsos retangulares cruzando larguras de banda com a dispersão do canal (demonstração prática de Flat Fading vs. Frequency Selective Fading).
* **Análise de Coerência:** Função de autocorrelação bidimensional do canal para extração da Banda de Coerência ($B_c$) e Tempo de Coerência ($T_c$).

## Estrutura do Projeto e Módulos

A implementação foi dividida modularmente para separar a geração estatística da extração de gráficos e análises. Abaixo estão os principais arquivos e suas responsabilidades:

* `main.py`
  * **Descrição:** Script orquestrador da simulação. É o ponto de entrada que instancia as classes do canal, define as variáveis de ambiente (alturas, distâncias, velocidade do usuário) e chama as funções de plotagem para gerar os resultados visuais.

* `modelo_canal.py`
  * **Descrição:** Contém o "motor" matemático do 3GPP TR 38.901. Implementa a lógica Orientada a Objetos para o sorteio das variáveis aleatórias, cálculos de probabilidade LoS/NLoS, distribuição de potência nos clusters (PDP), correção numérica dos ângulos de chegada (Azimute/Elevação) e a matriz de vetores 3D de incidência.

* `simulacao_sinal.py`
  * **Descrição:** Responsável por simular o sinal no domínio do tempo. Constrói o sinal recebido em banda básica (soma das componentes com suas respectivas potências, atrasos e fases variantes no tempo) e gera as matrizes de análise cruzando durações de pulso ($\delta t$) com espalhamentos de atraso ($\sigma_\tau$).

* `autocorrelacao.py`
  * **Descrição:** Implementa a equação de Autocorrelação Bidimensional Normalizada do canal ($\rho_{TT}$). Gera as superfícies 3D e extrai as fatias transversais para analisar o impacto do espalhamento de atraso na Banda de Coerência e o impacto da velocidade no Tempo de Coerência.

## Dependências

Certifique-se de instalar as bibliotecas necessárias antes de executar os módulos e utilitários:
```bash
pip install numpy matplotlib python-pptx
