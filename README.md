# Simulação do Canal 3GPP TR 38.901

Implementação em Python de um modelo estocástico de canal de propagação em frequências de 0.5 a 100 GHz, baseado nas diretrizes do relatório técnico 3GPP TR 38.901. O código abrange desde a parametrização de larga escala até a geração de desvanecimento rápido, dispersão Doppler e modelagem em banda básica do sinal recebido.

Este projeto foi desenvolvido para a disciplina de Comunicações Móveis (PPGEE2166) do Departamento de Engenharia Elétrica da Universidade de Brasília (UnB).

## Funcionalidades Simuladas (Cenário UMi - Urban Microcell)
* **Geometria 3D:** Cálculo de distâncias diretas e projetadas com probabilidade de Linha de Visada (LoS / NLoS).
* **Parâmetros de Larga Escala:** Sorteio estatístico de Espalhamento de Atraso (DS), Espalhamento Angular (AoA, ZoA, AoD, ZoD), Sombreamento e Fator de Rice.
* **Pequena Escala (Multipercurso):** Geração do Perfil de Atraso de Potência (PDP) e espectros angulares espaciais.
* **Efeito Doppler:** Inserção de mobilidade vetorial 3D e cálculo do desvanecimento rápido no tempo.
* **Seletividade do Canal:** Convolução de pulsos retangulares cruzando larguras de banda com a dispersão do canal.
* **Análise de Coerência:** Função de autocorrelação bidimensional do canal para extração da Banda de Coerência ($B_c$) e Tempo de Coerência ($T_c$).

## Estrutura do Projeto e Módulos

A implementação está dividida de forma modular. Na raiz do projeto encontram-se os guiões principais de execução, enquanto o subdiretório `Classes/` contém os módulos matemáticos específicos da norma.

### Ficheiros Principais (Raiz)
* `path_loss_3GPP_TR_38_901.py`: Principal classe da simulação. Integra as classes, calcula as perdas de percurso (path loss) e define o cenário base.
* `potencia_recebida.py`: Responsável por simular o sinal no domínio do tempo, reconstruindo o sinal recebido em banda básica e demonstrando a seletividade em frequência.
* `coerencia_canal.py`: Implementa a análise de autocorrelação bidimensional ($\rho_{TT}$), gerando as superfícies 3D e extraindo o impacto do espalhamento de atraso e da velocidade no canal.

### Módulos Auxiliares (Diretório `Classes/`)
* `gerar_param.py`: Define e inicializa os parâmetros geométricos e as variáveis globais do cenário (alturas, distâncias e probabilidades de visada).
* `gerar_param_larga_escala.py`: Realiza o sorteio log-normal e a atribuição dos parâmetros de larga escala (DS, espalhamentos angulares, fator de Rice e shadowing).
* `gerar_atrasos.py`: Focado na modelação de pequena escala, gera o Perfil de Atraso de Potência (PDP) exponencial para as componentes multipercurso.
* `gerar_direcoes_chegadas.py`: Executa o cálculo complexo das direções de chegada (Azimute e Elevação) baseando-se nas potências, e constrói a matriz de vetores unitários 3D.
* `calc_desv_doopler_fase.py`: Processa a interação entre a mobilidade do utilizador (vetor velocidade) e as componentes 3D para calcular o desvio Doppler e as respetivas flutuações de fase no tempo.

## Dependências

Certifique-se de instalar as bibliotecas necessárias antes de executar os módulos:
```bash
pip install numpy matplotlib
